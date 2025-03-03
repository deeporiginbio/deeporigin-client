"""this module contains various functions to run steps of an ABFE workflow"""

import ast
import os
from dataclasses import dataclass, field
from typing import Any, Callable, List, Optional

import numpy as np
import pandas as pd
from beartype import beartype
from box import Box
from deeporigin import chemistry
from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools import run
from deeporigin.tools.fep import Ligand, Protein, _load_params
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from deeporigin.tools.utils import query_run_status
from deeporigin.utils.config import construct_resource_url
from deeporigin.utils.notebook import render_mermaid
from IPython.display import HTML, display

# constants and types
ABFE_DB = "ABFE"
COL_DELTA_G = "FEP ΔG (kcal/mol)"
COL_LIGAND_FILE = "ligand_file"
COL_PROTEIN_FILE = "protein_file"
COL_COMPLEX_PREP_OUTPUT = "complex_prep_output"
COL_LIGAND_PREP_OUTPUT = "ligand_prep_output"
COL_EMEQ_OUTPUT = "emeq_output"
COL_SOLVATION_FEP_OUTPUT = "solvation_fep_output"
COL_MD_OUTPUT = "md_output"
COL_BINDING_FEP_OUTPUT = "binding_fep_output"
COL_END_TO_END_OUTPUT = "end_to_end_output"
COL_JOB_IDS = "job_ids"
COL_PROTEIN_NAME = "protein_name"
COL_FEP_RESULTS = "fep_results"


ABFE_DIR = os.path.join(os.path.expanduser("~"), ".deeporigin", "abfe")
os.makedirs(ABFE_DIR, exist_ok=True)


@beartype
def _load_all_params() -> Box:
    """utility function to load all parameters"""

    params = Box()
    params.binding_fep = _load_params("binding_fep")
    params.complex_prep = _load_params("complex_prep")
    params.emeq = _load_params("emeq")
    params.ligand_prep = _load_params("ligand_prep")
    params.simple_md = _load_params("simple_md")
    params.solvation_fep = _load_params("solvation_fep")
    params.end_to_end = _load_params("abfe_end_to_end")
    return params


@beartype
def _ensure_db_for_abfe() -> dict:
    """ensure that there is a database for FEP on Data Hub"""

    database = _ensure_database(ABFE_DB)

    try:
        api.add_smiles_column(
            name="Ligand",
            database_id=ABFE_DB,
        )
    except Exception:
        pass

    required_columns = [
        dict(name=COL_DELTA_G, type="float"),
        dict(name=COL_LIGAND_FILE, type="file"),
        dict(name=COL_PROTEIN_NAME, type="text"),
        dict(name=COL_PROTEIN_FILE, type="file"),
        dict(name=COL_COMPLEX_PREP_OUTPUT, type="file"),
        dict(name=COL_LIGAND_PREP_OUTPUT, type="file"),
        dict(name=COL_EMEQ_OUTPUT, type="file"),
        dict(name=COL_SOLVATION_FEP_OUTPUT, type="file"),
        dict(name=COL_MD_OUTPUT, type="file"),
        dict(name=COL_BINDING_FEP_OUTPUT, type="file"),
        dict(name=COL_END_TO_END_OUTPUT, type="file"),
        dict(name=COL_FEP_RESULTS, type="file"),
        dict(name=COL_JOB_IDS, type="text"),
    ]

    database = _ensure_columns(
        database=database,
        required_columns=required_columns,
    )

    return database


@dataclass
class ABFE:
    """ABFE class that can work with one protein and many ligands"""

    protein: Protein
    ligands: List[Ligand]

    delta_gs: List[float] = field(default_factory=list)
    row_ids: List[str] = field(default_factory=list)

    params: dict = field(default_factory=lambda: _load_all_params())

    database: Optional[Box] = None

    df: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(
            columns=[
                "ligand_file",
                "jobID",
                "step",
                "Status",
            ]
        )
    )

    results_df: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(
            columns=[
                "Ligand",
                "SMILES",
                "FEP ΔG (kcal/mol)",
                "Status",
            ]
        )
    )

    @classmethod
    def from_sessions(cls, sessions: List[int]) -> "ABFE":
        """initialize an ABFE class given a list of sessions

        Args:
            sessions (List[int]): list of sessions"""

        df = api.get_dataframe(ABFE_DB, use_file_names=False)
        idx = [ABFE_DB + "-" + str(session) for session in sessions]

        all_valid = set(idx).issubset(df.index)
        if not all_valid:
            raise DeepOriginException("Some sessions are not valid.")

        df = df.loc[idx]

        # check that rows correspond to only one protein
        if len(df["protein_name"].unique()) > 1:
            raise DeepOriginException(
                "List of sessions corresponds to more than 1 protein. Cannot continue"
            )

        protein_name = df["protein_name"].unique()[0]

        # make ligands
        file_ids = list(df["protein_file"]) + list(df["ligand_file"])
        files = api.list_files(file_ids=file_ids)

        ligands = []
        for row_id, row in df.iterrows():
            file_id = row["ligand_file"]
            file_name = [file.file.name for file in files if file.file.id == file_id][0]
            smiles_string = row["Ligand"]
            ligand = Ligand(file=file_name, smiles_string=smiles_string, n_molecules=1)
            ligands.append(ligand)

        files = api.list_files(file_ids=file_ids)

        # make protein
        protein_file_name = [file for file in files if file.file.name.endswith("pdb")][
            0
        ].file.name
        protein = Protein(
            file=protein_file_name,
            name=protein_name,
        )

        # Create the ABFE instance
        abfe = cls(
            ligands=ligands,
            protein=protein,
            delta_gs=[np.nan for _ in ligands],
            database=_ensure_db_for_abfe(),
        )

        abfe.row_ids = df.index

        # now make the status df
        rows = []

        # every ligand is initialized
        for lig in ligands:
            rows.append(
                {
                    "ligand_file": str(lig.file),
                    "jobID": None,
                    "step": "init",
                    "Status": "Succeeded",
                }
            )

        # now put the job IDs
        for ligand, (_, row) in zip(ligands, df.iterrows()):
            if pd.isna(row["job_ids"]):
                continue

            data = ast.literal_eval(row["job_ids"])
            for key in data.keys():
                if data[key] is None:
                    continue

                rows.append(
                    {
                        "ligand_file": str(ligand.file),
                        "jobID": data[key],
                        "step": key,
                        "Status": "Running",
                    }
                )

        abfe.df = pd.DataFrame(rows, columns=abfe.df.columns)

        return abfe

    @classmethod
    def from_dir(cls, directory: str) -> "ABFE":
        """initialize an ABFE class given some files in a directory

        Args:
            directory (str): directory containing ligand and protein files"""

        sdf_files = sorted(
            [
                os.path.join(directory, f)
                for f in os.listdir(directory)
                if f.lower().endswith(".sdf")
            ]
        )
        ligands = [Ligand(sdf_file) for sdf_file in sdf_files]

        pdb_files = [
            os.path.join(directory, f)
            for f in os.listdir(directory)
            if f.lower().endswith(".pdb")
        ]

        if len(pdb_files) != 1:
            raise ValueError(
                f"Expected exactly one PDB file in the directory, but found {len(pdb_files)}."
            )
        protein_file = pdb_files[0]
        protein = Protein(protein_file)

        # Create the ABFE instance
        abfe = cls(
            ligands=ligands,
            protein=protein,
            delta_gs=[np.nan for _ in ligands],
        )

        # Populate a row in abfe.df for each ligand
        rows = []
        for lig in ligands:
            rows.append(
                {
                    "ligand_file": str(lig.file),
                    "jobID": None,
                    "step": None,
                    "Status": "Ready to Upload",
                }
            )

        # Overwrite (or extend) the internal DataFrame
        abfe.df = pd.DataFrame(rows, columns=abfe.df.columns)

        return abfe

    def update(self):
        """update statuses of jobs in ABFE run"""

        for _, row in self.df.iterrows():
            if row["Status"] in ["Succeeded", "Failed"]:
                continue

            if row["jobID"] is None or len(row["jobID"]) == 0:
                continue

            row["Status"] = query_run_status(row["jobID"])

    def init(self):
        """Initialize an ABFE run. Upload ligand(s) and protein files to Deep Origin."""

        self.database = _ensure_db_for_abfe()

        # get params in
        self.params = _load_all_params()

        smiles_column_id = [
            col.id for col in self.database.cols if col.name == "Ligand"
        ][0]
        protein_name_column_id = [
            col.id for col in self.database.cols if col.name == COL_PROTEIN_NAME
        ][0]

        url = construct_resource_url(
            name=ABFE_DB,
            row_type="database",
        )

        print(f"Using database at: {url}")
        print("🧬 Uploading files to database...")

        self.row_ids = []

        # now upload the protein
        protein_file = api.upload_file(file_path=str(self.protein.file))

        self.delta_gs = [np.nan for _ in self.ligands]

        for ligand in self.ligands:
            # for each ligand, upload the ligand to a new row

            ligand_file = ligand.file

            response = api.upload_file_to_new_database_row(
                database_id=ABFE_DB,
                column_id=COL_LIGAND_FILE,
                file_path=ligand.file,
            )

            row_id = response.rows[0].hid

            self.row_ids.append(row_id)

            api.set_cell_data(
                ligand.smiles_string,
                database_id=ABFE_DB,
                row_id=row_id,
                column_id=smiles_column_id,
            )

            api.set_cell_data(
                self.protein.name,
                database_id=ABFE_DB,
                row_id=row_id,
                column_id=protein_name_column_id,
            )

            api.assign_files_to_cell(
                file_ids=[protein_file.id],
                database_id=ABFE_DB,
                column_id=COL_PROTEIN_FILE,
                row_id=row_id,
            )

            # update the dataframe
            self.df.loc[len(self.df)] = [
                ligand_file,
                "",
                "init",
                "Succeeded",
            ]

        print(f"🧬 Files uploaded to row {row_id}.")

    def get_results(self, image_size=(200, 100)):
        """
        Create and display a DataFrame in a Jupyter notebook with columns:
          - Molecule (an inline RDKit image from the SMILES)
          - SMILES (the raw SMILES string)
          - FEP ΔG (kcal/mol)
          - status
        """

        self.update()
        self.get_delta_gs()

        df = self.df.copy()

        # 1. Mark rows where jobID is non-empty
        df["has_jobID"] = df["jobID"].notna() & (df["jobID"] != "")

        # 2. Sort by 'has_jobID' in descending order,
        #    ensuring rows with a valid jobID appear first
        df = df.sort_values("has_jobID", ascending=False)

        # 3. Drop duplicates on 'ligand_file', keeping the first row for each group

        step_order = pd.CategoricalDtype(
            categories=["init", "end_to_end"], ordered=True
        )

        status_order = pd.CategoricalDtype(
            categories=["Queued", "Failed", "Running", "Succeeded"], ordered=True
        )
        df["step"] = df["step"].astype(step_order)
        df["Status"] = df["Status"].astype(status_order)

        df = df.groupby(["ligand_file", "step"], as_index=False, observed=True).agg(
            {"Status": "max"}
        )

        def pick_dominant_row(group):
            # If there's any end_to_end row for this ligand_file, keep only those
            if (group["step"] == "end_to_end").any():
                group = group[group["step"] == "end_to_end"]
            # Among the remaining rows, pick the single row with the highest Status
            # Since Status is an ordered categorical, 'idxmax()' gives the index of the best status
            best_idx = group["Status"].idxmax()
            return group.loc[best_idx]

        df = (
            df.groupby("ligand_file", group_keys=False)
            .apply(pick_dominant_row)
            .reset_index(drop=True)
        )

        # add the SMILES strings
        df["SMILES"] = [ligand.smiles_string for ligand in self.ligands]

        # only show basenames of files
        df["ligand_file"] = [os.path.basename(file) for file in df["ligand_file"]]

        # now add the delta Gs
        df[COL_DELTA_G] = self.delta_gs

        self.results_df = df
        return df

    def get_delta_gs(self):
        """return the delta Gs of the ABFE run and stored them in self.delta_gs"""

        # some early exits
        if len(self.row_ids) == 0:
            return

        # get rows of the dataframe for these runs
        df = api.get_dataframe("ABFE", use_file_names=False)
        df = df.loc[self.row_ids]

        file_ids = list(df["fep_results"])

        file_ids = [x for x in file_ids if pd.notna(x)]

        if len(file_ids) == 0:
            return

        api.download_files(
            file_ids=file_ids,
            save_to_dir=ABFE_DIR,
        )

        delta_gs = []

        for i, row in df.loc[self.row_ids].iterrows():
            file_id = row["fep_results"]
            if pd.isna(file_id):
                delta_gs.append(np.nan)
                continue

            file_id = file_id.replace("_file:", "")

            delta_g = float(
                pd.read_csv(os.path.join(ABFE_DIR, file_id))["Total"].iloc[0]
            )
            delta_gs.append(delta_g)

        self.delta_gs = delta_gs

    def show_results(self):
        """print the results of an ABFE run"""

        self.get_results()

        df = self.results_df.copy()

        # convert SMILES to aligned images
        smiles_list = list(df["SMILES"])
        df.drop("SMILES", axis=1, inplace=True)

        df["Ligands"] = chemistry.smiles_list_to_base64_png_list(smiles_list)

        # Use escape=False to allow the <img> tags to render as images
        display(HTML(df.to_html(escape=False)))

    def end_to_end(self):
        """run an end-to-end job"""

        for ligand, row_id in zip(self.ligands, self.row_ids):
            self.df = _run_job(
                ligand_file=ligand.file,
                row_id=row_id,
                df=self.df,
                params=self.params.end_to_end,
                func=end_to_end,
                database=self.database,
            )

    def status(self):
        """print the status of an ABFE run"""

        if len(self.df["ligand_file"].unique()) != 1:
            raise NotImplementedError(
                "Expected exactly one ligand file in the dataframe. I haven't coded this yet"
            )

        self.update()

        node_statuses = {
            "init": "NotStarted",
            "complex_prep": "NotStarted",
            "ligand_prep": "NotStarted",
            "solvation_FEP": "NotStarted",
            "simple_MD": "NotStarted",
            "binding_FEP": "NotStarted",
            "DeltaG": "NotStarted",
        }

        for _, row in self.df.iterrows():
            node_statuses[row["step"]] = row["Status"]

        render_mermaid_with_statuses(node_statuses)

    def emeq(self):
        """run complex prep step on each ligand"""

        for ligand, row_id in zip(self.ligands, self.row_ids):
            self.df = _run_job(
                ligand_file=ligand.file,
                row_id=row_id,
                df=self.df,
                params=self.params.emeq,
                func=emeq,
                database=self.database,
            )

    def complex_prep(self):
        """run complex prep step on each ligand"""

        for ligand, row_id in zip(self.ligands, self.row_ids):
            self.df = _run_job(
                ligand_file=ligand.file,
                row_id=row_id,
                df=self.df,
                params=self.params.complex_prep,
                func=complex_prep,
                database=self.database,
            )

    def ligand_prep(self):
        """run complex prep step on each ligand"""

        for ligand, row_id in zip(self.ligands, self.row_ids):
            self.df = _run_job(
                ligand_file=ligand.file,
                row_id=row_id,
                df=self.df,
                params=self.params.ligand_prep,
                func=ligand_prep,
                database=self.database,
            )

    def simple_md(self):
        """run complex prep step on each ligand"""

        for ligand, row_id in zip(self.ligands, self.row_ids):
            self.df = _run_job(
                ligand_file=ligand.file,
                row_id=row_id,
                df=self.df,
                params=self.params.simple_md,
                func=simple_md,
                database=self.database,
            )

    def solvation_fep(self):
        """run complex prep step on each ligand"""

        for ligand, row_id in zip(self.ligands, self.row_ids):
            self.df = _run_job(
                ligand_file=ligand.file,
                row_id=row_id,
                df=self.df,
                params=self.params.solvation_fep,
                func=solvation_fep,
                database=self.database,
            )

    def binding_fep(self):
        """run complex prep step on each ligand"""

        for ligand, row_id in zip(self.ligands, self.row_ids):
            self.df = _run_job(
                ligand_file=ligand.file,
                row_id=row_id,
                df=self.df,
                params=self.params.binding_fep,
                func=binding_fep,
                database=self.database,
            )


def end_to_end(
    *,
    row_id: str,
    params: Optional[dict] = None,
    database: Optional[dict] = None,
) -> str:
    """Function to run an end-to-end job"""

    tool_key = "deeporigin.abfe-end-to-end"

    if database is None:
        database = _ensure_database(ABFE_DB)

    if params is None:
        params = _load_params("abfe_end_to_end")

    params["protein"] = {
        "columnId": COL_PROTEIN_FILE,
        "rowId": row_id,
        "databaseId": database.hid,
    }

    params["ligand"] = {
        "columnId": COL_LIGAND_FILE,
        "rowId": row_id,
        "databaseId": database.hid,
    }

    outputs = {
        "output_file": {
            "columnId": COL_END_TO_END_OUTPUT,
            "rowId": row_id,
            "databaseId": database.hid,
        },
        "abfe_results_summary": {
            "columnId": COL_FEP_RESULTS,
            "rowId": row_id,
            "databaseId": database.hid,
        },
    }

    return run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )


@beartype
def emeq(
    *,
    row_id: str,
    params: Optional[dict] = None,
    database: Optional[dict] = None,
) -> str:
    """Run emeq on a ligand and protein pair, that exist as files on a row in the ABFE database. For this to work, the complex prep step must have been run first.

    Args:
        row_id (str): row id that contains the ligand and protein files.

    """

    tool_key = "deeporigin.abfe-emeq"

    if database is None:
        database = _ensure_database(ABFE_DB)

    if params is None:
        params = _load_params("emeq")

    params["input"] = {
        "columnId": COL_COMPLEX_PREP_OUTPUT,
        "rowId": row_id,
        "databaseId": database.hid,
    }

    outputs = {
        "output_file": {
            "columnId": COL_EMEQ_OUTPUT,
            "rowId": row_id,
            "databaseId": database.hid,
        }
    }

    return run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )


@beartype
def ligand_prep(
    *,
    row_id: str,
    params: Optional[dict] = None,
    database: Optional[dict] = None,
) -> str:
    """Function to prepare uploaded Ligand and protein files using Deep Origin MDSuite. Use this function to run system prep on a ligand and protein pair, that exist as files on a row in the ABFE database.

    Args:
        row_id (str): row id that contains the ligand and protein files.


    """

    tool_key = "deeporigin.abfe-prep"

    if database is None:
        database = _ensure_database(ABFE_DB)

    if params is None:
        params = _load_params("ligand_prep")

    params["ligand"] = {
        "columnId": COL_LIGAND_FILE,
        "rowId": row_id,
        "databaseId": database.hid,
    }

    outputs = {
        "output_file": {
            "columnId": COL_LIGAND_PREP_OUTPUT,
            "rowId": row_id,
            "databaseId": database.hid,
        }
    }

    return run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )


@beartype
def complex_prep(
    *,
    row_id: str,
    params: Optional[dict] = None,
    database: Optional[dict] = None,
) -> str:
    """Function to prepare uploaded Ligand and protein files using Deep Origin MDSuite. Use this function to run system prep on a ligand and protein pair, that exist as files on a row in the ABFE database.

    Args:
        row_id (str): row id that contains the ligand and protein files.


    """

    tool_key = "deeporigin.abfe-prep"

    if database is None:
        database = _ensure_database(ABFE_DB)

    if params is None:
        params = _load_params("complex_prep")

    params["ligand"] = {
        "columnId": COL_LIGAND_FILE,
        "rowId": row_id,
        "databaseId": database.hid,
    }
    params["protein"] = {
        "columnId": COL_PROTEIN_FILE,
        "rowId": row_id,
        "databaseId": database.hid,
    }

    outputs = {
        "output_file": {
            "columnId": COL_COMPLEX_PREP_OUTPUT,
            "rowId": row_id,
            "databaseId": database.hid,
        }
    }

    job_id = run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )

    return job_id


@beartype
def solvation_fep(
    *,
    row_id: str,
    params: Optional[dict] = None,
    database: Optional[dict] = None,
) -> str:
    """Run a solvation simulation


    Args:
        row_id (str): row id of the ligand and protein files.
    """

    tool_key = "deeporigin.abfe-solvation-fep"

    if database is None:
        database = _ensure_database(ABFE_DB)

    if params is None:
        params = _load_params("solvation_fep")

    params["input"] = {
        "columnId": COL_LIGAND_PREP_OUTPUT,
        "rowId": row_id,
        "databaseId": ABFE_DB,
    }

    outputs = {
        "output_file": {
            "columnId": COL_SOLVATION_FEP_OUTPUT,
            "rowId": row_id,
            "databaseId": database.hid,
        }
    }

    return run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )


@beartype
def simple_md(
    *,
    row_id: str,
    params: Optional[dict] = None,
    database: Optional[dict] = None,
) -> str:
    """Run a simple MD simulation

    Args:
        row_id (str): row id of the ligand and protein files.
    """

    tool_key = "deeporigin.abfe-simple-md"

    if database is None:
        database = _ensure_database(ABFE_DB)

    if params is None:
        params = _load_params("simple_md")

    params["input"] = {
        "columnId": COL_EMEQ_OUTPUT,
        "rowId": row_id,
        "databaseId": ABFE_DB,
    }

    outputs = {
        "output_file": {
            "columnId": COL_MD_OUTPUT,
            "rowId": row_id,
            "databaseId": database.hid,
        }
    }

    return run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )


@beartype
def binding_fep(
    *,
    row_id: str,
    params: Optional[dict] = None,
    database: Optional[dict] = None,
) -> str:
    """Run an ABFE simulation

    Args:
        row_id (str): row id of the ligand and protein files.

    """

    tool_key = "deeporigin.abfe-binding-fep"

    if database is None:
        database = _ensure_database(ABFE_DB)

    if params is None:
        params = _load_params("binding_fep")

    params["input"] = {
        "columnId": COL_MD_OUTPUT,
        "rowId": row_id,
        "databaseId": ABFE_DB,
    }

    outputs = {
        "output_file": {
            "columnId": COL_BINDING_FEP_OUTPUT,
            "rowId": row_id,
            "databaseId": database.hid,
        }
    }

    return run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )


@beartype
def render_mermaid_with_statuses(
    node_status: dict[str, str] = None,
) -> None:
    """
    Render a Mermaid diagram where each node is drawn as a rounded rectangle
    with a color indicating its status.

    Any node not specified in the node_status dict will default to "notStarted".

    """
    # Define the fixed nodes in the diagram.
    nodes = [
        "init",
        "emeq",
        "complex_prep",
        "ligand_prep",
        "solvation_fep",
        "simple_md",
        "binding_fep",
        "DeltaG",
    ]

    # Use an empty dictionary if none is provided.
    if node_status is None:
        node_status = {}

    # Build node definitions. For each node, use the provided status or default to "notStarted".
    node_defs = ""
    for node in nodes:
        status = node_status.get(node, "NotStarted")
        node_defs += f"    {node}({node}):::{status};\n"

    # Define the fixed edges of the diagram.
    edges = """
    init --> complex_prep;
    init --> ligand_prep;
    ligand_prep ----> solvation_fep;
    solvation_fep --> DeltaG;
    complex_prep --> emeq --> simple_md --> binding_fep --> DeltaG;
    """

    # Build the complete Mermaid diagram definition.
    mermaid_code = f"""
graph LR;
    %% Define styles for statuses:
    classDef NotStarted fill:#cccccc,stroke:#333,stroke-width:2px;
    classDef Queued fill:#cccccc,stroke:#222,stroke-width:2px;
    classDef Succeeded    fill:#90ee90,stroke:#333,stroke-width:2px;
    classDef Running       fill:#87CEFA,stroke:#333,stroke-width:2px;
    classDef Failed     fill:#ff7f7f,stroke:#333,stroke-width:2px;

{node_defs}
{edges}
    """

    # Render the diagram using your helper function.
    render_mermaid(mermaid_code)

    # Define HTML for the legend. Each status is displayed as a colored span.
    legend_html = """
    <div style="margin-top: 20px; font-family: sans-serif;">
      <span style="background-color:#cccccc; color: black; padding:2px 4px; margin: 0 8px;">NotStarted</span>
      <span style="background-color:#90ee90; color: black; padding:2px 4px; margin: 0 8px;">Suceedeed</span>
      <span style="background-color:#87CEFA; color: black; padding:2px 4px; margin: 0 8px;">Running</span>
      <span style="background-color:#ff7f7f; color: black; padding:2px 4px; margin: 0 8px;">Failed</span>
    </div>
    """
    # Display the legend below the Mermaid diagram.
    display(HTML(legend_html))


@beartype
def _run_job(
    *,
    ligand_file: str,
    row_id: str,
    df: pd.DataFrame,
    params: dict,
    func: Callable[..., Any],
    database: dict,
) -> pd.DataFrame:
    """utility function that runs a job if needed"""

    step = func.__name__

    existing = df[
        (df["ligand_file"] == ligand_file)
        & (df["step"] == step)
        & (df["Status"] != "Failed")
    ]
    if not existing.empty:
        print(f"Aborting for {ligand_file}: '{step}' step already exists.")
        return df

    data = api.get_cell_data(
        row_id=row_id,
        column_name=COL_JOB_IDS,
    )

    if data is None:
        data = {}

    else:
        data = ast.literal_eval(data)

    job_id = func(
        row_id=row_id,
        params=params,
        database=database,
    )

    data[step] = job_id

    api.set_cell_data(
        data,
        row_id=row_id,
        column_id=COL_JOB_IDS,
        database_id=ABFE_DB,
    )

    df.loc[len(df)] = [ligand_file, job_id, step, "Queued"]

    return df
