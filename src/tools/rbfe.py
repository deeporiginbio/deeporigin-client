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
from IPython.display import HTML, display

# constants and types
RBFE_DB = "RBFE"
COL_DELTA_G = "FEP ΔG (kcal/mol)"
COL_LIGAND1_FILE = "ligand1_file"
COL_LIGAND2_FILE = "ligand2_file"
COL_PROTEIN_NAME = "protein_name"
COL_PROTEIN_FILE = "protein_file"
COL_END_TO_END_OUTPUT = "end_to_end_output"
COL_JOB_IDS = "job_ids"


RBFE_DIR = os.path.join(os.path.expanduser("~"), ".deeporigin", "abfe")
os.makedirs(RBFE_DIR, exist_ok=True)


@beartype
def _load_all_params() -> Box:
    """utility function to load all parameters"""

    params = Box()
    params.end_to_end = _load_params("rbfe_end_to_end")
    return params


@beartype
def _ensure_db_for_abfe() -> dict:
    """ensure that there is a database for FEP on Data Hub"""

    database = _ensure_database(RBFE_DB)

    try:
        api.add_smiles_column(
            name="Ligand",
            database_id=RBFE_DB,
        )
    except Exception:
        pass

    required_columns = [
        dict(name=COL_DELTA_G, type="float"),
        dict(name=COL_LIGAND1_FILE, type="file"),
        dict(name=COL_LIGAND2_FILE, type="file"),
        dict(name=COL_PROTEIN_NAME, type="text"),
        dict(name=COL_PROTEIN_FILE, type="file"),
        dict(name=COL_END_TO_END_OUTPUT, type="file"),
        dict(name=COL_JOB_IDS, type="text"),
    ]

    database = _ensure_columns(
        database=database,
        required_columns=required_columns,
    )

    return database


@dataclass
class RBFE:
    """RBFE class that can work with one protein and many ligands"""

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

    @classmethod
    def from_dir(cls, directory: str) -> "RBFE":
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

        # Create the RBFE instance
        rbfe = cls(
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
        rbfe.df = pd.DataFrame(rows, columns=rbfe.df.columns)

        return rbfe

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
            name=RBFE_DB,
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
                database_id=RBFE_DB,
                column_id=COL_LIGAND_FILE,
                file_path=ligand.file,
            )

            row_id = response.rows[0].hid

            self.row_ids.append(row_id)

            api.set_cell_data(
                ligand.smiles_string,
                database_id=RBFE_DB,
                row_id=row_id,
                column_id=smiles_column_id,
            )

            api.set_cell_data(
                self.protein.name,
                database_id=RBFE_DB,
                row_id=row_id,
                column_id=protein_name_column_id,
            )

            api.assign_files_to_cell(
                file_ids=[protein_file.id],
                database_id=RBFE_DB,
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

    def results(self, image_size=(200, 100)):
        """
        Create and display a DataFrame in a Jupyter notebook with columns:
          - Molecule (an inline RDKit image from the SMILES)
          - SMILES (the raw SMILES string)
          - FEP ΔG (kcal/mol)
          - status
        """

        self.update()

        df = self.df.copy()

        # 1. Mark rows where jobID is non-empty
        df["has_jobID"] = df["jobID"].notna() & (df["jobID"] != "")

        # 2. Sort by 'has_jobID' in descending order,
        #    ensuring rows with a valid jobID appear first
        df = df.sort_values("has_jobID", ascending=False)

        # 3. Drop duplicates on 'ligand_file', keeping the first row for each group
        df = df.drop_duplicates(subset=["ligand_file"], keep="first")

        # 4. Drop the helper column
        df = df.drop(columns=["has_jobID"])

        rows = []
        for _, row in df.iterrows():
            ligand_file = row["ligand_file"]

            idx = [lig.file for lig in self.ligands].index(ligand_file)

            ligand = self.ligands[idx]

            status = row["Status"]

            row = {
                "Molecule": chemistry.smiles_to_base64_png(
                    ligand.smiles_string,
                    size=image_size,
                ),
                COL_DELTA_G: self.delta_gs[idx],
                "Status": status,
            }
            rows.append(row)

        df = pd.DataFrame(
            rows,
            columns=[
                "Molecule",
                COL_DELTA_G,
                "Status",
            ],
        )

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


def end_to_end(
    *,
    row_id: str,
    params: Optional[dict] = None,
    database: Optional[dict] = None,
) -> str:
    """Function to run an end-to-end job"""

    tool_key = "deeporigin.rbfe-end-to-end"

    if database is None:
        database = _ensure_database(RBFE_DB)

    if params is None:
        params = _load_params("rbfe_end_to_end")

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
        }
    }

    return run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )


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
        database_id=RBFE_DB,
    )

    df.loc[len(df)] = [ligand_file, job_id, step, "Queued"]

    return df
