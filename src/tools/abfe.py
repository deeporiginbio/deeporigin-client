"""this module contains various functions to run steps of an ABFE workflow"""

import importlib.resources
import json
import os
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import List, Optional, Union

import numpy as np
import pandas as pd
from beartype import beartype
from box import Box
from deeporigin import chemistry
from deeporigin.data_hub import api
from deeporigin.tools import run
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from deeporigin.tools.utils import query_run_status
from deeporigin.utils.config import construct_resource_url
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


class PrettyDict(Box):
    """A dict subclass with a custom pretty-print representation."""

    def __repr__(self):
        """pretty print a dict"""
        return json.dumps(
            dict(self),
            indent=2,
            ensure_ascii=False,
        )

    def _repr_html_(self):
        """pretty print a dict"""
        self.__repr__()


@beartype
def _load_params(step: str) -> Box:
    """load default values for abfe end to end run"""

    with importlib.resources.open_text("deeporigin.json", f"{step}.json") as f:
        return PrettyDict(json.load(f))


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
        dict(name=COL_PROTEIN_FILE, type="file"),
        dict(name=COL_COMPLEX_PREP_OUTPUT, type="file"),
        dict(name=COL_LIGAND_PREP_OUTPUT, type="file"),
        dict(name=COL_EMEQ_OUTPUT, type="file"),
        dict(name=COL_SOLVATION_FEP_OUTPUT, type="file"),
        dict(name=COL_MD_OUTPUT, type="file"),
        dict(name=COL_BINDING_FEP_OUTPUT, type="file"),
        dict(name=COL_END_TO_END_OUTPUT, type="file"),
    ]

    database = _ensure_columns(
        database=database,
        required_columns=required_columns,
    )

    return database


@dataclass
class Ligand:
    """class to represent a ligand on disk"""

    file: Union[str, Path]
    smiles_string: Optional[str] = None
    is_protonated: Optional[bool] = False

    def __post_init__(self):
        """generates a SMILES if it doesn't exist"""

        if self.smiles_string is None:
            self.smiles_string = chemistry.sdf_to_smiles(self.file)

        if self.is_protonated is None:
            self.is_protonated = chemistry.is_ligand_protonated(self.file)

    def _repr_pretty_(self, p, cycle):
        """pretty print a ligand"""

        if cycle:
            p.text("Ligand(...)")
        else:
            p.text("Ligand(")

            with p.group(2, "\n  ", "\n"):
                all_fields = fields(self)
                for idx, field in enumerate(all_fields):
                    value = getattr(self, field.name)
                    p.text(f"{field.name}: {value!r}")
                    # Only add a breakable if this isn't the last field.
                    if idx < len(all_fields) - 1:
                        p.breakable()
            p.text(")")


@dataclass
class Protein:
    """class to represent a protein on disk"""

    file: str | Path


@dataclass
class ABFE:
    """ABFE class that can work with one protein and many ligands"""

    protein: Protein
    ligands: List[Ligand]

    delta_gs: List[float] = field(default_factory=list)
    row_ids: List[str] = field(default_factory=list)

    params: Optional[dict] = None

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
    def from_dir(cls, directory: str) -> "ABFE":
        """initialize an ABFE class given some files in a directory"""

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

    def init(self):
        """Initialize an ABFE run. Upload ligand(s) and protein files to Data Hub."""

        database = _ensure_db_for_abfe()

        # get params in
        self.params = dict()
        self.params.binding_fep = _load_params("binding_fep")
        self.params.complex_prep = _load_params("complex_prep")
        self.params.emeq = _load_params("emeq")
        self.params.ligand_prep = _load_params("ligand_prep")
        self.params.simple_md = _load_params("simple_md")
        self.params.solvation_fep = _load_params("solvation_fep")

        smiles_column_id = [col.id for col in database.cols if col.name == "Ligand"][0]

        url = construct_resource_url(
            name=ABFE_DB,
            row_type="database",
        )

        print(f"Using database at: {url}")
        print("🧬 Uploading files to database...")

        self.row_ids = []

        # now upload the protein
        protein_file = api.upload_file(file_path=self.protein.file)

        for ligand in self.ligands:
            # for each ligand, upload the ligand to a new row
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

            api.assign_files_to_cell(
                file_ids=[protein_file.id],
                database_id=ABFE_DB,
                column_id=COL_PROTEIN_FILE,
                row_id=row_id,
            )

            # update the dataframe
            self.df.loc[self.df["ligand_file"] == ligand.file, "Status"] = "Uploaded"

        print(f"🧬 Files uploaded to row {row_id}.")

    def results(self, image_size=(200, 200)):
        """
        Create and display a DataFrame in a Jupyter notebook with columns:
          - Molecule (an inline RDKit image from the SMILES)
          - SMILES (the raw SMILES string)
          - FEP ΔG (kcal/mol)
          - status
        """

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

        job_ids = []
        row_id = self.row_ids[0]
        job_ids.append(
            _run_e2e(
                row_id=row_id,
                steps=self.params.steps,
                output_column_name=COL_END_TO_END_OUTPUT,
            )
        )

    def status(self):
        """print the status of an ABFE run"""

        node_statuses = {
            "init": "NotStarted",
            "complex_prep": "NotStarted",
            "ligand_prep": "NotStarted",
            "solvation_FEP": "NotStarted",
            "simple_MD": "NotStarted",
            "binding_FEP": "NotStarted",
            "DeltaG": "NotStarted",
        }

        keys = self._job_ids.keys()

        for key in keys:
            if self._status[key] == "Succeeded" or self._status[key] == "Failed":
                continue

            # IN NON-terminal state. update
            self._status[key] = query_run_status(self._job_ids[key])

        node_statuses.update(self._status)

    def emeq(self):
        """run EMEQ step on each ligand"""

        for row_id in self.row_ids:
            emeq(
                row_id=row_id,
                params=self.params.emeq,
            )

    def complex_prep(self):
        """run complex prep step on each ligand"""

        for row_id in self.row_ids:
            complex_prep(
                row_id=row_id,
                params=self.params.complex_prep,
            )

    def ligand_prep(self):
        """run ligand prep step on each ligand"""

        for row_id in self.row_ids:
            ligand_prep(
                row_id=row_id,
                params=self.params.ligand_prep,
            )

    def simple_md(self):
        """run simple MD step on each ligand"""

        for row_id in self.row_ids:
            simple_md(
                row_id=row_id,
                params=self.params.simple_md,
            )

    def solvation_fep(self):
        """run solvation FEP step on each ligand"""

        for row_id in self.row_ids:
            solvation_fep(
                row_id=row_id,
                params=self.params.solvation_fep,
            )

    def binding_fep(self):
        """run binding FEP step on each ligand"""

        for row_id in self.row_ids:
            binding_fep(
                row_id=row_id,
                params=self.params.binding_fep,
            )


def _run_e2e(
    *,
    row_id: str,
    steps: dict,
    output_column_name: str,
) -> str:
    """Function to run an end-to-end job"""

    database = _ensure_db_for_abfe()

    url = construct_resource_url(
        name=ABFE_DB,
        row_type="database",
    )

    print(f"Using row {row_id} in database at: {url}")

    tool_key = "deeporigin.md-suite-abfe-e2e"

    inputs = {
        "ligand": {
            "columnId": "ligand_file",
            "rowId": row_id,
            "databaseId": database.hid,
        },
        "protein": {
            "columnId": "protein_file",
            "rowId": row_id,
            "databaseId": database.hid,
        },
        "steps": steps,
    }

    outputs = {
        "output_file": {
            "columnId": output_column_name,
            "rowId": row_id,
            "databaseId": database.hid,
        }
    }

    job_id = run._process_job(
        inputs=inputs,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )

    return job_id


@beartype
def emeq(
    *,
    row_id: str,
    params: Optional[dict] = None,
) -> None:
    """Run emeq on a ligand and protein pair, that exist as files on a row in the ABFE database. For this to work, the complex prep step must have been run first.

    Args:
        row_id (str): row id that contains the ligand and protein files.

    """

    tool_key = "deeporigin.md-suite-emeq"

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

    run._process_job(
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
) -> None:
    """Function to prepare uploaded Ligand and protein files using Deep Origin MDSuite. Use this function to run system prep on a ligand and protein pair, that exist as files on a row in the ABFE database.

    Args:
        row_id (str): row id that contains the ligand and protein files.


    """
    database = _ensure_db_for_abfe()

    tool_key = "deeporigin.md-suite-prep"

    if params is None:
        params = _load_params("ligand_prep")

    params["ligand"] = {
        "columnId": COL_LIGAND_FILE,
        "rowId": row_id,
        "databaseId": database.hid,
    }
    params["protein"] = {
        "columnId": "protein_file",
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

    run._process_job(
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
) -> None:
    """Function to prepare uploaded Ligand and protein files using Deep Origin MDSuite. Use this function to run system prep on a ligand and protein pair, that exist as files on a row in the ABFE database.

    Args:
        row_id (str): row id that contains the ligand and protein files.


    """
    database = _ensure_db_for_abfe()

    tool_key = "deeporigin.md-suite-prep"

    if params is None:
        params = _load_params("ligand_prep")

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

    run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )


def solvation_fep(
    *,
    row_id: str,
    params: Optional[dict] = None,
) -> None:
    """Run a solvation simulation


    Args:
        row_id (str): row id of the ligand and protein files.
    """

    database = _ensure_db_for_abfe()

    tool_key = "deeporigin.md-suite-solvation"

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

    run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )


def simple_md(
    *,
    row_id: str,
    params: Optional[dict] = None,
):
    """Run a simple MD simulation

    Args:
        row_id (str): row id of the ligand and protein files.
    """
    kwargs = locals()
    database = _ensure_db_for_abfe()

    tool_key = "deeporigin.md-suite-md"

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

    run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )


def binding_fep(
    *,
    row_id: str,
    params: Optional[dict] = None,
):
    """Run an ABFE simulation

    Args:
        row_id (str): row id of the ligand and protein files.

    """

    database = _ensure_db_for_abfe()

    tool_key = "deeporigin.md-suite-abfe"

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

    run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )
