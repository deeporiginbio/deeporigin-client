"""this module contains various functions to run steps of an ABFE workflow"""

import importlib.resources
import json
import os
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import List, Literal, Optional, Union

import numpy as np
import pandas as pd
from beartype import beartype
from box import Box
from deeporigin import chemistry
from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools import run
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from deeporigin.tools.utils import query_run_status
from deeporigin.utils.config import construct_resource_url
from deeporigin.utils.notebook import render_mermaid
from IPython.display import HTML, display

# constants and types
ABFE_DB = "ABFE"

charge_methods = Literal["gas", "bcc"]
"""Available charge methods"""

ligand_force_fields = Literal["gaff", "gaff2", "openff"]
"""Available ligand force fields. `gaff` is General Amber Force Field, gaff2 is an updated version of gaff. """

force_fields = Literal["ff14SB", "ff99SB-ildn"]
"""Available force fields"""

integrators = Literal[
    "BAOABIntegrator",
    "LangevinIntegrator",
    "SteepestDescentMinimizer",
    "LeapFrogIntegrator",
    "SoluteSolventSplittingIntegrator",
]
"""Integrator available for simulation"""


class PrettyDict(Box):
    """A dict subclass with a custom pretty-print representation."""

    def __repr__(self):
        return json.dumps(
            dict(self),
            indent=2,
            ensure_ascii=False,
        )

    def _repr_html_(self):
        self.__repr__()


@beartype
def _load_abfe_params() -> Box:
    """load default values for abfe end to end run"""
    with importlib.resources.open_text("deeporigin.json", "abfe.json") as f:
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
        dict(name="FEP ΔG (kcal/mol)", type="float"),
        dict(name="ligand_file", type="file"),
        dict(name="protein_file", type="file"),
        dict(name="complex_prep_output", type="file"),
        dict(name="ligand_prep_output", type="file"),
        dict(name="emeq_output", type="file"),
        dict(name="solvation_output", type="file"),
        dict(name="md_output", type="file"),
        dict(name="abfe_output", type="file"),
        dict(name="end_to_end_output", type="file"),
    ]

    database = _ensure_columns(
        database=database,
        required_columns=required_columns,
    )

    return database


@dataclass
class MDParams:
    integrator: integrators = "BAOABIntegrator"
    Δt: float = 0.004
    T: float = 298.15
    cutoff: float = 0.9
    fourier_spacing: float = 0.12
    hydrogen_mass: int = 2
    barostat: str = "MonteCarloBarostat"
    barostat_exchange_interval: int = 500


@dataclass
class SystemPrepParams:
    charge_method: charge_methods = "bcc"
    do_loop_modelling: bool = False
    force_field: force_fields = "ff14SB"
    is_lig_protonated: bool = True
    is_protein_protonated: bool = True
    keep_waters: bool = False
    lig_force_field: ligand_force_fields = "gaff2"
    padding: float = 1.0  # nm
    save_gmx_files: bool = False


@dataclass
class Ligand:
    file: Union[str, Path]
    smiles_string: Optional[str] = None
    is_protonated: Optional[bool] = False

    def __post_init__(self):
        if self.smiles_string is None:
            self.smiles_string = chemistry.sdf_to_smiles(self.file)

        if self.is_protonated is None:
            self.is_protonated = chemistry.is_ligand_protonated(self.file)

    def _repr_pretty_(self, p, cycle):
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
    file: str | Path


@dataclass
class ABFE:
    """ABFE class that can work with one protein and many ligands"""

    protein: Protein
    ligands: List[Ligand]

    delta_gs: List[float] = field(default_factory=list)
    row_ids: List[str] = field(default_factory=list)

    params: PrettyDict = field(default_factory=_load_abfe_params)

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
                column_id="ligand_file",
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
                column_id="protein_file",
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
                "FEP ΔG (kcal/mol)": self.delta_gs[idx],
                "Status": status,
            }
            rows.append(row)

        df = pd.DataFrame(
            rows,
            columns=[
                "Molecule",
                "FEP ΔG (kcal/mol)",
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
                output_column_name="end_to_end_output",
            )
        )

    def status(self):
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
