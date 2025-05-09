"""This module contains utility functions for the Drug Discovery module"""

import importlib.resources
import json
import os
from typing import Any, Literal, Optional

from beartype import beartype

from deeporigin.data_hub import api
from deeporigin.utils.core import PrettyDict, hash_strings

# constants
DB_ABFE = "ABFE"
DB_RBFE = "RBFE"
DB_DOCKING = "Docking"
DB_PROTEINS = "Proteins"
DB_LIGANDS = "Ligands"

VALID_TOOLS = Literal[DB_ABFE, DB_RBFE, DB_DOCKING]

COL_DELTA_DELTA_G = "FEP ΔΔG (kcal/mol)"
COL_DELTA_G = "FEP ΔG (kcal/mol)"
COL_JOBID = "JobID"
COL_LIGAND1 = "Ligand1"  # for ABFE/RBFE
COL_LIGAND2 = "Ligand2"  # for RBFE
COL_CSV_OUTPUT = "OutputFile"  # this will be a CSV
COL_PROTEIN = "Protein"
COL_LIGAND = "Ligand"
COL_RESULT = "ResultFile"
COL_SMILES_HASH = "SMILESHash"
COL_COMPLEX_HASH = "ComplexHash"
COL_STEP = "Step"


DATA_DIRS = dict()

DATA_DIRS[DB_ABFE] = os.path.join(os.path.expanduser("~"), ".deeporigin", DB_ABFE)
DATA_DIRS[DB_RBFE] = os.path.join(os.path.expanduser("~"), ".deeporigin", DB_RBFE)
DATA_DIRS[DB_DOCKING] = os.path.join(os.path.expanduser("~"), ".deeporigin", DB_DOCKING)


os.makedirs(DATA_DIRS[DB_ABFE], exist_ok=True)
os.makedirs(DATA_DIRS[DB_RBFE], exist_ok=True)
os.makedirs(DATA_DIRS[DB_DOCKING], exist_ok=True)


@beartype
def _load_params(tool: str) -> PrettyDict:
    """load params for various tools, reading from JSON files"""

    with importlib.resources.open_text("deeporigin.json", f"{tool}.json") as f:
        return PrettyDict(json.load(f))


@beartype
def _start_tool_run(
    *,
    params: dict,
    database_columns: list,
    tool: VALID_TOOLS,
    protein_id: str,
    complex_hash: str,
    ligand1_id: Optional[str] = None,
    ligand2_id: Optional[str] = None,
) -> str:
    """starts a single run of ABFE end to end and logs it in the ABFE database. Internal function. Do not use.

    Args:
        protein_id (str): protein ID
        ligand_id (str): ligand ID
        params (dict): parameters for the ABFE end-to-end job
        database_columns (list): list of database columns dicts

    """

    # input validation
    if tool == "ABFE" and ligand1_id is None:
        raise ValueError("ligand1_id is required for ABFE")

    if tool == "RBFE" and (ligand1_id is None or ligand2_id is None):
        raise ValueError("ligand1_id and ligand2_id is required for RBFE")

    # tool key mapper
    tool_key_mapper = dict(
        ABFE="deeporigin.abfe-end-to-end/0.2.5",
        RBFE="deeporigin.rbfe-end-to-end",
        Docking="deeporigin.bulk-docking/0.3.0",
    )

    from deeporigin.tools import run

    # make a new row
    response = api.make_database_rows(tool, n_rows=1)
    row_id = response.rows[0].hid

    # a protein is needed for ABFE, RBFE, and docking
    params["protein"] = {
        "columnId": COL_PROTEIN,
        "rowId": protein_id,
        "databaseId": DB_PROTEINS,
    }

    # input ligand files
    if tool == "RBFE":
        params["ligand1"] = {
            "columnId": COL_LIGAND,
            "rowId": ligand1_id,
            "databaseId": DB_LIGANDS,
        }

        params["ligand2"] = {
            "columnId": COL_LIGAND,
            "rowId": ligand2_id,
            "databaseId": DB_LIGANDS,
        }
    elif tool == "ABFE":
        params["ligand"] = {
            "columnId": COL_LIGAND,
            "rowId": ligand1_id,
            "databaseId": DB_LIGANDS,
        }

    # output files
    if tool == "RBFE":
        outputs = {
            "output_file": {
                "columnId": COL_RESULT,
                "rowId": row_id,
                "databaseId": DB_RBFE,
            },
            "rbfe_results_summary": {
                "columnId": COL_CSV_OUTPUT,
                "rowId": row_id,
                "databaseId": DB_RBFE,
            },
        }
    elif tool == "ABFE":
        outputs = {
            "output_file": {
                "columnId": COL_RESULT,
                "rowId": row_id,
                "databaseId": DB_ABFE,
            },
            "abfe_results_summary": {
                "columnId": COL_CSV_OUTPUT,
                "rowId": row_id,
                "databaseId": DB_ABFE,
            },
        }
    elif tool == "Docking":
        outputs = {
            "data_file": {
                "columnId": COL_CSV_OUTPUT,
                "rowId": row_id,
                "databaseId": DB_DOCKING,
            },
            "results_sdf": {
                "columnId": COL_RESULT,
                "rowId": row_id,
                "databaseId": DB_DOCKING,
            },
        }

    if is_test_run(params):
        print("⚠️ Warning: test_run=1 in these parameters. Results may not be accurate.")

    metadata = dict(
        protein_id=protein_id,
        complex_hash=complex_hash,
        ligand1_id=ligand1_id,
        ligand2_id=ligand2_id,
    )

    job_id = run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key_mapper[tool],
        cols=database_columns,
        metadata=metadata,
    )

    # write job ID
    api.set_cell_data(
        job_id,
        column_id=COL_JOBID,
        row_id=row_id,
        database_id=tool,
    )

    # write ligand1 ID
    if ligand1_id is not None:
        api.set_cell_data(
            ligand1_id,
            column_id=COL_LIGAND1,
            row_id=row_id,
            database_id=tool,
        )

    # write ligand2 ID
    if ligand2_id is not None:
        api.set_cell_data(
            ligand2_id,
            column_id=COL_LIGAND2,
            row_id=row_id,
            database_id=tool,
        )

    # write protein ID
    api.set_cell_data(
        protein_id,
        column_id=COL_PROTEIN,
        row_id=row_id,
        database_id=tool,
    )

    # write complex hash
    api.set_cell_data(
        complex_hash,
        column_id=COL_COMPLEX_HASH,
        row_id=row_id,
        database_id=tool,
    )

    # write SMILEShash
    if tool == "Docking":
        smiles_hash = smiles_hash = hash_strings(params["smiles_list"])
        api.set_cell_data(
            smiles_hash,
            column_id=COL_SMILES_HASH,
            row_id=row_id,
            database_id=tool,
        )

    # write step
    if tool in ["ABFE", "RBFE"]:
        api.set_cell_data(
            "End-to-end",
            column_id=COL_STEP,
            row_id=row_id,
            database_id=tool,
        )

    return job_id


@beartype
def is_test_run(data: Any) -> bool:
    """check if test_run=1 in a dict"""

    if isinstance(data, dict):
        if data.get("test_run") == 1:
            return True
        for value in data.values():
            if is_test_run(value):
                return True
    elif isinstance(data, list):
        for item in data:
            if is_test_run(item):
                return True
    return False
