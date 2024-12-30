"""module to run Receptor Prep (pdbfixer) on Deep Origin"""

from typing import Optional

from beartype import beartype
from deeporigin.data_hub import api
from deeporigin.tools.utils import make_payload, run_tool


@beartype
def start_run(
    *,
    database_id: str,
    row_id: str,
    output_column_name: str,
    receptor_column_name: str,
    add_missing_residues: Optional[bool] = True,
    add_missing_atoms: Optional[bool] = True,
    add_missing_hydrogens: Optional[bool] = True,
    ph: Optional[float] = 7.4,
    remove_heterogens: Optional[bool] = True,
    remove_water: Optional[bool] = True,
) -> str:
    """starts a run of Receptor Prep (PDBFixer) on Deep Origin

    Args:
        database_id (str): ID or HID of the database
        row_id (str): ID or HID of the row
        output_column_name (str): name of the column to write results to (the prepared ligand)
        receptor_column_name (str): name of the column in the database that contains the input ligand
        add_missing_residues (bool, optional): whether to add missing residues. Defaults to True.
        add_missing_atoms (bool, optional): whether to add missing atoms. Defaults to True.
        add_missing_hydrogens (bool, optional): whether to add missing hydrogens. Defaults to True.
        ph (float, optional): pH value. Defaults to 7.4.
        remove_heterogens (bool, optional): whether to remove heterogens. Defaults to True.
        remove_water (bool, optional): whether to remove water. Defaults to True.

    Returns:
        str: ID of the job

    """

    db = api.describe_database(database_id=database_id)
    cols = db.cols

    inputs = {
        "receptor_pdb": dict(
            columnId=receptor_column_name,
            databaseId=database_id,
            rowId=row_id,
        ),
        "addMissingResidues": add_missing_residues,
        "addMissingAtoms": add_missing_atoms,
        "addMissingHydrogens": add_missing_hydrogens,
        "pH": ph,
        "removeHeterogens": remove_heterogens,
        "removeWater": remove_water,
        "nonStandardResidues": [],
    }

    outputs = dict(
        output_file=dict(
            rowId=row_id,
            columnId=output_column_name,
            databaseId=database_id,
        )
    )

    payload = make_payload(
        inputs=inputs,
        outputs=outputs,
        tool_id="deeporigin/receptor-prep",
        cols=cols,
    )
    response = run_tool(payload)

    execution_id = response.attributes.executionId
    job_id = response.id

    print(f"🧬 Job started with ID: {job_id}, execution ID: {execution_id}")
    return job_id
