"""This module contains functions to start runs of various first-party tools on the Deep Origin platform"""

from typing import Optional

from beartype import beartype
from deeporigin.data_hub import api
from deeporigin.tools.utils import make_payload, run_tool


@beartype
def ligand_prep(
    *,
    database_id: str,
    row_id: str,
    output_column_name: str,
    ligand_column_name: str,
) -> str:
    """starts a run of Ligand Prep (Meeko) on the Deep Origin platform.

    Args:
        database_id (str): ID or HID of the database
        row_id (str): ID or HID of the row
        output_column_name (str): name of the column to write results to (the prepared ligand)
        ligand_column_name (str): name of the column in the database that contains the input ligand

    Returns:
        str: ID of the job. This ID can be used to query status.

    """

    db = api.describe_database(database_id=database_id)
    cols = db.cols

    inputs = dict(
        ligand=dict(
            columnId=ligand_column_name,
            databaseId=database_id,
            rowId=row_id,
        )
    )

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
        tool_id="deeporigin/ligand-prep",
        cols=cols,
    )

    response = run_tool(payload)

    execution_id = response.attributes.executionId
    job_id = response.id

    print(f"🧬 Job started with ID: {job_id}, execution ID: {execution_id}")
    return job_id


@beartype
def receptor_prep(
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
    """starts a run of Receptor Prep (PDBFixer) on the Deep Origin platform.

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
        str: ID of the job. This ID can be used to query status.

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


@beartype
def autodock_vina(
    *,
    database_id: str,
    row_id: str,
    search_space: dict,
    docking: dict,
    output_column_name: str,
    receptor_column_name: str,
    ligand_column_name: str,
) -> str:
    """starts an run of AutoDock Vina on the Deep Origin platform.

    Args:
        database_id (str): database ID or name of the database to source inputs from and write outputs to
        row_id (str, optional): row ID or name of the row to source inputs from and write outputs to.
        search_space (dict): search space parameters. Must include keys 'center_x', 'center_y', 'center_z', 'size_x', 'size_y', and 'size_z'
        docking (dict): docking parameters. Must include keys 'energy_range', 'exhaustiveness', and 'num_modes'
        output_column_name (str): name of the column to write output file to
        receptor_column_name (str): name of the column to source the receptor file from
        ligand_column_name (str): name of the column to source the ligand file from

    Returns:
        str: ID of the job. This ID can be used to query status.


    """

    if docking.keys() != {"energy_range", "exhaustiveness", "num_modes"}:
        raise ValueError(
            "docking must be a dictionary with keys 'energy_range', 'exhaustiveness', and 'num_modes'"
        )

    if search_space.keys() != {
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
    }:
        raise ValueError(
            "search_space must be a dictionary with keys 'center_x', 'center_y', 'center_z', 'size_x', 'size_y', and 'size_z'"
        )

    tool_id = "deeporigin/autodock-vina"
    inputs = dict(
        receptor={
            "rowId": row_id,
            "columnId": receptor_column_name,
            "databaseId": database_id,
        },
        ligand={
            "rowId": row_id,
            "columnId": ligand_column_name,
            "databaseId": database_id,
        },
        searchSpace=search_space,
        docking=docking,
    )
    outputs = dict(
        output_file={
            "rowId": row_id,
            "columnId": output_column_name,
            "databaseId": database_id,
        },
    )

    db = api.describe_database(database_id=database_id)
    cols = db.cols

    payload = make_payload(
        outputs=outputs,
        inputs=inputs,
        tool_id=tool_id,
        cols=cols,
    )

    response = run_tool(payload)

    execution_id = response.attributes.executionId
    job_id = response.id

    print(f"🧬 Job started with ID: {job_id}, execution ID: {execution_id}")
    return job_id
