"""module to run AutoDock Vina on Deep Origin"""

from beartype import beartype
from deeporigin.data_hub import api
from deeporigin.tools.utils import make_payload, run_tool


@beartype
def start_run(
    *,
    database_id: str,
    row_id: str,
    search_space: dict,
    docking: dict,
    output_column_name: str,
    receptor_column_name: str,
    ligand_column_name: str,
) -> str:
    """starts an AutoDock Vina run

    Args:
        database_id (str): database ID or name of the database to source inputs from and write outputs to
        row_id (str, optional): row ID or name of the row to source inputs from and write outputs to.
        search_space (dict): search space parameters. Must include keys 'center_x', 'center_y', 'center_z', 'size_x', 'size_y', and 'size_z'
        docking (dict): docking parameters. Must include keys 'energy_range', 'exhaustiveness', and 'num_modes'
        output_column_name (str): name of the column to write output file to
        receptor_column_name (str): name of the column to source the receptor file from
        ligand_column_name (str): name of the column to source the ligand file from

    Returns:
        str: job ID of the run


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
