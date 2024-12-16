"""module to run AutoDock Vina on Deep Origin"""

from beartype import beartype
from deeporigin.config import get_value
from deeporigin.data_hub import api
from deeporigin.platform import tools
from deeporigin.tools.utils import _get_cluster_id, _resolve_column_name


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

    if not database_id.startswith("_database"):
        data = api.convert_id_format(hids=[database_id])
        database_id = data[0].id

    db = api.describe_database(database_id=database_id)
    cols = db.cols

    # resolve columns
    output_column_id = _resolve_column_name(output_column_name, cols)
    receptor_column_id = _resolve_column_name(receptor_column_name, cols)
    ligand_column_id = _resolve_column_name(ligand_column_name, cols)

    if not row_id.startswith("_row"):
        data = data = api.convert_id_format(hids=[row_id])
        row_id = data[0].id

    json_data = {
        "toolId": "deeporigin/autodock-vina",
        "inputs": {
            "receptor": {
                "rowId": row_id,
                "columnId": receptor_column_id,
                "databaseId": database_id,
            },
            "ligand": {
                "rowId": row_id,
                "columnId": ligand_column_id,
                "databaseId": database_id,
            },
            "searchSpace": search_space,
            "docking": docking,
        },
        "outputs": {
            "output_file": {
                "rowId": row_id,
                "columnId": output_column_id,
                "databaseId": database_id,
            },
        },
        "clusterId": _get_cluster_id(),
    }

    response = tools.execute_tool(
        org_friendly_id=get_value()["organization_id"],
        execute_tool_dto=json_data,
    )

    job_id = response.id

    print(f"🧬 Job started with ID: {job_id}")
    return job_id
