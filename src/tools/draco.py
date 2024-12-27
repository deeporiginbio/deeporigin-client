"""module to run Draco on Deep Origin"""

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
    input_file_column_name: str,
    output_column_name: str,
) -> str:
    """starts an Draco run

    Args:
        database_id (str): database ID or name of the database to source inputs from and write outputs to
        row_id (str, optional): row ID or name of the row to source inputs from and write outputs to.
        input_file_column_name (str): name of the column to source the input file from
        output_column_name (str): name of the column to write output file to

    Returns:
        str: job ID of the run


    """

    if not database_id.startswith("_database"):
        data = api.convert_id_format(hids=[database_id])
        database_id = data[0].id

    db = api.describe_database(database_id=database_id)
    cols = db.cols

    # resolve columns
    output_column_id = _resolve_column_name(output_column_name, cols)
    input_file_column_id = _resolve_column_name(input_file_column_name, cols)

    json_data = {
        "toolId": "deeporigin/draco",
        "inputs": {
            "input_file": {
                "rowId": row_id,
                "columnId": input_file_column_id,
                "databaseId": database_id,
            },
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
