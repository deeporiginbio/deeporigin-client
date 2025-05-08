"""This module contains functions to start runs of various first-party tools on the Deep Origin platform"""

from beartype import beartype

from deeporigin.data_hub import api
from deeporigin.tools import utils


@beartype
def draco(
    *,
    database_id: str,
    row_id: str,
    input_file_column_name: str,
    output_column_name: str,
) -> str:
    """starts an Draco run

    Args:
        database_id (str): database ID or name of the database to source inputs from and write outputs to
        row_id (str): row ID or name of the row to source inputs from and write outputs to.
        input_file_column_name (str): name of the column to source the input file from
        output_column_name (str): name of the column to write output file to

    Returns:
        str: job ID of the run


    """

    inputs = dict(
        input_file={
            "rowId": row_id,
            "columnId": input_file_column_name,
            "databaseId": database_id,
        },
    )
    outputs = dict(
        output_file={
            "rowId": row_id,
            "columnId": output_column_name,
            "databaseId": database_id,
        },
    )

    db = api.describe_database(database_id=database_id)

    return utils._process_job(
        cols=db.cols,
        inputs=inputs,
        outputs=outputs,
        tool_key="deeporigin.draco",
    )
