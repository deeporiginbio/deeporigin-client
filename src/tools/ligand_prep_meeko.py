"""module to run Ligand Prep (Meeko) on Deep Origin"""

from beartype import beartype
from deeporigin.data_hub import api
from deeporigin.tools.utils import make_payload, run_tool


@beartype
def start_run(
    *,
    database_id: str,
    row_id: str,
    output_column_name: str,
    ligand_column_name: str,
) -> str:
    """starts a run of Ligand Prep (Meeko) on Deep Origin

    Args:
        database_id (str): ID or HID of the database
        row_id (str): ID or HID of the row
        output_column_name (str): name of the column to write results to (the prepared ligand)
        ligand_column_name (str): name of the column in the database that contains the input ligand

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
