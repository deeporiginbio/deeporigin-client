"""module to run AutoDock Vina on Deep Origin"""

import functools
import time

from beartype import beartype
from deeporigin.config import get_value
from deeporigin.data_hub import api
from deeporigin.platform import clusters, tools


@functools.cache
@beartype
def _get_cluster_id() -> str:
    """gets a valid cluster ID to run tools on

    this defaults to pulling us-west-2"""

    available_clusters = clusters.list_clusters(
        org_friendly_id=get_value()["organization_id"]
    )

    cluster = [
        cluster
        for cluster in available_clusters
        if "us-west-2" in cluster.attributes.name
    ]
    cluster = cluster[0]
    cluster_id = cluster.id
    return cluster_id


@beartype
def _resolve_column_name(column_name: str, cols: list) -> str:
    """resolve column IDs from column name

    Args:
        column_name (str): column name
        cols (list): list of columns

    Returns:
        str: column ID
    """

    column_ids = [col.id for col in cols]
    column_names = [col.name for col in cols]

    if column_name not in column_names and column_name not in column_ids:
        raise ValueError(f"column_name must be one of {column_names} or {column_ids}")
    elif column_name in column_names:
        column_id = [col.id for col in cols if col.name == column_name][0]
    else:
        column_id = column_name

    return column_id


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
        org_friendly_id="likely-aardvark-ewo",
        execute_tool_dto=json_data,
    )

    job_id = response.id

    print(f"🧬 Job started with ID: {job_id}")
    return job_id


def query_run_status(job_id: str):
    """determine the status of a run, identified by job ID

    Args:
        job_id (str): job ID

    Returns:
        One of "Created", "Queued", "Running", "Succeeded", or "Failed"

    """

    data = tools.get_tool_execution(
        org_friendly_id="likely-aardvark-ewo", resource_id=job_id
    )

    return data.attributes.status


def wait_for_job(job_id: str, *, poll_interval: int = 4):
    """Run while job is in non-terminal state, polling repeatedly"""

    status = "Running"
    txt_length = 0
    while not (status == "Succeeded" or status == "Failed"):
        status = query_run_status(job_id)
        txt_length = len(status)
        print(status, end="", flush=True)
        time.sleep(poll_interval)
        bs = "".join(["\b" for _ in range(txt_length)])
        print(bs, end="", flush=True)
