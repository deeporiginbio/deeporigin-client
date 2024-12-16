"""this module contains utility functions used by tool execution"""

import functools
import time

from beartype import beartype
from deeporigin.config import get_value
from deeporigin.platform import clusters, tools


def query_run_status(job_id: str):
    """determine the status of a run, identified by job ID

    Args:
        job_id (str): job ID

    Returns:
        One of "Created", "Queued", "Running", "Succeeded", or "Failed"

    """

    data = tools.get_tool_execution(
        org_friendly_id=get_value()["organization_id"], resource_id=job_id
    )

    return data.attributes.status


def wait_for_job(job_id: str, *, poll_interval: int = 4):
    """Run while job is in non-terminal state, polling repeatedly"""

    status = "Running"
    txt_length = 0
    bs = "".join(["\b" for _ in range(txt_length)])
    while not (status == "Succeeded" or status == "Failed"):
        print(bs, end="", flush=True)
        status = query_run_status(job_id)
        txt_length = len(status)
        print(status, end="", flush=True)
        time.sleep(poll_interval)
        bs = "".join(["\b" for _ in range(txt_length)])


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
