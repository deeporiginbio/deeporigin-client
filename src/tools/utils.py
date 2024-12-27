"""this module contains utility functions used by tool execution"""

import functools
import time
from typing import Optional

from beartype import beartype
from deeporigin.config import get_value
from deeporigin.platform import clusters, tools
from deeporigin.platform.tools import execute_tool


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


@beartype
def run_tool(data: dict):
    """run any tool using provided DTO"""

    if "clusterId" not in data.keys():
        data["clusterId"] = _get_cluster_id()

    return execute_tool(
        org_friendly_id=get_value()["organization_id"],
        execute_tool_dto=data,
    )


@beartype
def make_payload(
    *,
    inputs: dict,
    outputs: dict,
    tool_id: str,
    cluster_id: Optional[str] = None,
    cols: Optional[list] = None,
) -> dict:
    """helper function to create payload for tool execution"""

    if cluster_id is None:
        cluster_id = _get_cluster_id()

    payload = dict(
        inputs=inputs,
        outputs=outputs,
        clusterId=cluster_id,
        toolId=tool_id,
    )

    if cols:
        payload = _column_name_to_column_id(payload, cols)

    return payload


@beartype
def _column_name_to_column_id(data: dict, cols: list) -> dict:
    """
    Recursively update all values for the key 'columnId' in a nested dictionary.

    Args:
        d (dict): The dictionary to process.
        func (callable): A function that takes the original value and returns the new value.

    Returns:
        dict: The modified dictionary.
    """
    for key, value in data.items():
        if key == "columnId":
            data[key] = _resolve_column_name(value, cols)
        elif isinstance(value, dict):
            data[key] = _column_name_to_column_id(
                value, cols
            )  # Recurse into nested dictionaries
        elif isinstance(value, list):
            data[key] = [
                _column_name_to_column_id(item, cols)
                if isinstance(item, dict)
                else item
                for item in value
            ]
    return data
