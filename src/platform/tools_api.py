"""bridge module to interact with the platform tools api"""

import concurrent.futures
import functools
import sys
from typing import Optional

from beartype import beartype

from deeporigin.config import get_value
from deeporigin.platform import clusters_api, tools_api
from deeporigin.platform.utils import _add_functions_to_module

__all__ = _add_functions_to_module(
    module=sys.modules[__name__],
    api_name="ToolsApi",
)


TERMINAL_STATES = {"Succeeded", "Failed"}
NON_TERMINAL_STATES = {"Created", "Queued", "Running"}


@beartype
def get_status_and_progress(execution_id: str) -> dict:
    """Determine the status of a run, identified by job ID

    Args:
        execution_id (str): execution_id


    """

    data = tools_api.get_tool_execution(
        execution_id=execution_id,
    )

    return dict(
        job_id=execution_id,
        status=data.attributes.status,
        progress=data.attributes.progressReport,
        execution_id=data.attributes.executionId,  # confusingly, this is not the execution_id coming in
        inputs=data.attributes.userInputs,
        attributes=data.attributes,
    )


@beartype
def get_statuses_and_progress(job_ids: list[str]) -> list:
    """get statuses and progress reports for multiple jobs in parallel

    Args:
        job_ids (list[str]): list of job IDs

    """
    results = []

    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Submit all jobs and create a mapping from future to job_id
        future_to_job_id = {
            executor.submit(get_status_and_progress, job_id): job_id
            for job_id in job_ids
        }

        # As each future completes, store the result in the status dictionary
        for future in concurrent.futures.as_completed(future_to_job_id):
            try:
                results.append(future.result())
            except Exception:
                pass

    return results


@beartype
def cancel_runs(job_ids: list[str]) -> None:
    """Cancel multiple jobs in parallel.

    Args:
        job_ids: List of job IDs to cancel.
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(cancel_run, job_id) for job_id in job_ids]
        concurrent.futures.wait(futures)


@beartype
def cancel_run(execution_id: str) -> None:
    """cancel a run

    Args:
        execution_id (str): execution ID
    """

    data = get_status_and_progress(execution_id)
    if data["status"] in ["Cancelled", "Failed", "Succeeded"]:
        return

    tools_api.action_tool_execution(
        execution_id=execution_id,
        action="cancel",
    )


@beartype
def query_run_statuses(job_ids: list[str]) -> dict:
    """get statuses for multiple jobs in parallel

    Args:
        job_ids (list[str]): list of job IDs

    """
    status = {}

    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Submit all jobs and create a mapping from future to job_id
        future_to_job_id = {
            executor.submit(query_run_status, job_id): job_id for job_id in job_ids
        }

        # As each future completes, store the result in the status dictionary
        for future in concurrent.futures.as_completed(future_to_job_id):
            job_id = future_to_job_id[future]
            try:
                status[job_id] = future.result()
            except Exception:
                status[job_id] = None

    return status


@beartype
def query_run_status(execution_id: str) -> str:
    """Determine the status of a run, identified by execution_id ID

    Args:
        execution_id (str): execution_id ID

    Returns:
        One of "Created", "Queued", "Running", "Succeeded", or "Failed"

    """

    data = tools_api.get_tool_execution(execution_id=execution_id)

    return data.attributes.status


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
        raise ValueError(
            f"column_name must be one of {column_names} or {column_ids}. Instead it was: {column_name}"
        )
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

    available_clusters = clusters_api.list_clusters(
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
def run_tool(
    *,
    data: dict,
    tool_key: str,
):
    """run any tool using provided data transfer object (DTO)

    Args:
        data (dict): data transfer object. This is typically generated by the `make_payload` function."""

    if "clusterId" not in data.keys():
        data["clusterId"] = _get_cluster_id()

    return tools_api.execute_tool(
        tool_key=tool_key,
        execute_tool_dto=data,
    )


@beartype
def make_payload(
    *,
    inputs: dict,
    outputs: dict,
    cluster_id: Optional[str] = None,
    cols: Optional[list] = None,
    metadata: Optional[dict] = None,
) -> dict:
    """helper function to create payload for tool execution. This helper function is used by all wrapper functions in the run module to create the payload.

    Args:
        inputs (dict): inputs
        outputs (dict): outputs
        cluster_id (Optional[str], optional): cluster ID. Defaults to None. If not provided, the default cluster (us-west-2) is used.
        cols: (Optional[list], optional): list of columns. Defaults to None. If provided, column names (in inputs or outputs) are converted to column IDs.
        metadata: (Optional[dict], optional): metadata to be added to the payload. Defaults to None.

    Returns:
        dict: correctly formatted payload, ready to be passed to run_tool
    """

    if cluster_id is None:
        cluster_id = _get_cluster_id()

    payload = dict(
        inputs=inputs,
        outputs=outputs,
        clusterId=cluster_id,
        metadata=metadata,
    )

    if cols:
        payload = _column_name_to_column_id(payload, cols)

    payload = add_provider_if_databaseid_found(payload)

    return payload


def add_provider_if_databaseid_found(data):
    """
    Recursively traverse a data structure of nested dictionaries/lists.
    If a dict contains the key 'databaseId', add a peer key '$provider' = 'datahub-cell'.
    Return the modified data structure.
    """
    if isinstance(data, dict):
        # If this dict has 'databaseId', add '$provider'
        if "databaseId" in data:
            data["$provider"] = "datahub-cell"

        # Recursively check all values that are dicts or lists
        for _, value in data.items():
            if isinstance(value, dict):
                add_provider_if_databaseid_found(value)
            elif isinstance(value, list):
                for item in value:
                    add_provider_if_databaseid_found(item)

    elif isinstance(data, list):
        # If data is a list, recurse into each element if it is dict/list
        for item in data:
            add_provider_if_databaseid_found(item)

    return data


@beartype
def _column_name_to_column_id(
    data: dict,
    cols: list,
) -> dict:
    """
    Recursively update all values for the key 'columnId' in a nested dictionary.

    Args:
        data (dict): The dictionary to process.
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


@beartype
def _process_job(
    *,
    inputs: dict,
    outputs: dict,
    tool_key: str,
    cols: Optional[list] = None,
    metadata: Optional[dict] = None,
) -> str:
    """helper function that uses inputs and outputs to construct a payload and run a tool"""

    payload = make_payload(
        outputs=outputs,
        inputs=inputs,
        cols=cols,
        metadata=metadata,
    )

    response = run_tool(
        data=payload,
        tool_key=tool_key,
    )
    job_id = response.id

    return job_id
