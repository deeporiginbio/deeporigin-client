"""This module contains utility functions used by tool execution. In general, you will not need to use many of these functions directly."""

import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
import functools
import json
import os
import time
from typing import Any, Optional

from beartype import beartype
from box import Box

from deeporigin.config import get_value
from deeporigin.data_hub import api
from deeporigin.platform import clusters, tools
from deeporigin.platform.tools import execute_tool
from deeporigin.utils.core import _ensure_do_folder

JOBS_CACHE_DIR = _ensure_do_folder() / "jobs"

TERMINAL_STATES = {"Succeeded", "Failed"}
NON_TERMINAL_STATES = {"Created", "Queued", "Running"}


@beartype
def get_status_and_progress(execution_id: str) -> dict:
    """Determine the status of a run, identified by execution_id ID

    Args:
        execution_id (str): execution_id ID


    """

    data = tools.get_tool_execution(
        org_friendly_id=get_value()["organization_id"],
        execution_id=execution_id,
    )

    return dict(
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

    tools.action_tool_execution(
        org_friendly_id=get_value()["organization_id"],
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

    data = tools.get_tool_execution(
        org_friendly_id=get_value()["organization_id"],
        execution_id=execution_id,
    )

    return data.attributes.status


@beartype
def wait_for_job(
    execution_id: str,
    *,
    poll_interval: int = 4,
) -> None:
    """Repeatedly poll Deep Origin for the job status, till the status is "Succeeded" or "Failed (a terminal state)

    This function is useful for blocking execution of your code till a specific task is complete.

    Args:
        execution_id (str): execution_id ID. This is typically printed to screen and returned when a job is initialized.
        poll_interval (int, optional): number of seconds to wait between polling. Defaults to 4.

    """

    status = "Running"
    txt_length = 0
    bs = "".join(["\b" for _ in range(txt_length)])
    while not (status == "Succeeded" or status == "Failed"):
        print(bs, end="", flush=True)
        status = query_run_status(execution_id)
        txt_length = len(status)
        print(status, end="", flush=True)
        time.sleep(poll_interval)
        bs = "".join(["\b" for _ in range(txt_length)])


@beartype
def wait_for_jobs(
    refresh_time: int = 3,
    hide_succeeded: bool = True,
) -> Any:
    """Wait for all jobs started via this client to complete

    Args:
        refresh_time (int, optional): number of seconds to wait between polling. Defaults to 3.
        hide_succeeded (bool, optional): whether to hide jobs that have already completed. Defaults to True.

    Note that this function signature is explicitly not annotated with a return type to avoid importing pandas outside this function

    Returns:
        pd.DataFrame: dataframe of all jobs.

    """
    df = get_job_dataframe(update=True)

    if hide_succeeded:
        df = df[df["Status"] != "Succeeded"]

    from IPython.display import clear_output, display

    try:
        while len(set(df["Status"]).difference(TERMINAL_STATES)) != 0:
            df = get_job_dataframe(update=True)

            if hide_succeeded:
                df = df[df["Status"] != "Succeeded"]

            display(df)
            time.sleep(refresh_time)

            clear_output(wait=True)

    except KeyboardInterrupt:
        return

    print("✔️ All jobs completed")
    df = get_job_dataframe()
    return df


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

    return execute_tool(
        tool_key=tool_key,
        org_friendly_id=get_value()["organization_id"],
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
        dict: correctly formatted payload, ready to be passed to execute_tool
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
def get_job_dataframe(update: bool = False) -> Any:
    """returns a dataframe of all jobs and statuses, reading from local cache

    Args:
        update (bool, optional): Whether to check for updates on non-terminal jobs. Defaults to False.

    Note that this function is deliberately not annotated with an output type because pandas is imported internally to this funciton.

    Returns:
        pd.DataFrame: A dataframe containing job information"""
    jobs = tools.get_tool_executions(
        page=1,
        page_size=100,
        org_friendly_id="likely-aardvark-ewo",
        filter={},
        order="executionId",
    )

    if update:
        _update_all_jobs(jobs)

    import pandas as pd

    df = pd.DataFrame(
        {
            "Job ID": [job.id for job in jobs],
            "Execution ID": [job.attributes.executionId for job in jobs],
            "Status": [job.attributes.status for job in jobs],
            "Tool": [job.attributes.tool for job in jobs],
        }
    )
    return df


@beartype
def _update_all_jobs(jobs: list) -> None:
    """Update all job response files with the latest status, in parallel.

    This is an internal function. Do not use.
    """

    def should_update(job):
        return job.attributes.status in NON_TERMINAL_STATES

    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(query_run_status, job.id)
            for job in jobs
            if should_update(job)
        ]

        for future in futures:
            try:
                future.result()
            except Exception as e:
                print(f"An error occurred: {e}")


@beartype
def read_jobs() -> list:
    """read jobs from files in the jobs cache directory"""

    jobs = []

    if os.path.exists(JOBS_CACHE_DIR):
        for file_name in os.listdir(JOBS_CACHE_DIR):
            if file_name.endswith(".json"):
                file_path = os.path.join(JOBS_CACHE_DIR, file_name)
                try:
                    with open(file_path, "r") as f:
                        job_data = json.load(f)
                        jobs.append(Box(job_data))
                except Exception as e:
                    print(f"Failed to load {file_path}: {e}")
    else:
        print(f"Directory {JOBS_CACHE_DIR} does not exist.")

    return jobs


@beartype
def _ensure_database(name: str) -> dict:
    """ensure that a database exists with the given name. If it doesn't exist, create it"""

    databases = api.list_rows(row_type="database")

    database = [db for db in databases if db["hid"] == name]

    if len(database) == 0:
        # make a new DB
        print(f"🧬 Creating a database called {name}...")
        api.create_database(
            hid=name,
            hid_prefix=name,
            name=name,
        )

    database = api.describe_database(database_id=name)
    return database


@beartype
def _ensure_columns(
    *,
    database: dict,
    required_columns: list[dict],
):
    """ensure that columns exist with the given names (and types). If they don't exist, create them"""

    existing_column_names = []
    if "cols" in list(database.keys()):
        existing_column_names = [col["name"] for col in database.cols]

    # check if we need to make columns
    for item in required_columns:
        column_name = item["name"]
        column_type = item["type"]

        if column_name in existing_column_names:
            continue
        print(f"🧬 Making column named: {column_name} in {database.hid}")

        api.add_database_column(
            cardinality="one",
            database_id=database.hid,
            name=column_name,
            required=False,
            type=column_type,
        )

    database = api.describe_database(database_id=database.id)
    return database
