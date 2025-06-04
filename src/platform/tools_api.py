"""bridge module to interact with the platform tools api"""

import concurrent.futures
from functools import wraps
import inspect
import sys
from typing import Any, Callable, Literal, Optional

from beartype import beartype
import pandas as pd

from deeporigin.platform import clusters_api
from deeporigin.platform.utils import _add_functions_to_module

__all__ = _add_functions_to_module(
    module=sys.modules[__name__],
    api_name="ToolsApi",
)


TERMINAL_STATES = {"Succeeded", "Failed", "Cancelled"}
NON_TERMINAL_STATES = {"Created", "Queued", "Running"}

NON_FAILED_STATES = {"Succeeded", "Running", "Queued", "Created"}

# possible providers for files that work with the tools API
PROVIDER = Literal["ufa", "s3"]


@beartype
def get_status_and_progress(
    execution_id: str,
    *,
    client=None,
    org_friendly_id: Optional[str] = None,
) -> dict:
    """Determine the status of a run, identified by job ID

    Args:
        execution_id (str): execution_id


    """

    data = get_tool_execution(  # noqa: F821
        execution_id=execution_id,
        org_friendly_id=org_friendly_id,
        client=client,
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
def get_statuses_and_progress(
    job_ids: list[str],
    *,
    client=None,
    org_friendly_id: Optional[str] = None,
) -> list:
    """get statuses and progress reports for multiple jobs in parallel

    Args:
        job_ids (list[str]): list of job IDs

    """
    results = []
    errors = []

    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Submit all jobs and create a mapping from future to job_id
        future_to_job_id = {
            executor.submit(
                get_status_and_progress,
                job_id,
                client=client,
                org_friendly_id=org_friendly_id,
            ): job_id
            for job_id in job_ids
        }

        # As each future completes, store the result in the status dictionary
        for future in concurrent.futures.as_completed(future_to_job_id):
            job_id = future_to_job_id[future]
            try:
                results.append(future.result())
            except Exception as e:
                errors.append((job_id, e))

    if errors:
        error_msgs = "\n".join([f"Job {jid}: {str(err)}" for jid, err in errors])
        raise RuntimeError(
            f"Some jobs failed in get_status_and_progress:\n{error_msgs}"
        )

    return results


@beartype
def cancel_runs(
    job_ids: list[str] | pd.core.series.Series | pd.core.frame.DataFrame,
    *,
    client=None,
    org_friendly_id: Optional[str] = None,
) -> None:
    """Cancel multiple jobs in parallel.

    Args:
        job_ids: List of job IDs to cancel.
    """

    if isinstance(job_ids, pd.core.series.Series):
        job_ids = job_ids.tolist()
    elif isinstance(job_ids, pd.core.frame.DataFrame):
        job_ids = job_ids["id"].tolist()

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(
                cancel_run,
                job_id,
                client=client,
                org_friendly_id=org_friendly_id,
            )
            for job_id in job_ids
        ]
        concurrent.futures.wait(futures)
        # Raise the first exception encountered
        for future in futures:
            exc = future.exception()
            if exc is not None:
                raise exc


@beartype
def cancel_run(
    execution_id: str,
    *,
    client=None,
    org_friendly_id: Optional[str] = None,
) -> None:
    """cancel a run

    Args:
        execution_id (str): execution ID
    """

    data = get_status_and_progress(
        execution_id,
        client=client,
        org_friendly_id=org_friendly_id,
    )
    if data["status"] in ["Cancelled", "Failed", "Succeeded"]:
        return

    action_tool_execution(  # noqa: F821
        execution_id=execution_id,
        action="cancel",
        client=client,
        org_friendly_id=org_friendly_id,
    )


@beartype
def query_run_statuses(
    job_ids: list[str],
    *,
    client=None,
    org_friendly_id: Optional[str] = None,
) -> dict:
    """get statuses for multiple jobs in parallel

    Args:
        job_ids (list[str]): list of job IDs

    """
    status = {}

    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Submit all jobs and create a mapping from future to job_id
        future_to_job_id = {
            executor.submit(
                query_run_status,
                job_id,
                client=client,
                org_friendly_id=org_friendly_id,
            ): job_id
            for job_id in job_ids
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
def query_run_status(
    execution_id: str,
    *,
    client=None,
    org_friendly_id: Optional[str] = None,
) -> str:
    """Determine the status of a run, identified by execution_id ID

    Args:
        execution_id (str): execution_id ID

    Returns:
        One of "Created", "Queued", "Running", "Succeeded", or "Failed"

    """

    data = get_tool_execution(  # noqa: F821
        execution_id=execution_id,
        client=client,
        org_friendly_id=org_friendly_id,
    )

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


@beartype
def run_tool(
    *,
    data: dict,
    tool_key: str,
    client=None,
    org_friendly_id: Optional[str] = None,
    clusters_client=None,
):
    """run any tool using provided data transfer object (DTO)

    Args:
        data (dict): data transfer object. This is typically generated by the `make_payload` function."""

    if "clusterId" not in data.keys():
        data["clusterId"] = clusters_api._get_cluster_id(
            client=clusters_client,
            org_friendly_id=org_friendly_id,
        )

    return execute_tool(  # noqa: F821
        tool_key=tool_key,
        execute_tool_dto=data,
        client=client,
        org_friendly_id=org_friendly_id,
    )


@beartype
def _process_job(
    *,
    inputs: dict,
    outputs: dict,
    tool_key: str,
    cols: Optional[list] = None,
    metadata: Optional[dict] = None,
    client=None,  # this is the tools client
    org_friendly_id: Optional[str] = None,
    clusters_client=None,
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
        client=client,
        org_friendly_id=org_friendly_id,
        clusters_client=clusters_client,
    )
    job_id = response.id

    return job_id


@beartype
def make_payload(
    *,
    inputs: dict,
    outputs: dict,
    cols: Optional[list] = None,
    metadata: Optional[dict] = None,
) -> dict:
    """helper function to create payload for tool execution. This helper function is used by all wrapper functions in the run module to create the payload.

    Args:
        inputs (dict): inputs
        outputs (dict): outputs
        cols: (Optional[list], optional): list of columns. Defaults to None. If provided, column names (in inputs or outputs) are converted to column IDs.
        metadata: (Optional[dict], optional): metadata to be added to the payload. Defaults to None.

    Returns:
        dict: correctly formatted payload, ready to be passed to run_tool
    """

    # note that clusterId is not included here
    # it need to be included for this payload to be valid
    payload = dict(
        inputs=inputs,
        outputs=outputs,
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


def generate_tool_function(
    *,
    tool_key: str,
    tool_version: str,
    org_friendly_id: Optional[str] = None,
) -> Callable:
    """utility function that generates a function that can be used to run a too, from schema

    Example usage:

        docking_func = tools_api.generate_tool_function(tool_key="deeporigin.bulk-docking", tool_version="0.3.0")
        docking_func(
            protein_path="path/to/protein.pdb",
            ligand_path="path/to/ligand.sdf",
            ...
        )

    """

    data = get_tool(tool_key=tool_key)  # noqa: F821
    data = [item for item in data if item.version == tool_version]

    if len(data) == 0:
        raise ValueError(
            f"No tool found for tool_key={tool_key} and tool_version={tool_version}"
        )
    elif len(data) > 1:
        raise ValueError(
            f"Multiple tools found for tool_key={tool_key} and tool_version={tool_version}"
        )

    schema = data[0]

    inputs_schema = schema["inputs"]
    outputs_schema = schema["outputs"]

    required_inputs = inputs_schema.get("required", [])
    required_outputs = outputs_schema.get("required", [])

    input_props = inputs_schema["properties"]

    # --- Construct argument specification
    parameters = []

    for name in required_inputs:
        prop = input_props[name]
        annotation = (
            str
            if prop.get("format") == "file"
            else (
                list[float]
                if prop.get("type") == "array" and prop["items"]["type"] == "number"
                else (
                    list[str]
                    if prop.get("type") == "array" and prop["items"]["type"] == "string"
                    else Any
                )
            )
        )
        parameters.append(
            inspect.Parameter(
                name, inspect.Parameter.POSITIONAL_OR_KEYWORD, annotation=annotation
            )
        )

    for name in required_outputs:
        parameters.append(
            inspect.Parameter(
                name, inspect.Parameter.POSITIONAL_OR_KEYWORD, annotation=str
            )
        )

    parameters.append(
        inspect.Parameter(
            "metadata",
            inspect.Parameter.KEYWORD_ONLY,
            default=None,
            annotation=dict[str, Any] | None,
        )
    )
    parameters.append(
        inspect.Parameter(
            "client",
            inspect.Parameter.KEYWORD_ONLY,
            default=None,
            annotation=Any,
        )
    )

    sig = inspect.Signature(parameters)

    @wraps(lambda *args, **kwargs: None)
    def generated_function(*args, **kwargs):
        bound = sig.bind(*args, **kwargs)
        bound.apply_defaults()
        arguments = bound.arguments

        inputs = {}
        for name in required_inputs:
            val = arguments[name]
            if input_props[name].get("format") == "file":
                inputs[name] = {"$provider": "ufa", "key": val}
            else:
                inputs[name] = val

        outputs = {
            name: {"$provider": "ufa", "key": arguments[name]}
            for name in required_outputs
        }

        metadata = arguments.get("metadata")
        client = arguments.get("client")

        return _process_job(
            inputs=inputs,
            outputs=outputs,
            tool_key=tool_key,
            metadata=metadata,
            client=client,
            org_friendly_id=org_friendly_id,
        )

    # Sanitize function name to be a valid Python identifier
    func_name = tool_key.replace(".", "_").replace("-", "_")
    generated_function.__name__ = func_name
    generated_function.__signature__ = sig

    return generated_function
