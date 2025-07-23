"""bridge module to interact with the platform tools api"""

import concurrent.futures
from functools import wraps
import inspect
import sys
from typing import Any, Callable, Literal, Optional

from beartype import beartype
import pandas as pd

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
    org_key: Optional[str] = None,
) -> dict:
    """Determine the status of a run, identified by job ID

    Args:
        execution_id (str): execution_id


    """

    data = get_tool_execution(  # noqa: F821
        execution_id=execution_id,
        org_key=org_key,
        client=client,
    )

    return data


@beartype
def get_statuses_and_progress(
    job_ids: list[str],
    *,
    client=None,
    org_key: Optional[str] = None,
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
                org_key=org_key,
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
    org_key: Optional[str] = None,
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
                org_key=org_key,
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
    org_key: Optional[str] = None,
) -> None:
    """cancel a run

    Args:
        execution_id (str): execution ID
    """

    data = get_status_and_progress(
        execution_id,
        client=client,
        org_key=org_key,
    )
    if data["status"] in ["Cancelled", "Failed", "Succeeded"]:
        return

    action_tool_execution(  # noqa: F821
        execution_id=execution_id,
        action="cancel",
        client=client,
        org_key=org_key,
    )


@beartype
def query_run_statuses(
    job_ids: list[str],
    *,
    client=None,
    org_key: Optional[str] = None,
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
                org_key=org_key,
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
    org_key: Optional[str] = None,
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
        org_key=org_key,
    )

    return data.attributes.status


@beartype
def run_tool(
    *,
    data: dict,
    tool_key: str,
    client=None,
    org_key: Optional[str] = None,
):
    """run any tool using provided data transfer object (DTO)

    Args:
        data (dict): data transfer object. This is typically generated by the `make_payload` function."""

    if "clusterId" not in data.keys():
        clusters = list_clusters(client=client, org_key=org_key)  # noqa: F821
        if len(clusters) == 0:
            raise RuntimeError("No clusters found.")
        cluster_id = clusters.data[0].id
        data["clusterId"] = cluster_id

    return execute_tool(  # noqa: F821
        tool_key=tool_key,
        execute_tool_schema_dto=data,
        client=client,
        org_key=org_key,
    )


@beartype
def _process_job(
    *,
    inputs: dict,
    outputs: dict,
    tool_key: str,
    metadata: Optional[dict] = None,
    client=None,
    org_key: Optional[str] = None,
) -> str:
    """helper function that uses inputs and outputs to construct a payload and run a tool"""

    payload = dict(
        outputs=outputs,
        inputs=inputs,
        metadata=metadata,
    )

    response = run_tool(
        data=payload,
        tool_key=tool_key,
        client=client,
        org_key=org_key,
    )
    job_id = response.id

    return job_id


def generate_tool_function(
    *,
    tool_key: str,
    tool_version: str,
    org_key: Optional[str] = None,
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
            org_key=org_key,
        )

    # Sanitize function name to be a valid Python identifier
    func_name = tool_key.replace(".", "_").replace("-", "_")
    generated_function.__name__ = func_name
    generated_function.__signature__ = sig

    return generated_function
