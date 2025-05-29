"""This module contains utility functions for the Drug Discovery module"""

import importlib.resources
import json
import os
from typing import Any, Optional

from beartype import beartype

from deeporigin.drug_discovery.constants import tool_mapper, valid_tools
from deeporigin.files import FilesClient
from deeporigin.platform import tools_api
from deeporigin.platform.utils import PlatformClients
from deeporigin.utils.core import PrettyDict

DATA_DIRS = dict()

for tool in tool_mapper.keys():
    DATA_DIRS[tool] = os.path.join(os.path.expanduser("~"), ".deeporigin", tool)
    os.makedirs(DATA_DIRS[tool], exist_ok=True)


@beartype
def _load_params(param_file: str) -> PrettyDict:
    """load params for various tools, reading from JSON files"""

    with importlib.resources.open_text("deeporigin.json", f"{param_file}.json") as f:
        return PrettyDict(json.load(f))


@beartype
def _start_tool_run(
    *,
    params: dict,
    metadata: dict,
    protein_path: str,
    tool: valid_tools,
    tool_version: str,
    ligand1_path: Optional[str] = None,
    ligand2_path: Optional[str] = None,
    provider: tools_api.PROVIDER = "ufa",
    _platform_clients: Optional[PlatformClients] = None,
    _output_dir_path: Optional[str] = None,
) -> str:
    """
    Starts a single run of an end-to-end tool (such as ABFE) and logs it in the ABFE database.

    This is an internal function that prepares input and output file parameters, sets up the job metadata,
    and submits the job to the platform's tools API. Only ABFE is currently supported.

    Args:
        params (dict): Parameters for the tool run, including input and configuration options.
        metadata (dict): Metadata to be logged with the job.
        protein_path (str): Remote path to the protein file to be used in the run.
        ligand1_path (str): Remnote path to the first ligand file.
        ligand2_path (Optional[str]): Remote path to the second ligand file (required for RBFE).
        tool (valid_tools): The tool to run (e.g., 'ABFE', 'RBFE').
        tool_version (str): Version of the tool to use.
        provider (tools_api.PROVIDER, optional): File provider for input/output files. Defaults to 'ufa'.
        _platform_clients (Optional[PlatformClients]): Platform client objects for API access.
        _output_dir_path (Optional[str]): Custom output directory path (on remote storage). If None, a default is constructed.

    Returns:
        str: The job ID of the started tool run.

    Raises:
        NotImplementedError: If a tool other than ABFE is specified.
    """

    if _output_dir_path is None:
        if tool == "ABFE":
            _output_dir_path = (
                "tool-runs/"
                + tool
                + "/"
                + os.path.basename(protein_path)
                + "/"
                + os.path.basename(ligand1_path)
                + "/"
            )
        else:
            raise NotImplementedError("Tools other than ABFE are not implemented yet")

    # a protein is needed for ABFE, RBFE, and docking
    params["protein"] = {
        "$provider": provider,
        "key": protein_path,
    }

    # input ligand files
    if tool == "RBFE":
        params["ligand1"] = {
            "$provider": provider,
            "key": ligand1_path,
        }

        params["ligand2"] = {
            "$provider": provider,
            "key": ligand2_path,
        }
    elif tool == "ABFE":
        params["ligand"] = {
            "$provider": provider,
            "key": ligand1_path,
        }

    # output files
    if tool == "RBFE":
        outputs = {
            "output_file": {
                "$provider": provider,
                "key": _output_dir_path + "output/",
            },
            "rbfe_results_summary": {
                "$provider": provider,
                "key": _output_dir_path + "results.csv",
            },
        }
    elif tool == "ABFE":
        outputs = {
            "output_file": {
                "$provider": provider,
                "key": _output_dir_path + "output/",
            },
            "abfe_results_summary": {
                "$provider": provider,
                "key": _output_dir_path + "results.csv",
            },
        }
    elif tool == "Docking":
        outputs = {
            "data_file": {
                "$provider": provider,
                "key": _output_dir_path + "results.csv",
            },
            "results_sdf": {
                "$provider": provider,
                "key": _output_dir_path + "results.sdf",
            },
        }

    if is_test_run(params):
        print(
            "⚠️ Warning: test_run=1 in these parameters. Results will not be accurate."
        )

    tools_client = getattr(_platform_clients, "ToolsApi", None)
    clusters_client = getattr(_platform_clients, "ClustersApi", None)

    if _platform_clients is None:
        from deeporigin.config import get_value

        org_friendly_id = get_value()["organization_id"]
    else:
        org_friendly_id = _platform_clients.org_friendly_id

    job_id = tools_api._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_mapper[tool],
        metadata=metadata,
        client=tools_client,
        org_friendly_id=org_friendly_id,
        clusters_client=clusters_client,
    )

    print(f"Job started with Job ID: {job_id}")

    return job_id


@beartype
def is_test_run(data: Any) -> bool:
    """check if test_run=1 in a dict"""

    if isinstance(data, dict):
        if data.get("test_run") == 1:
            return True
        for value in data.values():
            if is_test_run(value):
                return True
    elif isinstance(data, list):
        for item in data:
            if is_test_run(item):
                return True
    return False


@beartype
def _set_test_run(data, value: int = 1) -> None:
    """recursively iterate over a dict and set test_run=1 for all keys"""

    if isinstance(data, dict):
        for key, val in data.items():
            if key == "test_run":
                data[key] = value
            else:
                _set_test_run(val, value)
    elif isinstance(data, list):
        for item in data:
            _set_test_run(item, value)


@beartype
def find_files_on_ufa(
    *,
    tool: str,
    protein: str,
    ligand: Optional[str] = None,
    client: Optional[FilesClient] = None,
) -> list:
    """
    Find files on the UFA (Unified File API) storage for a given tool run.

    This function searches for files associated with a specific tool, protein, and optionally ligand
    in the UFA storage system. It constructs the appropriate search path based on the provided arguments
    and returns a list of file paths found under that directory.

    Args:
        tool (str): The name of the tool (e.g., 'ABFE', 'RBFE', etc.).
        protein (str): The protein identifier or filename used in the tool run.
        ligand (Optional[str]): The ligand identifier or filename used in the tool run. If not provided,
            the search will be performed at the protein level only.

    Returns:
        List[str]: A list of file paths found in the specified UFA directory.
    """

    if client is None:
        from deeporigin.files import FilesClient

        client = FilesClient()

    if ligand is not None:
        search_str = f"tool-runs/{tool}/{protein}/{ligand}/"
    else:
        search_str = f"tool-runs/{tool}/{protein}/"
    files = client.list_folder(search_str, recursive=True)
    files = list(files.keys())

    return files
