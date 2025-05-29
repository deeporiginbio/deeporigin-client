"""this module contains the Job class"""

import asyncio
from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
import time
from typing import Any, Callable, Optional
import uuid

from beartype import beartype
from dateutil import parser
import humanize
from IPython.display import HTML, display, update_display
from jinja2 import Environment, FileSystemLoader
import nest_asyncio
import pandas as pd

from deeporigin.platform import tools_api
from deeporigin.platform.utils import PlatformClients
from deeporigin.utils.core import elapsed_minutes

# Enable nested event loops for Jupyter
nest_asyncio.apply()


@dataclass
class Job:
    """A class representing a job that can be monitored and managed.

    This class provides functionality to track the status and progress of jobs,
    with support for real-time monitoring in Jupyter notebooks.

    Attributes:
        name: A string representing the name of the job.
        _ids: A list of job IDs being tracked.
        _viz_func: Optional function to customize job visualization.
        _parse_func: Optional function to parse job data.
        _progress_reports: List of progress reports for each job.
        _status: List of statuses for each job.
        _task: The asyncio task handling job monitoring.
    """

    name: str
    _ids: list[str]

    # functions
    _viz_func: Optional[Callable[["Job"], Any]] = None
    _parse_func: Optional[Callable[["Job"], Any]] = None
    _name_func: Optional[Callable[["Job"], Any]] = field(
        default_factory=lambda: lambda self: "Job"
    )

    _progress_reports: list = field(default_factory=list)
    _status: list = field(default_factory=list)
    _inputs: list = field(default_factory=list)
    _task = None
    _attributes: list = field(default_factory=list)
    _execution_ids: list = field(default_factory=list)
    _metadata: list = field(default_factory=list)

    # clients
    _platform_clients: Optional[PlatformClients] = None

    @classmethod
    def from_ids(
        cls,
        ids: list[str],
        *,
        _platform_clients: Optional[PlatformClients] = None,
    ) -> "Job":
        """Create a Job instance from a list of IDs.

        Args:
            ids: List of job IDs to track.

        Returns:
            A new Job instance with the given IDs.
        """
        return cls(
            name="job",
            _ids=ids,
            _platform_clients=_platform_clients,
        )

    @classmethod
    def from_id(
        cls,
        id: str,
        *,
        _platform_clients: Optional[PlatformClients] = None,
    ) -> "Job":
        """Create a Job instance from a single ID.

        Args:
            id: Job ID to track.

        Returns:
            A new Job instance with the given ID.
        """
        return cls(
            name="job",
            _ids=[id],
            _platform_clients=_platform_clients,
        )

    def sync(self):
        """Synchronize the job status and progress reports.

        This method updates the internal state by fetching the latest status
        and progress reports for each job ID. It skips jobs that have already
        reached a terminal state (Succeeded or Failed).
        """

        # get the org_friendly_id from the platform clients if it's not None
        org_friendly_id = getattr(self._platform_clients, "org_friendly_id", None)
        tools_client = getattr(self._platform_clients, "ToolsApi", None)

        # use
        results = tools_api.get_statuses_and_progress(
            self._ids,
            client=tools_client,
            org_friendly_id=org_friendly_id,
        )

        self._status = [result["status"] for result in results]
        self._progress_reports = [result["progress"] for result in results]
        self._execution_ids = [result["execution_id"] for result in results]
        self._inputs = [result["inputs"] for result in results]
        self._attributes = [result["attributes"] for result in results]
        self._metadata = [result["attributes"]["metadata"] for result in results]

    def _get_running_time(self) -> list:
        """Get the running time of the job.

        Returns:
            The running time of the job in minutes.
        """
        running_time = []
        for item in self._attributes:
            if item.completedAt is None or item.startedAt is None:
                running_time.append(0)
            else:
                running_time.append(elapsed_minutes(item.startedAt, item.completedAt))

        return running_time

    def _render_job_view(self):
        """Display the current job status and progress reports.

        This method renders and displays the current state of all jobs
        using the visualization function if set, or a default HTML representation.
        """
        # Get the template directory
        template_dir = Path(__file__).parent.parent / "templates"
        # Create Jinja2 environment with auto-escaping disabled
        # Note: Auto-escaping is disabled because the template needs to render HTML content
        # from _viz_func and properly formatted JSON data. The |safe filter is used
        # only for trusted content (JSON data and HTML from _viz_func).
        # All other template variables are properly escaped by the template itself.
        env = Environment(
            loader=FileSystemLoader(str(template_dir)),
            autoescape=False,  # Disabled for proper HTML and JSON rendering
        )

        from deeporigin.utils.notebook import get_notebook_environment

        if get_notebook_environment() == "marimo":
            template = env.get_template("job.html")
        else:
            template = env.get_template("job_jupyter.html")

        try:
            status_html = self._viz_func(self)
        except Exception as e:
            status_html = (
                f"No visualization function provided, or there was an error. Error: {e}"
            )

        try:
            card_title = self._name_func(self)
        except Exception:
            card_title = "No name function provided."

        started_at = []
        for item in self._attributes:
            if item.startedAt is None:
                started_at.append(None)
                continue
            dt = parser.isoparse(item.startedAt).astimezone(timezone.utc)

            # Compare to now (also in UTC)
            now = datetime.now(timezone.utc)
            started_at.append(humanize.naturaltime(now - dt))

        running_time = self._get_running_time()

        # Prepare template variables
        template_vars = {
            "status_html": status_html,
            "last_updated": time.strftime("%Y-%m-%d %H:%M:%S"),
            "outputs_json": json.dumps(
                [attribute.userOutputs for attribute in self._attributes], indent=2
            ),
            "inputs_json": json.dumps(
                [attribute.userInputs for attribute in self._attributes], indent=2
            ),
            "job_ids": self._ids,
            "execution_ids": self._execution_ids,
            "statuses": self._status,
            "started_at": started_at,
            "running_time": running_time,
            "card_title": card_title,
            "unique_id": str(uuid.uuid4()),
        }

        # Determine overall status based on priority: Failed > Cancelled > Succeeded
        if not self._status:
            template_vars["status"] = "Unknown"
        elif any(status == "Failed" for status in self._status):
            template_vars["status"] = "Failed"
        elif any(status == "Cancelled" for status in self._status):
            template_vars["status"] = "Cancelled"
        elif all(status == "Succeeded" for status in self._status):
            template_vars["status"] = "Succeeded"
        elif all(status == "Created" for status in self._status):
            template_vars["status"] = "Created"
        else:
            template_vars["status"] = "Running"

        # Try to parse progress reports as JSON, fall back to raw text if it fails
        try:
            # Try to parse each item in the list as JSON
            parsed_reports = []
            for report in self._progress_reports:
                parsed_reports.append(json.loads(str(report)))

            template_vars["raw_progress_json"] = json.dumps(parsed_reports, indent=2)
        except Exception:
            # If something goes wrong with the list processing, fall back to raw text
            template_vars["raw_progress_json"] = str(self._progress_reports)
            template_vars["raw_progress_json"].replace("\n", "<br>")

        # Render the template
        return template.render(**template_vars)

    def show(self):
        """Display the job view in a Jupyter notebook.

        This method renders the job view and displays it in a Jupyter notebook.
        """
        rendered_html = self._render_job_view()
        display(HTML(rendered_html))

    def watch(self):
        """Start monitoring job progress in real-time.

        This method initiates a background task that periodically updates
        and displays the job status. It will automatically stop when all
        jobs reach a terminal state (Succeeded or Failed). If there are no
        active jobs to monitor, it will display a message and show the current
        state once.
        """
        # Check if there are any active jobs (not Failed or Succeeded)
        if not any(status not in ["Failed", "Succeeded"] for status in self._status):
            display(
                HTML(
                    "<div style='color: gray;'>No active jobs to monitor. This display will not update.</div>"
                )
            )
            self.show()
            return

        # Stop any existing task before starting a new one
        self.stop_watching()

        initial_html = HTML("<div style='color: gray;'>Initializing...</div>")
        display_id = "timestamp_display"
        display(initial_html, display_id=display_id)

        async def update_progress_report():
            """Update and display job progress at regular intervals.

            This coroutine runs in the background, updating the display
            with the latest job status and progress every 5 seconds.
            It automatically stops when all jobs reach a terminal state.
            """
            while True:
                self.sync()
                html = self._render_job_view()

                update_display(HTML(html), display_id=display_id)

                # Check if all jobs are in terminal states
                if all(status in ["Failed", "Succeeded"] for status in self._status):
                    self.stop_watching()
                    break

                await asyncio.sleep(5)

        # Schedule the task.
        self._task = asyncio.create_task(update_progress_report())

    def stop_watching(self):
        """Stop the background monitoring task.

        This method safely cancels and cleans up any running monitoring task.
        It is called automatically when all jobs reach a terminal state,
        or can be called manually to stop monitoring.
        """
        if self._task is not None:
            self._task.cancel()  # note that this is not job.cancel() -- we're cancelling the asyncio task
            self._task = None

    def _repr_html_(self) -> str:
        """Return HTML representation for Jupyter notebooks.

        This method is called by Jupyter to display the job object in a notebook.
        It uses the visualization function if set, otherwise returns a basic
        HTML representation of the job's state.

        Returns:
            HTML string representing the job object.
        """

        return self._render_job_view()

    def cancel(self):
        """Cancel all jobs being tracked by this instance.

        This method sends cancellation requests for all job IDs tracked by this instance
        using the utils.cancel_runs function.

        Returns:
            The result of the cancellation operation from utils.cancel_runs.
        """

        org_friendly_id = getattr(self._platform_clients, "org_friendly_id", None)
        tools_client = getattr(self._platform_clients, "ToolsApi", None)

        tools_api.cancel_runs(
            self._ids,
            client=tools_client,
            org_friendly_id=org_friendly_id,
        )


@beartype
def get_dataframe(
    *,
    tool_key: Optional[str] = None,
    only_with_status: Optional[list[str] | set[str]] = None,
    include_metadata: bool = False,
    include_inputs: bool = False,
    include_outputs: bool = False,
    resolve_user_names: bool = False,
    _platform_clients: Optional[PlatformClients] = None,
) -> pd.DataFrame:
    """Get a dataframe of the job statuses and progress reports.

    Returns:
        A dataframe of the job statuses and progress reports.
    """

    if only_with_status is None:
        only_with_status = ["Succeeded", "Running", "Queued", "Failed", "Created"]

    if isinstance(only_with_status, set):
        only_with_status = list(only_with_status)

    _filter = {
        "status": {"$in": only_with_status},
        "metadata": {
            "$exists": True,
            "$ne": None,
        },
    }

    if tool_key is not None:
        _filter["tool"] = {
            "toolManifest": {
                "key": tool_key,
            },
        }

    if _platform_clients is None:
        from deeporigin.config import get_value

        org_friendly_id = get_value()["organization_id"]
    else:
        org_friendly_id = _platform_clients.org_friendly_id

    tools_client = getattr(_platform_clients, "ToolsApi", None)

    response = tools_api.get_tool_executions(
        org_friendly_id=org_friendly_id,
        filter=_filter,
        page_size=10000,
        client=tools_client,
    )
    jobs = response["data"]

    if resolve_user_names:
        from deeporigin.platform import organizations_api

        orgs_client = getattr(_platform_clients, "OrganizationsApi", None)

        users = organizations_api.get_organization_users(
            org_friendly_id=org_friendly_id,
            client=orgs_client,
        )

        # Create a mapping of user IDs to user names
        user_id_to_name = {user["id"]: user["attributes"]["name"] for user in users}

    def find_boolean_value(d, key):
        """Recursively search for a key in a nested dictionary and return its value as boolean."""
        if isinstance(d, dict):
            if key in d:
                return d[key] == 1
            for v in d.values():
                result = find_boolean_value(v, key)
                if result is not None:
                    return result
        elif isinstance(d, list):
            for item in d:
                result = find_boolean_value(item, key)
                if result is not None:
                    return result
        return False

    # Initialize lists to store data
    data = {
        "id": [],
        "createdAt": [],
        "executionId": [],
        "completedAt": [],
        "startedAt": [],
        "protein_id": [],
        "status": [],
        "tool_key": [],
        "tool_version": [],
        "user_name": [],
        "thread_pinning": [],
        "test_run": [],
        "run_duration_minutes": [],
        "n_ligands": [],
    }

    if include_metadata:
        data["metadata"] = []

    if include_inputs:
        data["user_inputs"] = []

    if include_outputs:
        data["user_outputs"] = []

    for job in jobs:
        attributes = job["attributes"]

        # Add basic fields
        data["id"].append(job["id"])
        data["createdAt"].append(attributes["createdAt"])
        data["executionId"].append(attributes["executionId"])
        data["completedAt"].append(attributes["completedAt"])
        data["startedAt"].append(attributes["startedAt"])
        data["status"].append(attributes["status"])
        data["tool_key"].append(attributes["tool"]["key"])
        data["tool_version"].append(attributes["tool"]["version"])

        if resolve_user_names:
            data["user_name"].append(
                user_id_to_name.get(attributes["user"]["id"], "Unknown")
            )
        else:
            data["user_name"].append(attributes["user"]["id"])

        # Check for thread_pinning and test_run in userInputs
        user_inputs = attributes.get("userInputs", {})
        data["thread_pinning"].append(find_boolean_value(user_inputs, "thread_pinning"))
        data["test_run"].append(find_boolean_value(user_inputs, "test_run"))

        if include_inputs:
            data["user_inputs"].append(user_inputs)

        if "smiles_list" in user_inputs:
            data["n_ligands"].append(len(user_inputs["smiles_list"]))
        else:
            data["n_ligands"].append(1)

        # Handle protein_id (may not exist or metadata may be None)
        metadata = attributes.get("metadata")
        this_protein_id = metadata.get("protein_id") if metadata else None

        if this_protein_id is None:
            this_protein_id = (
                attributes.get("userInputs", {}).get("protein", {}).get("rowId")
                or attributes.get("userInputs", {}).get("protein", {}).get("key")
                or "No ID"
            )

        data["protein_id"].append(this_protein_id)

        # Calculate run duration in minutes and round to nearest integer
        if attributes["completedAt"] and attributes["startedAt"]:
            start = parser.isoparse(attributes["startedAt"])
            end = parser.isoparse(attributes["completedAt"])
            duration = round((end - start).total_seconds() / 60)
            data["run_duration_minutes"].append(duration)
        else:
            data["run_duration_minutes"].append(None)

        if include_metadata:
            data["metadata"].append(metadata)

        if include_outputs:
            data["user_outputs"].append(attributes.get("userOutputs", {}))

    # Create DataFrame
    df = pd.DataFrame(data)

    # Convert datetime columns
    datetime_cols = ["createdAt", "completedAt", "startedAt"]
    for col in datetime_cols:
        df[col] = (
            pd.to_datetime(df[col], errors="coerce", utc=True)  # parse → tz-aware
            .dt.tz_localize(None)  # drop the UTC tz-info
            .astype("datetime64[us]")  # truncate to microseconds
        )

    return df
