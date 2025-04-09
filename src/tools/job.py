"""this module contains the Job class"""

import asyncio
from dataclasses import dataclass, field
import json
import time
from typing import Any, Callable, Optional

from IPython.display import HTML, display, update_display
import nest_asyncio

from deeporigin.tools import utils

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
    _viz_func: Optional[Callable[["Job"], Any]] = None
    _parse_func: Optional[Callable[["Job"], Any]] = None
    _progress_reports: list = field(default_factory=list)
    _status: list = field(default_factory=list)
    _inputs: list = field(default_factory=list)
    _task = None
    _attributes: list = field(default_factory=list)
    _execution_ids: list = field(default_factory=list)
    _metadata: list = field(default_factory=list)

    @classmethod
    def from_ids(cls, ids: list[str]) -> "Job":
        """Create a Job instance from a list of IDs.

        Args:
            ids: List of job IDs to track.

        Returns:
            A new Job instance with the given IDs.
        """
        return cls(name="job", _ids=ids)

    def sync(self):
        """Synchronize the job status and progress reports.

        This method updates the internal state by fetching the latest status
        and progress reports for each job ID. It skips jobs that have already
        reached a terminal state (Succeeded or Failed).
        """

        # use
        results = utils.get_statuses_and_progress(self._ids)

        self._status = [result["status"] for result in results]
        self._progress_reports = [result["progress"] for result in results]
        self._execution_ids = [result["execution_id"] for result in results]
        self._inputs = [result["inputs"] for result in results]
        self._attributes = [result["attributes"] for result in results]
        self._metadata = [result["attributes"]["metadata"] for result in results]

    def show(self):
        """Display the current job status and progress reports.

        This method renders and displays the current state of all jobs
        using the visualization function if set, or a default HTML representation.
        """
        # Display the output widget in the current cell
        display(HTML(self._render_progress_html()))

    def _render_progress_html(self):
        """Render HTML representation of job progress.

        Returns:
            HTML string representing the current job status and progress.
        """
        return self._viz_func(self)

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
            display(HTML("<div style='color: gray;'>No active jobs to monitor</div>"))
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
                current_time = time.strftime("%Y-%m-%d %H:%M:%S")

                self.sync()
                new_html = self._render_progress_html()

                new_html += (
                    f"<div style='color: gray;'>Last updated: {current_time}</div>"
                )

                update_display(HTML(new_html), display_id=display_id)

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
            self._task.cancel()
            self._task = None

    def _repr_html_(self) -> str:
        """Return HTML representation for Jupyter notebooks.

        This method is called by Jupyter to display the job object in a notebook.
        It uses the visualization function if set, otherwise returns a basic
        HTML representation of the job's state.

        Returns:
            HTML string representing the job object.
        """
        if self._viz_func is not None:
            return self._render_progress_html()
        else:
            return f"""
            <div style="border: 1px solid #ccc; padding: 10px; border-radius: 5px;">
                <h3>Job: {self.name}</h3>
                <p>Number of IDs: {len(self._ids)}</p>
                <p>Status: {self._status}</p>
                <p>Progress Reports: {self._progress_reports}</p>
            </div>
            """

    def cancel(self):
        """Cancel all jobs being tracked by this instance.

        This method sends cancellation requests for all job IDs tracked by this instance
        using the utils.cancel_runs function.

        Returns:
            The result of the cancellation operation from utils.cancel_runs.
        """
        utils.cancel_runs(self._ids)
