"""this module contains the Job class"""

import asyncio
from dataclasses import dataclass, field
import json
import time
from typing import Any, Callable, Optional

from IPython.display import HTML, clear_output, display
from ipywidgets import Output
import nest_asyncio

from deeporigin.tools import utils

# Enable nested event loops for Jupyter
nest_asyncio.apply()


@dataclass
class Job:
    """Job class"""

    name: str
    _ids: list[str]
    _viz_func: Optional[Callable[["Job"], Any]] = None
    _parse_func: Optional[Callable[["Job"], Any]] = None
    _progress_reports: list = field(default_factory=list)
    _status: list = field(default_factory=list)
    _task = None

    @classmethod
    def from_ids(cls, ids: list[str]) -> "Job":
        """Create a Job instance from a list of IDs.

        Args:
            ids: List of job IDs

        Returns:
            Job instance
        """
        return cls(name="job", _ids=ids)

    def sync(self):
        for i, job_id in enumerate(self._ids):
            # Skip if this job has already succeeded
            if i < len(self._status) and self._status[i] == "Succeeded":
                continue

            data = utils.get_status_and_progress(job_id)
            if data["progress"]:
                if i < len(self._progress_reports):
                    self._progress_reports[i] = json.loads(data["progress"])
                else:
                    self._progress_reports.append(json.loads(data["progress"]))
            elif i >= len(self._progress_reports):
                self._progress_reports.append(None)

            if data["status"]:
                if i < len(self._status):
                    self._status[i] = data["status"]
                else:
                    self._status.append(data["status"])
            elif i >= len(self._status):
                self._status.append(None)

    async def _poll_and_show(self):
        """Internal method that handles the polling and display loop"""
        while True:
            self.sync()
            with self._output:
                clear_output(wait=True)
                # First show the timestamp
                display(
                    HTML(
                        f"<div style='color: gray;'>Last updated: {time.strftime('%Y-%m-%d %H:%M:%S')}</div>"
                    )
                )
                # Then show the visualization
                self._viz_func(self)

            # Check if all jobs are either Succeeded or Failed
            if all(status in ["Succeeded", "Failed"] for status in self._status):
                break

            await asyncio.sleep(10)

    def show(self):
        """Show the job status and progress reports."""
        # Display the output widget in the current cell

        return self._viz_func(self)

    def start_polling(self):
        import asyncio
        import time
        from IPython.display import display, update_display, HTML

        # Create the initial display and store its display_id.
        initial_html = HTML("<div style='color: gray;'>Initializing...</div>")
        display_id = "timestamp_display"
        display(initial_html, display_id=display_id)

        async def update_timestamp():
            while True:
                current_time = time.strftime("%Y-%m-%d %H:%M:%S")
                new_html = (
                    f"<div style='color: gray;'>Last updated: {current_time}</div>"
                )

                new_html += self.show()

                update_display(HTML(new_html), display_id=display_id)

                await asyncio.sleep(5)

        # Schedule the task.
        self._task = asyncio.create_task(update_timestamp())

    def stop_polling(self):
        """Stop any ongoing background polling"""
        if self._task is not None:
            self._task.cancel()
            self._task = None

    def _repr_html_(self) -> str:
        """Return HTML representation for Jupyter notebooks.

        If a visualization function is set, it will be called.
        Otherwise, a basic HTML representation is returned.
        """
        if self._viz_func is not None:
            self.show()
            return ""
        else:
            return f"""
            <div style="border: 1px solid #ccc; padding: 10px; border-radius: 5px;">
                <h3>Job: {self.name}</h3>
                <p>Number of IDs: {len(self._ids)}</p>
                <p>Status: {self._status}</p>
                <p>Progress Reports: {self._progress_reports}</p>
            </div>
            """
