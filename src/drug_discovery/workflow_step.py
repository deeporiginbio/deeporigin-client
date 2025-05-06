"""Base class for workflow steps like ABFE, RBFE, and Docking."""

from typing import Optional

from beartype import beartype

from deeporigin.tools.job import Job, get_dataframe
from deeporigin.utils.core import PrettyDict


class WorkflowStep:
    """Base class for workflow steps that handle jobs."""

    """
    If True, the jobs will be fused into a single job.
    This is useful for workflow steps that are run in parallel in batches,
    such as Docking.
    """
    _fuse_jobs: bool = False
    _tool_key: str = ""  # To be overridden by derived classes
    parent = None
    jobs: list[Job] | None = None

    def __init__(self, parent):
        self.parent = parent
        self._params = PrettyDict()

    def show_jobs(self, summary: Optional[bool] = None):
        """Show the jobs for this workflow step."""

        if self.jobs is None:
            print("No jobs to show")
            return

        if summary is None:
            if len(self.jobs) > 5:
                summary = True
            else:
                summary = False

        if summary:
            df = get_dataframe()
            df = df[df["tool_key"].str.contains(self._tool_key)]

            # filter by protein
            df = df[df["protein_id"] == self.parent.protein._do_id]
            return df
        else:
            for job in self.jobs:
                job.sync()
                job.show()

    @beartype
    def _make_jobs_from_ids(self, job_ids: list[str]) -> None:
        """Get a list of jobs for this workflow step. When _fuse_jobs is True, the jobs will be fused into a single job."""

        if len(job_ids) == 0:
            # nothing to do
            return

        if self._fuse_jobs:
            # fuse all job IDs into a single job
            job = Job.from_ids(job_ids)
            job._viz_func = self._render_progress
            job._name_func = self._name_job
            job.sync()
            self.jobs = [job]

        else:
            # make a new job for each job ID
            jobs = []
            for job_id in job_ids:
                job = Job.from_ids([job_id])
                job._viz_func = self._render_progress
                job._name_func = self._name_job
                job.sync()
                jobs.append(job)

            self.jobs = jobs

    @beartype
    def _render_progress(self, job: Job) -> str:
        """Render progress visualization for a job. To be implemented by subclasses."""
        raise NotImplementedError

    @beartype
    def _name_job(self, job: Job) -> str:
        """Generate a name for a job. To be implemented by subclasses."""
        raise NotImplementedError

    @beartype
    def _connect(self):
        """connect to deep origin API and perform any other setup necessary. This is meant to be implemented in the subclass"""
        raise NotImplementedError
