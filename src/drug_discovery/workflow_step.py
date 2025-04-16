"""Base class for workflow steps like ABFE, RBFE, and Docking."""

from beartype import beartype

from deeporigin.tools.job import Job
from deeporigin.utils.core import PrettyDict


class WorkflowStep:
    """Base class for workflow steps that handle jobs."""

    """
    If True, the jobs will be fused into a single job.
    This is useful for workflow steps that are run in parallel in batches,
    such as Docking.
    """
    _fuse_jobs: bool = False
    parent = None
    jobs: list[Job] = None

    def __init__(self, parent):
        self.parent = parent
        self._params = PrettyDict()

    def show_jobs(self):
        """Show the jobs for this workflow step."""

        for job in self.jobs:
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
