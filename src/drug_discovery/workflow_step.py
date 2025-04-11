"""Base class for workflow steps like ABFE, RBFE, and Docking."""

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

    @property
    def jobs(self) -> list[Job]:
        """Get a list of jobs for this workflow step. When _fuse_jobs is True, the jobs will be fused into a single job."""

        if self._fuse_jobs and self._job_ids:
            job = Job.from_ids(self._job_ids)
            job._viz_func = self._render_progress
            job._name_func = self._name_job
            job.sync()
            return [job]

        jobs = []
        for job_id in self._job_ids:
            job = Job.from_ids([job_id])
            job._viz_func = self._render_progress
            job._name_func = self._name_job
            job.sync()
            jobs.append(job)
        return jobs

    def _render_progress(self, job: Job) -> str:
        """Render progress visualization for a job. To be implemented by subclasses."""
        raise NotImplementedError

    def _name_job(self, job: Job) -> str:
        """Generate a name for a job. To be implemented by subclasses."""
        raise NotImplementedError

    def _get_jobs(self, job_ids: list[str]) -> None:
        """Generate job objects from job IDs.

        Args:
            job_ids: List of job IDs to store
        """

        # the default behavior is to create a job object for each job ID
        # for other behaviors, subclasses can override this method
        self.jobs = [Job.from_ids([job_id]) for job_id in job_ids]
