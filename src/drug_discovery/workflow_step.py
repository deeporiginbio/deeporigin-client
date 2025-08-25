"""Base class for workflow steps like ABFE, RBFE, and Docking."""

from beartype import beartype
import pandas as pd

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

    def get_jobs_df(
        self,
        *,
        include_outputs: bool = True,
    ) -> pd.DataFrame:
        """Get the jobs for this workflow step as a dataframe"""
        df = get_dataframe(
            _platform_clients=self.parent._platform_clients,
            resolve_user_names=False,
            include_metadata=True,
            include_outputs=include_outputs,
        )

        if len(df) == 0:
            return df

        # filter by tool key
        df = df[df["tool_key"].str.contains(self._tool_key)]

        # filter by protein file
        df = df[
            df["metadata"].apply(
                lambda x: isinstance(x, dict)
                and x.get("protein_hash") == self.parent.protein.to_hash()
            )
        ]

        df = df.reset_index(drop=True)

        return df

    @beartype
    def get_jobs(self) -> None:
        """Get the jobs for this workflow step and save to self.jobs"""
        df = self.get_jobs_df()

        job_ids = df["id"].tolist()

        if self._fuse_jobs:
            self.jobs = [
                Job.from_ids(
                    job_ids,
                    _platform_clients=self.parent._platform_clients,
                )
            ]
        else:
            self.jobs = [
                Job.from_id(
                    job_id,
                    _platform_clients=self.parent._platform_clients,
                )
                for job_id in job_ids
            ]

        for job in self.jobs:
            job.sync()

    @beartype
    def _make_jobs_from_ids(self, job_ids: list[str]) -> None:
        """Get a list of jobs for this workflow step. When _fuse_jobs is True, the jobs will be fused into a single job."""

        if len(job_ids) == 0:
            # nothing to do
            return

        if self._fuse_jobs:
            # fuse all job IDs into a single job
            job = Job.from_ids(job_ids)
            self.jobs = [job]

        else:
            # make a new job for each job ID
            jobs = []
            for job_id in job_ids:
                job = Job.from_ids([job_id])
                jobs.append(job)

            self.jobs = jobs
