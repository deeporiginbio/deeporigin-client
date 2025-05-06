"""This module encapsulates methods to run docking and show docking results on Deep Origin"""

import concurrent.futures
import math
import os
from typing import Optional

from beartype import beartype
from deeporigin_molstar import JupyterViewer
import more_itertools
import numpy as np
import pandas as pd
from tqdm import tqdm

from deeporigin.data_hub import api
from deeporigin.drug_discovery import chemistry as chem
from deeporigin.drug_discovery import utils
from deeporigin.drug_discovery.structures.pocket import Pocket
from deeporigin.drug_discovery.workflow_step import WorkflowStep
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools.job import Job
from deeporigin.utils.core import ensure_file_extension

Number = float | int


class Docking(WorkflowStep):
    """class to handle Docking-related tasks within the Complex class.

    Objects instantiated here are meant to be used within the Complex class."""

    def __init__(self, parent):
        super().__init__(parent)
        self._fuse_jobs = True

    @beartype
    def _name_job(self, job: Job) -> str:
        """Generate a name for a job."""

        unique_smiles = set()
        for inputs in job._inputs:
            unique_smiles.update(inputs["smiles_list"])
        num_ligands = len(unique_smiles)
        return f"Docking <code>{job._metadata[0]['protein_id']}</code> to {num_ligands} ligands."

    @classmethod
    @beartype
    def _render_progress(cls, job: Job) -> str:
        """Render progress visualization for a job."""

        data = job._progress_reports

        total_ligands = sum([len(inputs["smiles_list"]) for inputs in job._inputs])
        total_docked = 0
        total_failed = 0

        for item in data:
            if item is None:
                continue
            total_docked += item.count("ligand docked")
            total_failed += item.count("ligand failed")

        total_running_time = sum(job._get_running_time())
        speed = total_docked / total_running_time if total_running_time > 0 else 0

        from deeporigin.utils.notebook import render_progress_bar

        return render_progress_bar(
            completed=total_docked,
            total=total_ligands,
            failed=total_failed,
            title="Docking Progress",
            body_text=f"Average speed: {speed:.2f} dockings/minute",
        )

    def show_results(self):
        """show results of bulk Docking run in a table, rendering 2D structures of molecules"""

        df = self.get_results()

        if (df is None) or (len(df) == 0):
            # no results available yet
            return

        from IPython.display import HTML, display

        smiles_list = list(df["SMILES"])
        images = chem.smiles_list_to_base64_png_list(smiles_list)
        df["Structure"] = images
        df.drop("SMILES", axis=1, inplace=True)

        # Reorder columns to put Structure first
        cols = df.columns.tolist()
        cols.remove("Structure")
        df = df[["Structure"] + cols]

        display(HTML(df.to_html(escape=False)))

    def show_poses(self):
        """show docked ligands with protein in 3D"""

        file_paths = self.get_poses()

        if file_paths is None:
            # no results available yet
            return

        from deeporigin_molstar.src.viewers import DockingViewer

        docking_viewer = DockingViewer()
        html_content = docking_viewer.render_with_seperate_crystal(
            protein_data=str(self.parent.protein.file_path),
            protein_format="pdb",
            ligands_data=file_paths,
            ligand_format="sdf",
            paginate=True,
        )
        JupyterViewer.visualize(html_content)

    @beartype
    def _get_result_files(
        self,
        column_name: str,
        file_extension: str,
    ) -> list[str] | None:
        """Helper function to get result files of a specific type.

        Args:
            column_name: Name of the column containing file IDs (e.g. "OutputFile" or "ResultFile")
            file_extension: File extension to use (e.g. ".csv" or ".sdf")

        Returns:
            List of file paths to the result files
        """
        jobs = self._get_jobs(only_with_status=["Succeeded"])
        row_ids = [job.attributes.userOutputs.results_sdf.rowId for job in jobs]

        df = pd.DataFrame(
            api.get_dataframe(
                "Docking",
                return_type="dict",
                use_file_names=False,
            )
        )

        df = df[df["ID"].isin(row_ids)]

        if len(df) == 0:
            print("No results available yet")
            return None

        file_ids = df[column_name].tolist()

        existing_files = os.listdir(utils.DATA_DIRS["Docking"])
        existing_files = [
            "_file:" + file.replace(file_extension, "") for file in existing_files
        ]
        missing_files = list(set(file_ids) - set(existing_files))

        if len(missing_files) > 0:
            api.download_files(
                file_ids=missing_files,
                use_file_names=False,
                save_to_dir=utils.DATA_DIRS["Docking"],
            )

        file_paths = [
            os.path.join(utils.DATA_DIRS["Docking"], file_id.replace("_file:", ""))
            for file_id in file_ids
        ]

        file_paths = ensure_file_extension(
            file_paths=file_paths, extension=file_extension
        )

        return file_paths

    @beartype
    def get_results(self) -> pd.DataFrame | None:
        """return a list of paths to CSV files that contain the results from docking"""

        file_paths = self._get_result_files("OutputFile", ".csv")

        if file_paths is None:
            # no results, nothing to do
            return None

        for file in file_paths:
            from deeporigin.utils.core import fix_embedded_newlines_in_csv

            fix_embedded_newlines_in_csv(file)

        all_df = []
        for file in file_paths:
            df = pd.read_csv(file)
            all_df.append(df)

        df = pd.concat(all_df)
        return df

    @beartype
    def get_poses(self) -> list[str] | None:
        """return a list of paths to SDF files that contain the poses of all ligands after docking"""
        return self._get_result_files("ResultFile", ".sdf")

    def _connect(self):
        """fetch job IDs for this protein from DB"""

        jobs = self._get_jobs()
        job_ids = [job.id for job in jobs]
        self._make_jobs_from_ids(job_ids)

    def _get_jobs(
        self,
        *,
        pocket_center=None,
        box_size=None,
        only_with_status: Optional[list[str]] = None,
    ) -> list:
        """get all job IDs for this protein"""

        if only_with_status is None:
            only_with_status = ["Succeeded", "Running", "Queued", "Failed"]

        _filter = {
            "status": {"$in": only_with_status},
            "metadata": {
                "$exists": True,
                "$ne": None,
            },
            "tool": {
                "toolManifest": {
                    "key": "deeporigin.bulk-docking",
                },
            },
        }
        from deeporigin.config import get_value
        from deeporigin.platform import tools

        response = tools.get_tool_executions(
            org_friendly_id=get_value()["organization_id"],
            filter=_filter,
            page_size=10000,
        )
        jobs = response["data"]

        # remove jobs that have no metadata about protein id
        jobs = [job for job in jobs if "protein_id" in job.attributes.metadata.keys()]

        # remove jobs that don't match this protein
        jobs = [
            job
            for job in jobs
            if job.attributes.metadata.protein_id == self.parent.protein._do_id
        ]

        if pocket_center is not None:
            import numpy as np

            jobs = [
                job
                for job in jobs
                if bool(
                    np.all(
                        np.isclose(
                            pocket_center, job.attributes.userInputs.pocket_center
                        )
                    )
                )
            ]

        if box_size is not None:
            jobs = [
                job
                for job in jobs
                if bool(
                    np.all(np.isclose(box_size, job.attributes.userInputs.box_size))
                )
            ]

        # only keep jobs where at least one ligand in that job matches what we have in the current complex
        relevant_job_ids = []
        smiles_strings = [ligand.smiles for ligand in self.parent.ligands]

        for job in jobs:
            this_smiles = job.attributes.userInputs.smiles_list

            for smiles in this_smiles:
                if smiles in smiles_strings:
                    relevant_job_ids.append(job.id)
                    break

        jobs = [job for job in jobs if job.id in relevant_job_ids]

        return jobs

    @beartype
    def run(
        self,
        *,
        pocket: Optional[Pocket] = None,
        box_size: Optional[tuple[Number, Number, Number]] = None,
        pocket_center: Optional[tuple[Number, Number, Number]] = None,
        batch_size: Optional[int] = 32,
        n_workers: Optional[int] = None,
    ):
        """Run bulk docking on Deep Origin. Ligands will be split into batches based on the batch_size argument, and will run in parallel on Deep Origin clusters.

        Args:
            box_size (tuple[float, float, float]): box size
            pocket_center (tuple[float, float, float]): pocket center
            batch_size (int, optional): batch size. Defaults to 30.
            n_workers (int, optional): number of workers. Defaults to None.

        """

        if self.parent.protein._do_id is None:
            raise DeepOriginException(
                "Protein must be uploaded to Deep Origin before docking."
            )

        if batch_size is None and n_workers is None:
            raise DeepOriginException(
                "Either batch_size or n_workers must be specified."
            )
        elif batch_size is not None and n_workers is not None:
            print(
                "Both batch_size and n_workers are specified. Using n_workers to determine batch_size..."
            )

        if n_workers is not None:
            batch_size = math.ceil(len(self.parent.ligands) / n_workers)
            print(f"Using a batch size of {batch_size}")

        if pocket is None and box_size is None and pocket_center is None:
            raise DeepOriginException(
                "Specify a pocket, or a box size and pocket center."
            )

        if pocket is not None:
            box_size = float(2 * np.cbrt(pocket.props["volume"]))
            box_size = [box_size, box_size, box_size]
            pocket_center = pocket.get_center().tolist()

        smiles_strings = [ligand.smiles for ligand in self.parent.ligands]

        print(f"Docking {len(smiles_strings)} ligands...")

        jobs = self._get_jobs(pocket_center=pocket_center, box_size=box_size)

        already_docked_ligands = []
        for job in jobs:
            this_smiles = job.attributes.userInputs.smiles_list
            already_docked_ligands.extend(this_smiles)

        smiles_strings = set(smiles_strings) - set(already_docked_ligands)
        smiles_strings = list(smiles_strings)

        print(
            f"Docking {len(smiles_strings)} ligands, after filtering out already docked ligands..."
        )

        if len(smiles_strings) == 0:
            print("No new ligands to dock")
            return

        chunks = list(more_itertools.chunked(smiles_strings, batch_size))
        jobs = self._get_jobs()
        job_ids = [job.id for job in jobs]

        def process_chunk(chunk):
            params = dict(
                box_size=list(box_size),
                pocket_center=list(pocket_center),
                smiles_list=chunk,
            )

            return utils._start_tool_run(
                params=params,
                protein_id=self.parent.protein._do_id,
                database_columns=self.parent._db.proteins.cols
                + self.parent._db.docking.cols,
                complex_hash=self.parent._hash,
                tool="Docking",
            )

        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Submit all chunks to the executor
            future_to_chunk = {
                executor.submit(process_chunk, chunk): chunk for chunk in chunks
            }

            # Process results with progress bar
            for future in tqdm(
                concurrent.futures.as_completed(future_to_chunk),
                total=len(chunks),
                desc="Starting docking jobs",
            ):
                job_id = future.result()
                if job_id is not None:
                    job_ids.append(job_id)

        job = Job.from_ids(job_ids)
        job._viz_func = self._render_progress
        job._name_func = self._name_job
        job.sync()

        # for docking, we always have a single job
        self.jobs = [job]
