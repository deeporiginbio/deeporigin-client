"""This module encapsulates methods to run docking and show docking results on Deep Origin"""

import concurrent.futures
import hashlib
import math
import os
from pathlib import Path
from typing import Optional

from beartype import beartype
from deeporigin_molstar import JupyterViewer
import more_itertools
import numpy as np
import pandas as pd
from tqdm import tqdm

from deeporigin.drug_discovery import chemistry as chem
from deeporigin.drug_discovery import utils
from deeporigin.drug_discovery.constants import tool_mapper
from deeporigin.drug_discovery.structures.pocket import Pocket
from deeporigin.drug_discovery.workflow_step import WorkflowStep
from deeporigin.exceptions import DeepOriginException
from deeporigin.platform import files_api, tools_api
from deeporigin.tools.job import Job, get_dataframe

Number = float | int
LOCAL_BASE = Path.home() / ".deeporigin"


class Docking(WorkflowStep):
    """class to handle Docking-related tasks within the Complex class.

    Objects instantiated here are meant to be used within the Complex class."""

    """tool version to use for Docking"""
    tool_version = "0.4.0"
    _tool_key = "deeporigin.bulk-docking"  # Tool key for Docking jobs

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

        protein_file = os.path.basename(job._inputs[0]["protein"]["key"])

        return f"Docking <code>{protein_file}</code> to {num_ligands} ligands."

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

        file_paths = self.get_results(file_type="sdf")

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
        raise NotImplementedError("Not implemented yet")

    @beartype
    def get_results(self, *, file_type: str = "csv") -> pd.DataFrame | None | list[str]:
        """return a list of paths to CSV files that contain the results from docking"""

        files_client = getattr(self.parent._platform_clients, "FilesApi", None)

        files = utils.find_files_on_ufa(
            tool="Docking",
            protein=self.parent.protein.file_path.name,
            client=files_client,
        )

        if file_type == "csv":
            results_files = [file for file in files if file.endswith("/results.csv")]
        elif file_type == "sdf":
            results_files = [file for file in files if file.endswith("/results.sdf")]
        else:
            raise ValueError(f"Invalid file type: {file_type}")

        if len(results_files) == 0:
            print("No Docking results found for this protein.")
            return None

        files_api.download_files(
            results_files,
            client=files_client,
        )

        all_df = []

        home_dir = os.path.expanduser("~")

        local_paths = [
            os.path.join(home_dir, ".deeporigin", file) for file in results_files
        ]

        if file_type == "csv":
            for local_path in local_paths:
                from deeporigin.utils.core import fix_embedded_newlines_in_csv

                fix_embedded_newlines_in_csv(local_path)

                df = pd.read_csv(local_path)
                all_df.append(df)

            df = pd.concat(all_df)
            return df
        else:
            return local_paths

    @beartype
    def _get_jobs(
        self,
        *,
        pocket_center=None,
        box_size=None,
        only_with_status: Optional[list[str]] = None,
    ):
        """search for all jobs that match this protein and ligands in the Job DB, and return a dataframe of the results"""

        df = get_dataframe(
            tool_key=tool_mapper["Docking"],
            only_with_status=tools_api.NON_FAILED_STATES,
            include_metadata=True,
            include_inputs=True,
            include_outputs=True,
        )

        if pocket_center is not None:
            # Filter df rows where pocket_center matches row["user_inputs"]["pocket_center"]
            mask = df["user_inputs"].apply(
                lambda x: bool(np.all(np.isclose(pocket_center, x["pocket_center"])))
            )
            df = df[mask]

        if box_size is not None:
            # Filter df rows where box_size matches row["user_inputs"]["box_size"]
            mask = df["user_inputs"].apply(
                lambda x: bool(np.all(np.isclose(box_size, x["box_size"])))
            )
            df = df[mask]

        # filter to only keep jobs that match this protein.
        protein_basename = os.path.basename(self.parent.protein.file_path)
        mask = df["metadata"].apply(lambda x: x["protein_file"] == protein_basename)
        df = df[mask]

        # only keep jobs where at least one ligand in that job matches what we have in the current complex
        smiles_strings = [ligand.smiles for ligand in self.parent.ligands]
        mask = df["user_inputs"].apply(
            lambda x: any(s in smiles_strings for s in x["smiles_list"])
        )
        df = df[mask]
        return df

    @beartype
    def run(
        self,
        *,
        pocket: Optional[Pocket] = None,
        box_size: Optional[tuple[Number, Number, Number]] = None,
        pocket_center: Optional[tuple[Number, Number, Number]] = None,
        batch_size: Optional[int] = 32,
        n_workers: Optional[int] = None,
        _output_dir_path: Optional[str] = None,
        use_parallel: bool = True,
    ):
        """Run bulk docking on Deep Origin. Ligands will be split into batches based on the batch_size argument, and will run in parallel on Deep Origin clusters.

        Args:
            box_size (tuple[float, float, float]): box size
            pocket_center (tuple[float, float, float]): pocket center
            batch_size (int, optional): batch size. Defaults to 30.
            n_workers (int, optional): number of workers. Defaults to None.
            use_parallel (bool, optional): whether to run jobs in parallel. Defaults to True.
        """

        protein_basename = os.path.basename(self.parent.protein.file_path)

        if _output_dir_path is None:
            _output_dir_path = "tool-runs/Docking/" + protein_basename + "/"

        self.parent._sync_protein_and_ligands()

        metadata = dict(protein_file=protein_basename)

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

        df = self._get_jobs(pocket_center=pocket_center, box_size=box_size)

        already_docked_ligands = []
        for _, row in df.iterrows():
            this_smiles = row["user_inputs"]["smiles_list"]
            already_docked_ligands.extend(this_smiles)

        smiles_strings = set(smiles_strings) - set(already_docked_ligands)
        smiles_strings = list(smiles_strings)

        print(
            f"Docking {len(smiles_strings)} ligands, after filtering out already docked ligands..."
        )

        job_ids = df["id"].tolist()
        chunks = list(more_itertools.chunked(smiles_strings, batch_size))

        def process_chunk(chunk):
            params = dict(
                box_size=list(box_size),
                pocket_center=list(pocket_center),
                smiles_list=chunk,
            )

            # Create a stable hash for the chunk
            chunk_str = ",".join(chunk)
            chunk_hash = hashlib.md5(chunk_str.encode("utf-8")).hexdigest()
            this_output_dir_path = os.path.join(_output_dir_path, chunk_hash) + "/"

            return utils._start_tool_run(
                params=params,
                metadata=metadata,
                protein_path=self.parent.protein._remote_path,
                tool="Docking",
                tool_version=self.tool_version,
                _platform_clients=self.parent._platform_clients,
                _output_dir_path=this_output_dir_path,
            )

        if len(smiles_strings) > 0:
            if use_parallel:
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
            else:
                for chunk in tqdm(
                    chunks, total=len(chunks), desc="Starting docking jobs"
                ):
                    job_id = process_chunk(chunk)
                    if job_id is not None:
                        job_ids.append(job_id)
        else:
            print("No new ligands to dock")

        job = Job.from_ids(job_ids)
        job._viz_func = self._render_progress
        job._name_func = self._name_job
        job.sync()

        # for docking, we always have a single job
        self.jobs = [job]

        return job
