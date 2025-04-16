"""This module encapsulates methods to run docking and show docking results on Deep Origin"""

import math
from typing import Optional

from beartype import beartype
from deeporigin_molstar import DockingViewer, JupyterViewer
import more_itertools
import numpy as np
import pandas as pd

from deeporigin.data_hub import api
from deeporigin.drug_discovery import chemistry as chem
from deeporigin.drug_discovery import utils
from deeporigin.drug_discovery.structures.ligand import ligands_to_dataframe
from deeporigin.drug_discovery.structures.pocket import Pocket
from deeporigin.drug_discovery.workflow_step import WorkflowStep
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools.job import Job
from deeporigin.tools.utils import get_statuses_and_progress, query_run_statuses
from deeporigin.utils.core import hash_strings

Number = float | int


class Docking(WorkflowStep):
    """class to handle Docking-related tasks within the Complex class.

    Objects instantiated here are meant to be used within the Complex class."""

    def __init__(self, parent):
        super().__init__(parent)
        self._fuse_jobs = True

    @beartype
    def _render_progress(self, job: Job) -> str:
        """Render progress visualization for a job."""
        # TODO: Implement Docking-specific progress visualization
        return "Docking Progress Visualization"

    @beartype
    def _name_job(self, job: Job) -> str:
        """Generate a name for a job."""

        num_ligands = sum(len(inputs["smiles_list"]) for inputs in job._inputs)
        return f"Docking run. Docking protein: <code>{job._metadata[0]['protein_id']}</code> to {num_ligands} ligands."

    def get_results(self) -> pd.DataFrame:
        """Get docking results from Deep Origin"""

        df1 = self.parent.get_csv_results_for("Docking")

        df2 = ligands_to_dataframe(self.parent.ligands)
        df2["SMILES"] = df2["Ligand"]
        df2.drop(columns=["Ligand"], inplace=True)

        df = pd.merge(
            df1,
            df2,
            on="SMILES",
            how="inner",
        )
        return df

    def show_progress(self):
        """show progress of bulk Docking run"""

        data = get_statuses_and_progress(self._job_ids)

        total_docked = 0

        total_ligands = 0
        for item in data:
            progress = item["progress"]
            if progress is None:
                continue

            batch_docked = _parse_progress(progress)
            total_docked += batch_docked
            total_ligands += len(item["inputs"].smiles_list)

        if total_ligands == 0:
            print("Cannot show progress yet. Jobs are yet to start.")
            return

        from deeporigin.utils.notebook import show_progress_bar

        show_progress_bar(
            completed=total_docked,
            total=total_ligands,
            title="Docking Progress",
        )

    def show_results(self):
        """show results of bulk Docking run in a table, rendering 2D structures of molecules"""

        df = self.get_results()

        if len(df) == 0:
            print("No results found")
            return

        from IPython.display import HTML, display

        smiles_list = list(df["SMILES"])
        images = chem.smiles_list_to_base64_png_list(smiles_list)
        df["Structure"] = images
        df.drop("SMILES", axis=1, inplace=True)
        display(HTML(df.to_html(escape=False)))

    def show_poses(self):
        """show docked ligands with protein in 3D"""

        files = self.parent.get_result_files_for(tool="Docking")

        sdf_file = chem.merge_sdf_files(files)

        # use the Deep Origin Docking Viewer
        # to view the docking results
        docking_viewer = DockingViewer()
        html_content = docking_viewer.render_with_seperate_crystal(
            protein_data=str(self.parent.protein.file),
            protein_format="pdb",
            ligands_data=[sdf_file],
            ligand_format="sdf",
        )

        JupyterViewer.visualize(html_content)

    def get_poses(self, output_sdf_file: str) -> None:
        """generate a single SDF file containing all the poses of all ligands docked to the protein

        Args:
            output_sdf_file (str): path to output SDF file. All poses will be written to a SDF file in this location.

        """

        files = self.parent.get_result_files_for(tool="Docking")

        chem.merge_sdf_files(files, output_sdf_file)

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
            raise DeepOriginException(
                "Either batch_size or n_workers must be specified, but not both."
            )

        if n_workers is not None:
            batch_size = math.ceil(len(self.parent.ligands) / n_workers)

        if pocket is None and box_size is None and pocket_center is None:
            raise DeepOriginException(
                "Specify a pocket, or a box size and pocket center."
            )

        if pocket is not None:
            box_size = float(2 * np.cbrt(pocket.props["volume"]))
            box_size = [box_size, box_size, box_size]
            pocket_center = pocket.get_center().tolist()

        df = pd.DataFrame(
            api.get_dataframe(
                utils.DB_DOCKING,
                return_type="dict",
                use_file_names=False,
            )
        )

        df = df[self.parent._hash == df[utils.COL_COMPLEX_HASH]]

        # remove failed runs
        statuses = query_run_statuses(df["JobID"].tolist())
        df["Status"] = df["JobID"].replace(statuses)
        df = df["Failed" != df["Status"]]

        smiles_strings = [ligand.smiles for ligand in self.parent.ligands]

        chunks = list(more_itertools.chunked(smiles_strings, batch_size))

        params = dict(
            box_size=list(box_size),
            pocket_center=list(pocket_center),
        )

        database_columns = self.parent._db.proteins.cols + self.parent._db.docking.cols

        job_ids = []

        for chunk in chunks:
            params["smiles_list"] = chunk
            smiles_hash = hash_strings(params["smiles_list"])

            if smiles_hash in df[utils.COL_SMILES_HASH].tolist():
                print("Skipping this tranche because this has already been run...")
                continue

            job_id = utils._start_tool_run(
                params=params,
                protein_id=self.parent.protein._do_id,
                database_columns=database_columns,
                complex_hash=self.parent._hash,
                tool="Docking",
            )
            job_ids.append(job_id)

        if self.jobs is None:
            self.jobs = []

        job = Job.from_ids(job_ids)

        job._viz_func = self._render_progress
        job._name_func = self._name_job

        self.jobs.append(job)


@beartype
def _parse_progress(txt: str) -> int:
    """Parse Docking progress from raw progress text"""

    txt = txt.split("\n")
    num_docked_ligands = len(txt) - 1

    return num_docked_ligands
