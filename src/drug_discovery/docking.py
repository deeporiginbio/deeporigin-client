"""This module encapsulates methods to run docking and show docking results on Deep Origin"""

import math
from typing import Optional

import more_itertools
import pandas as pd
from beartype import beartype
from deeporigin.data_hub import api
from deeporigin.drug_discovery import chemistry as chem
from deeporigin.drug_discovery import utils
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools.utils import query_run_statuses
from deeporigin.utils.core import PrettyDict, hash_strings
from deeporigin_molstar import DockingViewer, JupyterViewer

Number = float | int


class Docking:
    """class to handle Docking-related tasks within the Complex class.

    Objects instantiated here are meant to be used within the Complex class."""

    def __init__(self, parent):
        self.parent = parent
        self._params = PrettyDict()

    def get_results(self) -> pd.DataFrame:
        """Get docking results from Deep Origin"""

        # to do -- some way to make sure that we handle failed runs, complete runs, etc.
        # status = self.get_status("Docking")

        df1 = self.parent.get_csv_results_for("Docking")

        df2 = chem.ligands_to_dataframe(self.parent.ligands)
        df2["SMILES"] = df2["Ligand"]
        df2.drop(columns=["Ligand"], inplace=True)

        df = pd.merge(
            df1,
            df2,
            on="SMILES",
            how="inner",
        )
        return df

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
        box_size: tuple[Number, Number, Number],
        pocket_center: tuple[Number, Number, Number],
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

        smiles_strings = [ligand.smiles_string for ligand in self.parent.ligands]

        chunks = list(more_itertools.chunked(smiles_strings, batch_size))

        params = dict(
            box_size=list(box_size),
            pocket_center=list(pocket_center),
        )

        database_columns = self.parent._db.proteins.cols + self.parent._db.docking.cols

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

            self.parent._job_ids["Docking"].append(job_id)
