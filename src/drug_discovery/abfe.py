"""This module encapsulates methods to run ABFE and show ABFE results on Deep Origin.

The ABFE object instantiated here is contained in the Complex class is meant to be used within that class."""

import os
from pathlib import Path
from typing import Literal, Optional

from beartype import beartype
import pandas as pd

from deeporigin.drug_discovery import utils
from deeporigin.drug_discovery.constants import tool_mapper
from deeporigin.drug_discovery.structures.ligand import Ligand, LigandSet
from deeporigin.drug_discovery.workflow_step import WorkflowStep
from deeporigin.exceptions import DeepOriginException
from deeporigin.platform import file_api
from deeporigin.tools.job import Job, get_dataframe
from deeporigin.utils.notebook import get_notebook_environment

LOCAL_BASE = Path.home() / ".deeporigin"


class ABFE(WorkflowStep):
    """class to handle ABFE-related tasks within the Complex class.

    Objects instantiated here are meant to be used within the Complex class."""

    """tool version to use for ABFE"""
    tool_version = "0.2.8"
    _tool_key = tool_mapper["ABFE"]

    _max_atom_count: int = 100_000

    def __init__(self, parent):
        super().__init__(parent)

        self._params.end_to_end = utils._load_params("abfe_end_to_end")

    def get_results(self) -> pd.DataFrame | None:
        """get ABFE results and return in a dataframe.

        This method returns a dataframe showing the results of ABFE runs associated with this simulation session. The ligand file name and ΔG are shown, together with user-supplied properties"""

        df = self.get_jobs(include_outputs=True)

        results_files = []

        for _, row in df.iterrows():
            file_path = row["user_outputs"]["abfe_results_summary"]["key"]
            results_files.append(file_path)

        return results_files

        files_client = getattr(self.parent._platform_clients, "FilesApi", None)

        files = file_api.list_files_in_dir(
            file_path=f"tool-runs/ABFE/{self.parent.protein.to_hash()}.pdb/",
            client=files_client,
        )

        results_files = [file for file in files if file.endswith("/results.csv")]
        results_files = dict.fromkeys(results_files, None)

        if len(results_files) == 0:
            print("No ABFE results found for this protein.")
            return None

        results_files = file_api.download_files(
            results_files,
            client=files_client,
        )

        # read all the CSV files using pandas and
        # set Ligand1 column to ligand name (parent dir of results.csv)
        dfs = []
        for file in results_files:
            df = pd.read_csv(LOCAL_BASE / file)
            ligand_name = os.path.basename(os.path.dirname(file))
            df["File"] = ligand_name
            dfs.append(df)
        df1 = pd.concat(dfs)

        df1.drop(columns=["Ligand1", "Ligand2"], inplace=True)

        df2 = self.parent.ligands.to_dataframe()
        df2["SMILES"] = df2["Ligand"]
        df2.drop(columns=["Ligand", "initial_smiles"], inplace=True)

        return df1, df2

        df = pd.merge(
            df1,
            df2,
            on="File",
            how="inner",
            validate="one_to_one",
        )

        # rename the File column to Ligand
        df.rename(columns={"File": "Ligand"}, inplace=True)

        df["Protein"] = self.parent.protein.file_path.name

        return df

    def show_results(self):
        """Show ABFE results in a dataframe.

        This method returns a dataframe showing the results of ABFE runs associated with this simulation session. The ligand file name, 2-D structure, and ΔG are shown."""

        df = self.get_results()

        if df is None or len(df) == 0:
            return

        from rdkit.Chem import PandasTools

        PandasTools.AddMoleculeColumnToFrame(df, smilesCol="SMILES", molCol="Structure")
        PandasTools.RenderImagesInAllDataFrames()

        # show structure first
        new_order = ["Structure"] + [col for col in df.columns if col != "Structure"]

        # re‑index DataFrame
        df = df[new_order]

        if get_notebook_environment() == "marimo":
            import marimo as mo

            return mo.plain(df)

        else:
            return df

    def get_jobs(
        self,
        *,
        include_metadata: bool = False,
        include_outputs: bool = False,
    ):
        """get jobs for this workflow step"""
        df = super().get_jobs_df(include_outputs=include_outputs)

        ligand_hashes = [ligand.to_hash() for ligand in self.parent.ligands]

        # filter df by ligand_hash
        df = df[df["metadata"].apply(lambda d: d.get("ligand_hash") in ligand_hashes)]

        # make a new column called ligand_smiles using the metadata column
        df["ligand_smiles"] = df["metadata"].apply(
            lambda d: d.get("ligand_smiles") if isinstance(d, dict) else None
        )

        # make a new column called protein_file using the metadata column
        df["protein_name"] = df["metadata"].apply(
            lambda d: d.get("protein_name") if isinstance(d, dict) else None
        )

        # make a new column called ligand_file using the metadata column
        df["ligand_name"] = df["metadata"].apply(
            lambda d: d.get("ligand_name") if isinstance(d, dict) else None
        )

        if not include_metadata:
            df.drop(columns=["metadata"], inplace=True)

        return df

    def show_jobs(self):
        """show jobs for this workflow step"""
        df = self.get_jobs()

        return utils.render_smiles_in_dataframe(df, smiles_col="ligand_smiles")

    @beartype
    def set_test_run(self, value: int = 1):
        """set test_run parameter in abfe parameters"""

        utils._set_test_run(self._params.end_to_end, value)

    @beartype
    def _get_ligands_to_run(
        self,
        *,
        ligands: list[Ligand] | LigandSet,
        re_run: bool,
    ) -> list[Ligand]:
        """Helper method to determine which ligands need to be run based on already run jobs and re_run flag."""

        if isinstance(ligands, LigandSet):
            ligands = ligands.ligands

        if re_run:
            # we're re-running, so we need to re-run all ligands
            return ligands

        df = get_dataframe(
            tool_key=tool_mapper["ABFE"],
            only_with_status=["Succeeded", "Running", "Queued", "Created"],
            include_metadata=True,
            resolve_user_names=False,
            _platform_clients=self.parent._platform_clients,
        )

        # Build set of ligand names that have already been run
        if len(df) > 0:
            ligand_hashes_already_run = {
                ligand_hash
                for ligand_hash in df["metadata"].apply(
                    lambda d: d.get("ligand_hash") if isinstance(d, dict) else None
                )
                if isinstance(ligand_hash, str) and ligand_hash
            }
        else:
            ligand_hashes_already_run = set()

        print(ligand_hashes_already_run)

        # no re-run, remove already run ligands
        ligands_to_run = [
            ligand
            for ligand in ligands
            if ligand.to_hash() not in ligand_hashes_already_run
        ]
        return ligands_to_run

    @beartype
    def run(
        self,
        *,
        ligands: Optional[list[Ligand] | LigandSet] = None,
        ligand: Optional[Ligand] = None,
        re_run: bool = False,
        _output_dir_path: Optional[str] = None,
    ) -> list[Job] | None:
        """Method to run an end-to-end ABFE run.

        Args:
            ligands: List of ligand to run. Defaults to None. When None, all ligands in the object will be run. To view a list of valid ligands, use the `.show_ligands()` method"""

        if ligands is None and ligand is None:
            ligands = self.parent.ligands
        elif ligands is None:
            ligands = [ligand]

        if isinstance(ligands, LigandSet):
            ligands = ligands.ligands

        # check that there is a prepared system for each ligand
        for ligand in ligands:
            if ligand.to_hash() not in self.parent._prepared_systems:
                raise DeepOriginException(
                    f"Complex with Ligand {ligand.name} is not prepared. Please prepare the system using the `prepare` method of Complex."
                ) from None

        # check that for every prepared system, the number of atoms is less than the max atom count
        for ligand_name, prepared_system in self.parent._prepared_systems.items():
            from deeporigin.drug_discovery.external_tools.utils import (
                count_atoms_in_pdb_file,
            )

            num_atoms = count_atoms_in_pdb_file(prepared_system)

            if num_atoms > self._max_atom_count:
                raise ValueError(
                    f"System with {ligand_name} has too many atoms. It has {num_atoms} atoms, but the maximum allowed is {self._max_atom_count}."
                )

        self.parent._sync_protein_and_ligands()

        ligands_to_run = self._get_ligands_to_run(ligands=ligands, re_run=re_run)

        if len(ligands_to_run) == 0:
            print(
                "All requested ligands have already been run, or are queued to run. To re-run, set re_run=True"
            )
            return

        if self.jobs is None:
            self.jobs = []

        jobs_for_this_run = []

        # TODO -- parallelize this
        for ligand in ligands_to_run:
            metadata = {
                "protein_hash": self.parent.protein.to_hash(),
                "ligand_hash": ligand.to_hash(),
                "ligand_smiles": ligand.smiles,
                "protein_name": self.parent.protein.name,
                "ligand_name": ligand.name,
            }

            job_id = utils._start_tool_run(
                metadata=metadata,
                ligand1_path=ligand._remote_path,
                protein_path=self.parent.protein._remote_path,
                params=self._params.end_to_end,
                tool="ABFE",
                tool_version=self.tool_version,
                _platform_clients=self.parent._platform_clients,
                _output_dir_path=_output_dir_path,
            )

            job = Job.from_id(job_id, _platform_clients=self.parent._platform_clients)

            self.jobs.append(job)
            jobs_for_this_run.append(job)

        return jobs_for_this_run

    @beartype
    def show_trajectory(
        self,
        *,
        ligand: Ligand,
        step: Literal["md", "binding"],
        window: int = 1,
    ):
        """Show the system trajectory FEP run.

        Args:
            ligand: The ligand to show the trajectory for.
            step (Literal["md", "abfe"]): The step to show the trajectory for.
            window (int, optional): The window number to show the trajectory for. Defaults to 1.
        """

        remote_base = Path(
            f"tool-runs/ABFE/{self.parent.protein.file_path.name}/{Path(ligand.file_path).name}"
        )

        remote_pdb_file = (
            remote_base / "output/protein/ligand/systems/complex/complex.pdb"
        )
        files_to_download = [remote_pdb_file]

        files_client = getattr(self.parent._platform_clients, "FilesApi", None)

        if step == "binding":
            # Check for valid windows

            # figure out valid windows
            files = utils.find_files_on_ufa(
                tool="ABFE",
                protein=self.parent.protein.file_path.name,
                ligand=Path(ligand.file_path).name,
                client=files_client,
            )
            xtc_files = [
                file
                for file in files
                if file.endswith("Prod_1/_allatom_trajectory_40ps.xtc")
                and "binding/binding" in file
            ]

            import re

            valid_windows = [
                int(re.search(r"window_(\d+)", path).group(1)) for path in xtc_files
            ]

            if window not in valid_windows:
                raise DeepOriginException(
                    f"Invalid window number: {window}. Valid windows are: {sorted(valid_windows)}"
                ) from None

            remote_xtc_file = (
                remote_base
                / f"output/protein/ligand/binding/binding/window_{window}/Prod_1/_allatom_trajectory_40ps.xtc"
            )

        else:
            remote_xtc_file = (
                remote_base
                / "output/protein/ligand/simple_md/simple_md/prod/_allatom_trajectory_40ps.xtc"
            )

        files_to_download.append(remote_xtc_file)
        files_to_download = dict.fromkeys(map(str, files_to_download), None)

        file_api.download_files(files_to_download, client=files_client)

        from deeporigin_molstar.src.viewers import ProteinViewer

        protein_viewer = ProteinViewer(
            data=str(LOCAL_BASE / remote_pdb_file), format="pdb"
        )
        html_content = protein_viewer.render_trajectory(
            str(LOCAL_BASE / remote_xtc_file)
        )

        from deeporigin_molstar import JupyterViewer

        JupyterViewer.visualize(html_content)
