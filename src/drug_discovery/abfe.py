"""This module encapsulates methods to run ABFE and show ABFE results on Deep Origin.

The ABFE object instantiated here is contained in the Complex class is meant to be used within that class."""

import json
import os
from pathlib import Path
from typing import Literal, Optional

from beartype import beartype
import pandas as pd

from deeporigin.drug_discovery import utils
from deeporigin.drug_discovery.constants import tool_mapper
from deeporigin.drug_discovery.structures.ligand import Ligand, ligands_to_dataframe
from deeporigin.drug_discovery.workflow_step import WorkflowStep
from deeporigin.exceptions import DeepOriginException
from deeporigin.platform import files_api
from deeporigin.tools.job import Job, get_dataframe
from deeporigin.utils.notebook import get_notebook_environment

LOCAL_BASE = Path.home() / ".deeporigin"


class ABFE(WorkflowStep):
    """class to handle ABFE-related tasks within the Complex class.

    Objects instantiated here are meant to be used within the Complex class."""

    """tool version to use for ABFE"""
    tool_version = "0.2.7"
    _tool_key = "deeporigin.abfe-end-to-end"  # Tool key for ABFE jobs

    _max_atom_count: int = 100_000

    def __init__(self, parent):
        super().__init__(parent)

        self._params.end_to_end = utils._load_params("abfe_end_to_end")

    def get_results(self) -> pd.DataFrame | None:
        """get ABFE results and return in a dataframe.

        This method returns a dataframe showing the results of ABFE runs associated with this simulation session. The ligand file name and ΔG are shown, together with user-supplied properties"""

        files_client = getattr(self.parent._platform_clients, "FilesApi", None)

        files = utils.find_files_on_ufa(
            tool="ABFE",
            protein=self.parent.protein.file_path.name,
            client=files_client,
        )

        results_files = [file for file in files if file.endswith("/results.csv")]

        if len(results_files) == 0:
            print("No ABFE results found for this protein.")
            return None

        files_api.download_files(
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

        df2 = ligands_to_dataframe(self.parent.ligands)
        df2["SMILES"] = df2["Ligand"]
        df2.drop(columns=["Ligand", "initial_smiles"], inplace=True)

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

        # re‑index your DataFrame
        df = df[new_order]

        if get_notebook_environment() == "marimo":
            import marimo as mo

            return mo.plain(df)

        else:
            return df

    def set_test_run(self, value: int = 1):
        """set test_run paramemter in abfe parameters"""

        utils._set_test_run(self._params.end_to_end, value)

    @beartype
    def run_end_to_end(
        self,
        *,
        ligands: Optional[list[Ligand]] = None,
        re_run: bool = False,
        _output_dir_path: Optional[str] = None,
    ) -> list[Job] | None:
        """Method to run an end-to-end ABFE run.

        Args:
            ligands: List of ligand to run. Defaults to None. When None, all ligands in the object will be run. To view a list of valid ligands, use the `.show_ligands()` method"""

        if ligands is None:
            ligands = self.parent.ligands

        # check that there is a prepared system for each ligand
        for ligand in ligands:
            if ligand.name not in self.parent._prepared_systems:
                raise ValueError(
                    f"Complex with Ligand {ligand.name} is not prepared. Please prepare the system using the `prepare` method of Complex."
                )

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

        ligand_names = [os.path.basename(ligand._remote_path) for ligand in ligands]

        df = get_dataframe(
            tool_key=tool_mapper["ABFE"],
            only_with_status=["Succeeded", "Running", "Queued", "Created"],
            include_metadata=True,
            resolve_user_names=False,
            _platform_clients=self.parent._platform_clients,
        )

        # filter to find relevant jobs
        df = df[
            df["metadata"].apply(
                lambda d: isinstance(d, dict) and d.get("ligand_file") in ligand_names
            )
        ]

        # get a list of ligands that have already been run
        if len(df) > 0:
            ligands_already_run = list(
                df["metadata"].apply(
                    lambda d: isinstance(d, dict) and d.get("ligand_file")
                )
            )
        else:
            ligands_already_run = []

        if re_run:
            # need to re-run, so don't remove already run ligands
            ligands_to_run = ligand_names
        else:
            # no re-run, remove already run ligands
            ligands_to_run = [
                ligand for ligand in ligand_names if ligand not in ligands_already_run
            ]

        if len(ligands_to_run) == 0 and not re_run:
            print(
                "All requested ligands have already been run, or are queued to run. To re-run, set re_run=True"
            )
            return

        if self.jobs is None:
            self.jobs = []

        jobs_for_this_run = []

        for ligand_path in ligands_to_run:
            metadata = dict(
                protein_file=os.path.basename(self.parent.protein._remote_path),
                ligand_file=os.path.basename(ligand_path),
            )

            job_id = utils._start_tool_run(
                metadata=metadata,
                ligand1_path=ligand_path,
                protein_path=self.parent.protein._remote_path,
                params=self._params.end_to_end,
                tool="ABFE",
                tool_version=self.tool_version,
                _platform_clients=self.parent._platform_clients,
                _output_dir_path=_output_dir_path,
            )

            job = Job.from_id(job_id, _platform_clients=self.parent._platform_clients)

            job._viz_func = self._render_progress
            job._name_func = self._name_job

            job.sync()

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
                )

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

        files_api.download_files(files_to_download, client=files_client)

        from deeporigin_molstar.src.viewers import ProteinViewer

        protein_viewer = ProteinViewer(
            data=str(LOCAL_BASE / remote_pdb_file), format="pdb"
        )
        html_content = protein_viewer.render_trajectory(
            str(LOCAL_BASE / remote_xtc_file)
        )

        from deeporigin_molstar import JupyterViewer

        JupyterViewer.visualize(html_content)

    @classmethod
    @beartype
    def parse_progress(cls, job: Job) -> dict:
        """parse progress from a ABFE job"""

        steps = [
            "init",
            "complex",
            "ligand",
            "simple_md",
            "solvation",
            "binding",
            "delta_g",
        ]

        if len(job._progress_reports) == 0:
            return {step: "NotStarted" for step in steps}

        try:
            data = job._progress_reports[0]

            if data is None:
                progress = {step: "NotStarted" for step in steps}
                progress["init"] = "Running"
                return progress
            else:
                data = json.loads(data)

            if "cmd" in data and data["cmd"] == "FEP Results":
                return {step: "Succeeded" for step in steps}

            if "status" in data and data["status"] == "Initiating":
                progress = {step: "NotStarted" for step in steps}
                progress["init"] = "Running"
                return progress

            status_value = job._status[0]

            # If the overall status is Succeeded, return a dictionary with every key set to "Succeeded".
            if status_value == "Succeeded":
                return {step: "Succeeded" for step in steps}

            current_step = data["run_name"]

            # Validate the input step
            if current_step not in steps:
                raise ValueError(
                    f"Invalid process step provided: {current_step}. "
                    f"Valid steps are: {', '.join(steps)}."
                )

            progress = {}
            for step in steps:
                if step == current_step:
                    progress[step] = "Running"
                    # Once we hit the current step, stop processing further steps.
                    break
                else:
                    progress[step] = "Succeeded"

            # if the job failed, mark the step that is running as failed
            if job._status[0] == "Failed":
                progress[current_step] = "Failed"

        except Exception:
            progress = {step: "Indeterminate" for step in steps}
            progress["init"] = "Indeterminate"

        return progress

    @classmethod
    @beartype
    def _name_job(cls, job: Job) -> str:
        """utility function to name a job using inputs to that job"""
        try:
            return f"ABFE run using <code>{job._metadata[0]['protein_file']}</code> and <code>{job._metadata[0]['ligand_file']}</code>"
        except Exception:
            return "ABFE run"

    @classmethod
    @beartype
    def _render_progress(cls, job: Job) -> str:
        """
        Render HTML for a Mermaid diagram where each node is drawn as arounded rectangle
        with a color indicating its status.

        Any node not specified in the node_status dict willdefault to "notStarted".
        """

        from deeporigin.utils.notebook import mermaid_to_html

        statuses = cls.parse_progress(job)

        # Define the fixed nodes in the diagram.
        nodes = [
            "init(Init)",
            "complex(Complex Prep)",
            "ligand(Ligand Prep)",
            "solvation(Solvation FEP)",
            "simple_md(Simple MD)",
            "binding(Binding FEP)",
            "delta_g(ΔG)",
        ]

        # Build node definitions. For each node, use the providedstatus or default to "notStarted".
        node_defs = ""
        for node in nodes:
            label = node.split("(")[0]
            status = statuses.get(label, "NotStarted")
            node_defs += f"    {node}:::{status};\n"

        # Define the fixed edges of the diagram.
        edges = """
        init --> complex;
        init --> ligand;
        ligand ----> solvation;
        solvation --> delta_g;
        complex ---> simple_md --> binding -->delta_g;
        """

        # Build the complete Mermaid diagram definition.
        mermaid_code = f"""
    graph LR;
        %% Define styles for statuses:
        classDef NotStarted   fill:#cccccc,stroke:#333,stroke-width:2px;
        classDef Queued    fill:#cccccc,stroke:#222,stroke-width:2px;
        classDef Succeeded   fill:#90ee90,stroke:#333,stroke-width:2px;
        classDef Running      fill:#87CEFA,stroke:#333,stroke-width:2px;
        classDef Failed    fill:#ff7f7f,stroke:#333,stroke-width:2px;

    {node_defs}
    {edges}
        """

        # Render the diagram using your helper function.
        mermaid_html = mermaid_to_html(mermaid_code)

        # Define HTML for the legend. Each status is displayed asa colored span.
        legend_html = """
        <div style="margin-top: 20px; font-family: sans-serif;">
          <span style="background-color:#cccccc; color: black;padding:2px 4px; margin: 0 8px;">NotStarted</span>
          <span style="background-color:#90ee90; color: black;padding:2px 4px; margin: 0 8px;">Suceedeed</span>
          <span style="background-color:#87CEFA; color: black;padding:2px 4px; margin: 0 8px;">Running</span>
          <span style="background-color:#ff7f7f; color: black;padding:2px 4px; margin: 0 8px;">Failed</span>
        </div>
        """
        # Display the legend below the Mermaid diagram.
        return mermaid_html + legend_html
