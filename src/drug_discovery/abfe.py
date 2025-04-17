"""This module encapsulates methods to run ABFE and show ABFE results on Deep Origin.

The ABFE object instantiated here is contained in the Complex class is meant to be used within that class."""

import json
import os
import pathlib
from typing import Literal, Optional
import zipfile

from beartype import beartype
import pandas as pd

from deeporigin.data_hub import api
from deeporigin.drug_discovery import chemistry as chem
from deeporigin.drug_discovery import utils
from deeporigin.drug_discovery.structures.ligand import ligands_to_dataframe
from deeporigin.drug_discovery.workflow_step import WorkflowStep
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools.job import Job


class ABFE(WorkflowStep):
    """class to handle ABFE-related tasks within the Complex class.

    Objects instantiated here are meant to be used within the Complex class."""

    def __init__(self, parent):
        super().__init__(parent)

        self._params.end_to_end = utils._load_params("abfe_end_to_end")

    def get_results(self) -> pd.DataFrame:
        """get ABFE results and return in a dataframe.

        This method returns a dataframe showing the results of ABFE runs associated with this simulation session. The ligand file name and ΔG are shown, together with user-supplied properties"""

        df1 = self.parent.get_csv_results_for(utils.DB_ABFE)

        if len(df1) == 0:
            print("No ABFE results to display.")
            return

        df1["ID"] = df1["Ligand1"]
        df1.drop(columns=["Ligand1", "SMILES"], inplace=True)

        df2 = ligands_to_dataframe(self.parent.ligands)
        df2["SMILES"] = df2["Ligand"]
        df2.drop(columns=["Ligand"], inplace=True)

        df = pd.merge(
            df1,
            df2,
            on="ID",
            how="inner",
            validate="one_to_one",
        )

        return df

    def show_results(self):
        """Show ABFE results in a dataframe.

        This method returns a dataframe showing the results of ABFE runs associated with this simulation session. The ligand file name, 2-D structure, and ΔG are shown."""

        df = self.get_results()

        if len(df) == 0:
            return

        # convert SMILES to aligned images
        smiles_list = list(df["SMILES"])
        df.drop("SMILES", axis=1, inplace=True)

        df["Structure"] = chem.smiles_list_to_base64_png_list(smiles_list)

        # show structure first
        new_order = ["Structure"] + [col for col in df.columns if col != "Structure"]

        # re‑index your DataFrame
        df = df[new_order]

        # Use escape=False to allow the <img> tags to render as images
        from IPython.display import HTML, display

        display(HTML(df.to_html(escape=False)))

    @beartype
    def run_end_to_end(
        self,
        *,
        ligand_ids: Optional[list[str]] = None,
    ):
        """Method to run an end-to-end ABFE run.

        Args:
            ligand_ids (Optional[str], optional): List of ligand IDs to run. Defaults to None. When None, all ligands in the object will be run. To view a list of valid ligand IDs, use the `.show_ligands()` method"""

        if ligand_ids is None:
            ligand_ids = [ligand._do_id for ligand in self.parent.ligands]

        # check that protein ID is valid
        if self.parent.protein._do_id is None:
            raise DeepOriginException(
                "Protein has not been uploaded yet. Use .connect() first."
            )

        # check that ligand IDs are valid
        valid_ligand_ids = [ligand._do_id for ligand in self.parent.ligands]

        if None in valid_ligand_ids:
            raise DeepOriginException(
                "Some ligands have not been uploaded yet. Use .connect() first."
            )

        if not set(ligand_ids).issubset(valid_ligand_ids):
            raise DeepOriginException(
                f"Some ligand IDs re not valid. Valid ligand IDs are: {valid_ligand_ids}"
            )

        database_columns = (
            self.parent._db.ligands.cols
            + self.parent._db.proteins.cols
            + self.parent._db.abfe.cols
        )

        # only run on ligands that have not been run yet
        # first check that there are no existing runs
        df = api.get_dataframe(utils.DB_ABFE)
        df = df[df[utils.COL_PROTEIN] == self.parent.protein._do_id]
        df = df[(df[utils.COL_LIGAND1].isin(ligand_ids))]

        already_run_ligands = set(df[utils.COL_LIGAND1])
        ligand_ids = set(ligand_ids) - already_run_ligands

        if len(ligand_ids) == 0:
            print("All requested ligands have already been run.")
            return

        for ligand_id in ligand_ids:
            job_id = utils._start_tool_run(
                protein_id=self.parent.protein._do_id,
                ligand1_id=ligand_id,
                database_columns=database_columns,
                params=self._params.end_to_end,
                tool=utils.DB_ABFE,
                complex_hash=self.parent._hash,
            )

            job = Job.from_ids([job_id])

            job._viz_func = self._render_progress
            job._name_func = self._name_job

            self.jobs.append(job)

    @beartype
    def show_trajectory(
        self,
        ligand_id: str,
        step: Literal["md", "binding"],
    ):
        """Show the system trajectory FEP run.

        Args:
            ligand_id (str): The ID of the ligand to show the trajectory for.
            step (Literal["md", "abfe"]): The step to show the trajectory for.
        """

        valid_ids = [ligand._do_id for ligand in self.parent.ligands]

        if None in valid_ids:
            self.parent.connect()

            valid_ids = [ligand._do_id for ligand in self.parent.ligands]

        if ligand_id not in valid_ids:
            raise DeepOriginException(
                f"Ligand ID {ligand_id} not found in the list of ligands. Should be one of {valid_ids}"
            )

        # get the files for the run
        files = self.parent.get_result_files_for(tool="ABFE", ligand_ids=[ligand_id])

        file = files[0]

        # Get the file path and create directory path with same name
        file_path = pathlib.Path(file)
        dir_name = f"{file_path.stem}-execution"
        dir_path = file_path.parent / dir_name

        # Check if directory already exists
        if not os.path.exists(dir_path):
            os.makedirs(dir_path, exist_ok=True)

            # Unzip the file into the directory
            print(f"Extracting trajectory files to {dir_path}...")
            with zipfile.ZipFile(file_path, "r") as zip_ref:
                zip_ref.extractall(dir_path)

        pdb_file = dir_path / "execution/protein/ligand/systems/complex/complex.pdb"

        if step == "binding":
            xtc_file = (
                dir_path
                / "execution/protein/ligand/binding/binding/window_1/Prod_1/_allatom_trajectory_40ps.xtc"
            )
        else:
            xtc_file = (
                dir_path
                / "execution/protein/ligand/simple_md/md/prod/_allatom_trajectory_40ps.xtc"
            )

        from deeporigin_molstar.src.viewers import ProteinViewer

        protein_viewer = ProteinViewer(data=str(pdb_file), format="pdb")
        html_content = protein_viewer.render_trajectory(str(xtc_file))

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

        data = job._progress_reports[0]

        if data is None:
            progress = {step: "NotStarted" for step in steps}
            progress["init"] = "Running"
            return progress
        else:
            data = json.loads(data)

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

        return progress

    @classmethod
    @beartype
    def _name_job(cls, job: Job) -> str:
        """utility function to name a job using inputs to that job"""
        try:
            return f"ABFE run using <code>{job._metadata[0]['protein_id']}</code> and <code>{job._metadata[0]['ligand1_id']}</code>"
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
            "init",
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
