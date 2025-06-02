"""This module encapsulates methods to run RBFE and show RBFE results on Deep Origin.

The RBFE object instantiated here is contained in the Complex class is meant to be used within that class."""

import os
from pathlib import Path
from typing import Optional

from beartype import beartype

from deeporigin.drug_discovery import chemistry as chem
from deeporigin.drug_discovery import utils
from deeporigin.drug_discovery.structures.ligand import Ligand
from deeporigin.drug_discovery.workflow_step import WorkflowStep
from deeporigin.tools.job import Job
from deeporigin.utils.core import PrettyDict

LOCAL_BASE = Path.home() / ".deeporigin"


class RBFE(WorkflowStep):
    """class to handle RBFE-related tasks within the Complex class.

    Objects instantiated here are meant to be used within the Complex class."""

    """tool version to use for RBFE"""
    tool_version = "0.2.7"
    _tool_key = "deeporigin.rbfe-end-to-end"  # Tool key for RBFE jobs

    _max_atom_count: int = 100_000

    def __init__(self, parent):
        super().__init__(parent)
        self._params = PrettyDict()
        self._params.end_to_end = utils._load_params("rbfe_end_to_end")

    def _render_progress(self, job) -> str:
        """Render progress visualization for a job."""
        # TODO: Implement RBFE-specific progress visualization
        return "RBFE Progress Visualization"

    def _name_job(self, job) -> str:
        """Generate a name for a job."""
        return f"RBFE run using <code>{job._metadata[0]['proteprotein_filein_id']}</code>, <code>{job._metadata[0]['ligand1_file']}</code>, and <code>{job._metadata[0]['ligand2_file']}</code>"

    def get_results(self):
        """Fetch RBFE results and return in a dataframe.

        This method returns a dataframe showing the results of RBFE runs associated with this simulation session. The ligand file name, SMILES string and ΔΔG are shown."""

        raise NotImplementedError("RBFE is not implemented yet.")

    def show_results(self):
        """Show RBFE results in a dataframe.

        This method returns a dataframe showing the results of RBFE runs associated with this simulation session. The ligand file name, 2-D structure, and ΔΔG are shown."""

        df = self.get_results()

        if len(df) == 0:
            return

        # convert SMILES to aligned images
        smiles1_list = list(df["SMILES1"])
        smiles2_list = list(df["SMILES2"])
        df.drop("SMILES1", axis=1, inplace=True)
        df.drop("SMILES2", axis=1, inplace=True)

        df["Structure1"] = chem.smiles_list_to_base64_png_list(smiles1_list)
        df["Structure2"] = chem.smiles_list_to_base64_png_list(smiles2_list)

        # Use escape=False to allow the <img> tags to render as images
        from IPython.display import HTML, display

        display(HTML(df.to_html(escape=False)))

    @beartype
    def run_ligand_pair(
        self,
        *,
        ligand1: Ligand,
        ligand2: Ligand,
        re_run: bool = False,
        _output_dir_path: Optional[str] = None,
    ) -> list[Job] | None:
        """Method to run an end-to-end RBFE run.

        Args:
            ligand1: First ligand to run.
            ligand2: Second ligand to run.
            re_run: Whether to re-run the job if it already exists.
            _output_dir_path: Path to the output directory.
        """

        # check that there is a prepared system for each ligand
        for ligand in [ligand1, ligand2]:
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

        if self.jobs is None:
            self.jobs = []

        metadata = dict(
            protein_file=os.path.basename(self.parent.protein._remote_path),
            ligand1_file=os.path.basename(ligand1._remote_path),
            ligand2_file=os.path.basename(ligand2._remote_path),
        )

        job_id = utils._start_tool_run(
            metadata=metadata,
            ligand1_path=ligand1._remote_path,
            ligand2_path=ligand2._remote_path,
            protein_path=self.parent.protein._remote_path,
            params=self._params.end_to_end,
            tool="RBFE",
            tool_version=self.tool_version,
            _platform_clients=self.parent._platform_clients,
            _output_dir_path=_output_dir_path,
        )

        job = Job.from_id(job_id, _platform_clients=self.parent._platform_clients)

        job._viz_func = self._render_progress
        job._name_func = self._name_job

        job.sync()

        self.jobs.append(job)

        return job
