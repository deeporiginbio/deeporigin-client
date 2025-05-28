"""This module encapsulates methods to run RBFE and show RBFE results on Deep Origin.

The RBFE object instantiated here is contained in the Complex class is meant to be used within that class."""

from beartype import beartype
import pandas as pd

from deeporigin.drug_discovery import chemistry as chem
from deeporigin.drug_discovery import utils
from deeporigin.drug_discovery.workflow_step import WorkflowStep
from deeporigin.utils.core import PrettyDict


class RBFE(WorkflowStep):
    """class to handle RBFE-related tasks within the Complex class.

    Objects instantiated here are meant to be used within the Complex class."""

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
        return f"RBFE run using <code>{job._metadata[0]['protein_id']}</code>, <code>{job._metadata[0]['ligand1_id']}</code>, and <code>{job._metadata[0]['ligand2_id']}</code>"

    def get_results(self):
        """Fetch RBFE results and return in a dataframe.

        This method returns a dataframe showing the results of RBFE runs associated with this simulation session. The ligand file name, SMILES string and ΔΔG are shown."""

        df = self.parent.get_csv_results_for(utils.DB_RBFE)

        if len(df) == 0:
            print("No RBFE results to display.")
            return pd.DataFrame()

        return df

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
    def run_end_to_end(
        self,
        *,
        ligand1_id: str,
        ligand2_id: str,
    ):
        """Run end-to-end ABFE run on a pair of ligands.

        Args:
            ligand1_id (str): ID of ligand 1
            ligand2_id (str): ID of ligand 2


        """

        raise NotImplementedError("RBFE is not implemented yet.")
