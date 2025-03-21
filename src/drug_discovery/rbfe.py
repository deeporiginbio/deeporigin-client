"""This module encapsulates methods to run RBFE and show RBFE results on Deep Origin.

The RBFE object instantiated here is contained in the Complex class is meant to be used within that class."""

import pandas as pd
from beartype import beartype
from deeporigin.data_hub import api
from deeporigin.drug_discovery import chemistry as chem
from deeporigin.drug_discovery import utils
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils.core import PrettyDict


class RBFE:
    """class to handle RBFE-related tasks within the Complex class.

    Objects instantiated here are meant to be used within the Complex class."""

    def __init__(self, parent):
        self.parent = parent
        self._params = PrettyDict()

        self._params.end_to_end = utils._load_params("rbfe_end_to_end")

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

        if ligand1_id is None or ligand2_id is None:
            raise DeepOriginException(
                "Both ligand1_id and ligand2_id must be specified."
            )

        if ligand1_id == ligand2_id:
            raise DeepOriginException("ligand1_id and ligand2_id cannot be the same.")

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

        if ligand1_id not in valid_ligand_ids or ligand2_id not in valid_ligand_ids:
            raise DeepOriginException(
                f"Some ligand IDs re not valid. Valid ligand IDs are: {valid_ligand_ids}"
            )

        database_columns = (
            self.parent._db.ligands.cols
            + self.parent._db.proteins.cols
            + self.parent._db.rbfe.cols
        )

        # only run on ligands that have not been run yet
        # first check that there are no existing runs
        df = api.get_dataframe(utils.DB_RBFE)
        df = df[df[utils.COL_PROTEIN] == self.parent.protein._do_id]
        df = df[
            (df[utils.COL_LIGAND1] == ligand1_id)
            & (df[utils.COL_LIGAND2] == ligand2_id)
        ]

        job_id = utils._start_tool_run(
            protein_id=self.parent.protein._do_id,
            ligand1_id=ligand1_id,
            ligand2_id=ligand2_id,
            database_columns=database_columns,
            params=self._params.end_to_end,
            tool=utils.DB_RBFE,
            complex_hash=self.parent._hash,
        )

        self.parent._job_ids[utils.DB_RBFE].append(job_id)
