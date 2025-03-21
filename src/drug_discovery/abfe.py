"""This module encapsulates methods to run ABFE and show ABFE results on Deep Origin.

The ABFE object instantiated here is contained in the Complex class is meant to be used within that class."""

from typing import Optional

import pandas as pd
from beartype import beartype
from deeporigin.data_hub import api
from deeporigin.drug_discovery import chemistry as chem
from deeporigin.drug_discovery import utils
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils.core import PrettyDict


class ABFE:
    """class to handle ABFE-related tasks within the Complex class.

    Objects instantiated here are meant to be used within the Complex class."""

    def __init__(self, parent):
        self.parent = parent
        self._params = PrettyDict()

        self._params.end_to_end = utils._load_params("rbfe_end_to_end")

    def get_results(self) -> pd.DataFrame:
        """get ABFE results and return in a dataframe.

        This method returns a dataframe showing the results of ABFE runs associated with this simulation session. The ligand file name and ΔG are shown, together with user-supplied properties"""

        df1 = self.parent.get_csv_results_for(utils.DB_ABFE)

        if len(df1) == 0:
            print("No ABFE results to display.")
            return

        df1["ID"] = df1["Ligand"]
        df1.drop(columns=["Ligand", "SMILES"], inplace=True)

        df2 = chem.ligands_to_dataframe(self.parent.ligands)
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

        for ligand_id in ligand_ids:
            job_id = utils._start_tool_run(
                protein_id=self.parent.protein._do_id,
                ligand1_id=ligand_id,
                database_columns=database_columns,
                params=self._params.end_to_end,
                tool=utils.DB_ABFE,
                complex_hash=self.parent._hash,
            )

            self.parent._job_ids[utils.DB_ABFE].append(job_id)
