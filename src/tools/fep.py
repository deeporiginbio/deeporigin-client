"""Module to help work on FEP calculations. This module provides the FEP class, that allows you to run FEP calculations on Deep Origin."""

import os
from dataclasses import dataclass

import pandas as pd
from beartype import beartype
from deeporigin import chemistry
from deeporigin.data_hub import api


@dataclass
class FEP:
    """Class for FEP simulations. This class can be used to run FEP calculations on Deep Origin. This class can contain N ligands and a single protein."""

    ligands: list[chemistry.Ligand]
    protein: chemistry.Protein

    def show_abfe_results(self):
        """Show ABFE results in a dataframe.

        This method returns a dataframe showing the results of ABFE runs associated with this simulation session. The ligand file name, 2-D structure, and ΔG are shown."""

        df = self.get_abfe_results()

        if len(df) == 0:
            print("No ABFE results to display. Start a run first.")
            return

        # convert SMILES to aligned images
        smiles_list = list(df["SMILES"])
        df.drop("SMILES", axis=1, inplace=True)

        df["Structure"] = chemistry.smiles_list_to_base64_png_list(smiles_list)

        # Use escape=False to allow the <img> tags to render as images
        from IPython.display import HTML, display

        display(HTML(df.to_html(escape=False)))


@beartype
def _start_rbfe_run_and_log(
    *,
    protein_id: str,
    ligand1_id: str,
    ligand2_id: str,
    params: dict,
    database_columns: list,
):
    """starts a single run of ABFE end to end and logs it in the ABFE database. Internal function. Do not use.

    Args:
        protein_id (str): protein ID
        ligand_id (str): ligand ID
        params (dict): parameters for the ABFE end-to-end job
        database_columns (list): list of database columns dicts

    """

    from deeporigin.tools import run

    tool_key = "deeporigin.rbfe-end-to-end"

    # make a new row
    response = api.make_database_rows(DB_RBFE, n_rows=1)
    row_id = response.rows[0].hid

    # write ligand1 ID
    api.set_cell_data(
        ligand1_id,
        column_id=COL_LIGAND1,
        row_id=row_id,
        database_id=DB_RBFE,
    )

    # write ligand2 ID
    api.set_cell_data(
        ligand2_id,
        column_id=COL_LIGAND2,
        row_id=row_id,
        database_id=DB_RBFE,
    )

    # write protein ID
    api.set_cell_data(
        protein_id,
        column_id=COL_PROTEIN,
        row_id=row_id,
        database_id=DB_RBFE,
    )

    # write step
    api.set_cell_data(
        "End-to-end",
        column_id=COL_STEP,
        row_id=row_id,
        database_id=DB_RBFE,
    )

    # start job
    params["protein"] = {
        "columnId": COL_PROTEIN,
        "rowId": protein_id,
        "databaseId": DB_PROTEINS,
    }

    params["ligand1"] = {
        "columnId": COL_LIGAND,
        "rowId": ligand1_id,
        "databaseId": DB_LIGANDS,
    }

    params["ligand2"] = {
        "columnId": COL_LIGAND,
        "rowId": ligand1_id,
        "databaseId": DB_LIGANDS,
    }

    outputs = {
        "output_file": {
            "columnId": COL_OUTPUT,
            "rowId": row_id,
            "databaseId": DB_RBFE,
        },
        "rbfe_results_summary": {
            "columnId": COL_RESULT,
            "rowId": row_id,
            "databaseId": DB_RBFE,
        },
    }

    job_id = run._process_job(
        inputs=params,
        outputs=outputs,
        tool_key=tool_key,
        cols=database_columns,
    )

    # write job ID
    api.set_cell_data(
        job_id,
        column_id=COL_JOBID,
        row_id=row_id,
        database_id=DB_RBFE,
    )
