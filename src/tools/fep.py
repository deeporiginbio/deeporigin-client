"""Module to help work on FEP calculations. This module provides the FEP class, that allows you to run FEP calculations on Deep Origin."""

import os
from dataclasses import dataclass
from typing import Optional

import pandas as pd
from beartype import beartype
from deeporigin import chemistry
from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database

# constants
DB_ABFE = "ABFE2"
DB_RBFE = "RBFE"
DB_PROTEINS = "Proteins"
DB_LIGANDS = "Ligands"

COL_PROTEIN = "Protein"
COL_LIGAND = "Ligand"  # for ABFE
COL_LIGAND1 = "Ligand1"  # for RBFE
COL_LIGAND2 = "Ligand2"  # for RBFE
COL_STEP = "Step"
COL_JOBID = "JobID"
COL_OUTPUT = "OutputFile"
COL_RESULT = "ResultFile"
COL_DELTA_G = "FEP ΔG (kcal/mol)"
COL_DELTA_DELTA_G = "FEP ΔΔG (kcal/mol)"

FEP_DIR = os.path.join(os.path.expanduser("~"), ".deeporigin", "fep")
os.makedirs(FEP_DIR, exist_ok=True)


@dataclass
class FEP:
    """Class for FEP simulations. This class can be used to run FEP calculations on Deep Origin. This class can contain N ligands and a single protein."""

    ligands: list[chemistry.Ligand]
    protein: chemistry.Protein

    _ligands_db: Optional[dict] = None
    _proteins_db: Optional[dict] = None
    _abfe_db: Optional[dict] = None
    _rbfe_db: Optional[dict] = None

    params_abfe_end_to_end = _load_params("abfe_end_to_end")
    params_rbfe_end_to_end = _load_params("rbfe_end_to_end")

    @classmethod
    def from_dir(cls, directory: str) -> "FEP":
        """initialize an FEP class given some files in a directory

        Args:
            directory (str): directory containing ligand and protein files.

        Protein file should be in PDB format. Ligand files should be in SDF format. Each SDF file should contain a single molecule. If your SDF files contain more than one molecule, use `deeporigin.chemistry.split_sdf_file` to split them into separate files.


        """

        sdf_files = sorted(
            [
                os.path.join(directory, f)
                for f in os.listdir(directory)
                if f.lower().endswith(".sdf")
            ]
        )
        ligands = [chemistry.Ligand(sdf_file) for sdf_file in sdf_files]

        pdb_files = [
            os.path.join(directory, f)
            for f in os.listdir(directory)
            if f.lower().endswith(".pdb")
        ]

        if len(pdb_files) != 1:
            raise ValueError(
                f"Expected exactly one PDB file in the directory, but found {len(pdb_files)}."
            )
        protein_file = pdb_files[0]
        protein = chemistry.Protein(protein_file)

        # Create the ABFE instance
        fep = cls(
            ligands=ligands,
            protein=protein,
        )

        return fep

    def get_abfe_results(self):
        """Get ABFE results.

        Returns a dataframe containing results of ABFE runs. The retuned dataframe has the following columns:

        - Ligand: contains the file name of the original SDF file
        - Step: Step in the ABFE flow. `End-to-end` indicates the end-to-end flow, starting from ligand and protein, and going all the way to the outputs.
        - Status: `Succeeded` or `Running`, etc. Status of job on Deep Origin.
        """

        df = pd.DataFrame(
            api.get_dataframe(
                DB_ABFE,
                return_type="dict",
                use_file_names=False,
            )
        )

        from deeporigin.tools.utils import query_run_status

        # we only care about proteins and ligands corresponding to this session
        df = df[df[COL_PROTEIN] == self.protein._do_id]
        valid_ligands = [ligand._do_id for ligand in self.ligands]
        df = df[df[COL_LIGAND].isin(valid_ligands)]

        df["Status"] = ""
        df.loc[~df["OutputFile"].isna(), "Status"] = "Succeeded"

        df.loc[df["OutputFile"].isna(), "Status"] = df.loc[
            df["OutputFile"].isna(), COL_JOBID
        ].apply(query_run_status)

        # for each row, fill in delta_gs if needed

        # download all files for delta_gs
        file_ids = list(df[COL_RESULT].dropna())
        existing_files = os.listdir(FEP_DIR)
        existing_files = ["_file:" + file for file in existing_files]
        missing_files = list(set(file_ids) - set(existing_files))
        if len(missing_files) > 0:
            api.download_files(
                file_ids=missing_files,
                use_file_names=False,
                save_to_dir=FEP_DIR,
            )

        # open each file, read the delta_g, write it to
        # the local dataframe
        for idx, row in df.iterrows():
            if not pd.isna(row[COL_RESULT]) and pd.isna(row[COL_DELTA_G]):
                file_id = row[COL_RESULT].replace("_file:", "")
                delta_g = float(
                    pd.read_csv(os.path.join(FEP_DIR, file_id))["Total"].iloc[0]
                )
                df.loc[idx, COL_DELTA_G] = delta_g

        # drop some columns
        df.drop("Validation Status", axis=1, inplace=True)
        df.drop("JobID", axis=1, inplace=True)
        df.drop(COL_OUTPUT, axis=1, inplace=True)
        df.drop(COL_RESULT, axis=1, inplace=True)
        df.drop("ID", axis=1, inplace=True)
        df.drop(COL_PROTEIN, axis=1, inplace=True)

        # map ligand IDs to ligand file names
        # Create a mapping dictionary: _do_id -> file
        mapping = {
            ligand._do_id: os.path.basename(ligand.file) for ligand in self.ligands
        }
        smiles_mapping = {
            ligand._do_id: ligand.smiles_string for ligand in self.ligands
        }

        # Replace the values in the 'Ligand' column with the corresponding file
        df["SMILES"] = df[COL_LIGAND].map(smiles_mapping)
        df[COL_LIGAND] = df[COL_LIGAND].map(mapping)

        return df

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
