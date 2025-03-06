"""Module to help work on FEP calculations. This module provides the FEP class, that allows you to run FEP calculations on Deep Origin."""

import importlib.resources
import json
import os
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from beartype import beartype
from box import Box
from deeporigin import chemistry
from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database

# constants
DB_ABFE = "ABFE2"
DB_PROTEINS = "Proteins"
DB_LIGANDS = "Ligands"
COL_PROTEIN = "Protein"
COL_LIGAND = "Ligand"
COL_STEP = "Step"
COL_JOBID = "JobID"
COL_OUTPUT = "OutputFile"
COL_RESULT = "ResultFile"
COL_DELTA_G = "FEP ΔG (kcal/mol)"

FEP_DIR = os.path.join(os.path.expanduser("~"), ".deeporigin", "fep")
os.makedirs(FEP_DIR, exist_ok=True)


class PrettyDict(Box):
    """A dict subclass with a custom pretty-print representation."""

    def __repr__(self):
        """pretty print a dict"""
        return json.dumps(
            dict(self),
            indent=2,
            ensure_ascii=False,
        )

    def _repr_html_(self):
        """pretty print a dict"""
        self.__repr__()


@beartype
def _load_params(step: str) -> Box:
    """load default values for abfe end to end run"""

    with importlib.resources.open_text("deeporigin.json", f"{step}.json") as f:
        return PrettyDict(json.load(f))


@dataclass
class Ligand:
    """class to represent a ligand (typically backed by a SDF file)"""

    file: Union[str, Path]
    smiles_string: Optional[str] = None
    n_molecules: Optional[int] = None

    # this ID keeps track of whether it is uploaded to deep origin or not
    _do_id: Optional[str] = None

    def __post_init__(self):
        """generates a SMILES if it doesn't exist"""

        # check that there's only one molecule here
        if self.n_molecules is None:
            if chemistry.count_molecules_in_sdf_file(self.file) > 1:
                raise ValueError(
                    "Too many molecules. Expected a single molecule in the SDF file, but got multiple"
                )
            self.n_molecules = 1

        if self.smiles_string is None:
            smiles_string = chemistry.sdf_to_smiles(self.file)
            if len(smiles_string) > 1:
                raise ValueError("Expected a single SMILES strings, but got multiple")
            self.smiles_string = smiles_string[0]

    def _repr_pretty_(self, p, cycle):
        """pretty print a ligand"""

        if cycle:
            p.text("Ligand(...)")
        else:
            p.text("Ligand(")

            with p.group(2, "\n  ", "\n"):
                all_fields = fields(self)
                for idx, field in enumerate(all_fields):
                    value = getattr(self, field.name)
                    p.text(f"{field.name}: {value!r}")
                    # Only add a breakable if this isn't the last field.
                    if idx < len(all_fields) - 1:
                        p.breakable()
            p.text(")")


@dataclass
class Protein:
    """class to represent a protein (typically backed by a PDB file)"""

    file: Union[str, Path]
    name: Optional[str] = None

    # this ID keeps track of whether it is uploaded to deep origin or not
    _do_id: Optional[str] = None

    def __post_init__(self):
        self.file = Path(self.file)
        if self.name is None:
            self.name = self.file.name


@dataclass
class FEP:
    """Class for FEP simulations. This class can be used to run FEP calculations on Deep Origin. This class can contain N ligands and a single protein."""

    ligands: list[Ligand]
    protein: Protein

    _ligands_db: Optional[dict] = None
    _proteins_db: Optional[dict] = None
    _abfe_db: Optional[dict] = None

    params_abfe_end_to_end = _load_params("abfe_end_to_end")

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
        ligands = [Ligand(sdf_file) for sdf_file in sdf_files]

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
        protein = Protein(protein_file)

        # Create the ABFE instance
        fep = cls(
            ligands=ligands,
            protein=protein,
        )

        return fep

    @beartype
    def _ensure_dbs_for_fep(self):
        """ensure that there are databases for FEP on the data hub

        Ensures the following DBs:

        - ligand (list of ligand files)
        - protein (list of protein files)
        - ABFE (ligand, protein, step, job_id, output, results, delta_g)
        - RBFE (ligand1, ligand2, protein, step, job_id, output, results, delta_g)

        """

        # ligands
        database = _ensure_database(DB_LIGANDS)
        required_columns = [
            dict(name=COL_LIGAND, type="file"),
        ]
        database = _ensure_columns(
            database=database,
            required_columns=required_columns,
        )
        self._ligands_db = database

        # proteins
        database = _ensure_database(DB_PROTEINS)
        required_columns = [
            dict(name=COL_PROTEIN, type="file"),
        ]
        database = _ensure_columns(
            database=database,
            required_columns=required_columns,
        )
        self._proteins_db = database

        # ABFE
        database = _ensure_database(DB_ABFE)
        required_columns = [
            dict(name=COL_PROTEIN, type="text"),
            dict(name=COL_LIGAND, type="text"),
            dict(name=COL_STEP, type="text"),
            dict(name=COL_JOBID, type="text"),
            dict(name=COL_OUTPUT, type="file"),
            dict(name=COL_RESULT, type="file"),
            dict(name=COL_DELTA_G, type="file"),
        ]
        database = _ensure_columns(
            database=database,
            required_columns=required_columns,
        )
        self._abfe_db = database

    def show_ligands(self):
        """show ligands in the FEP object in a dataframe. This method visualizes the ligands using core-aligned 2D visualizations. To simply obtain a listing of the ligands in this class, use `.ligands`."""

        # convert SMILES to aligned images
        smiles_list = [ligand.smiles_string for ligand in self.ligands]
        id_list = [ligand._do_id for ligand in self.ligands]
        file_list = [os.path.basename(ligand.file) for ligand in self.ligands]

        images = chemistry.smiles_list_to_base64_png_list(smiles_list)

        data = {
            "ID": id_list,
            "File": file_list,
            "Ligand": images,
        }

        df = pd.DataFrame(data)

        # Use escape=False to allow the <img> tags to render as images
        from IPython.display import HTML, display

        display(HTML(df.to_html(escape=False)))

    def connect(self):
        """Connects the local instantiation of the simulation to Deep Origin.

        If contained ligand or protein files do not exist on Deep Origin, they will be uploaded. The connect method also connects to existing runs, if any.

        """

        self._ensure_dbs_for_fep()

        # ensure that ligands are uploaded
        df = api.get_dataframe(DB_LIGANDS)

        for ligand in self.ligands:
            ligand_file = os.path.basename(ligand.file)
            matching_indices = df.index[df[COL_LIGAND] == ligand_file].tolist()
            if len(matching_indices) == 0:
                print(f"Uploading {ligand.file}...")
                response = api.upload_file_to_new_database_row(
                    database_id=DB_LIGANDS,
                    column_id=COL_LIGAND,
                    file_path=str(ligand.file),
                )

                ligand._do_id = response.rows[0].hid
            else:
                ligand._do_id = matching_indices[0]

        # ensure that protein is uploaded
        df = api.get_dataframe(DB_PROTEINS)
        protein_file = os.path.basename(self.protein.file)
        matching_indices = df.index[df[COL_PROTEIN] == protein_file].tolist()
        if len(matching_indices) == 0:
            print(f"Uploading {self.protein.file}...")
            response = api.upload_file_to_new_database_row(
                database_id=DB_PROTEINS,
                column_id=COL_PROTEIN,
                file_path=str(self.protein.file),
            )

            self.protein._do_id = response.rows[0].hid
        else:
            self.protein._do_id = matching_indices[0]

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

    def abfe_end_to_end(
        self,
        *,
        ligand_ids: Optional[str] = None,
    ) -> str:
        """Method to run an end-to-end ABFE run.

        Args:
            ligand_ids (Optional[str], optional): List of ligand IDs to run. Defaults to None. When None, all ligands in the object will be run. To view a list of valid ligand IDs, use the `.show_ligands()` method"""

        if ligand_ids is None:
            ligand_ids = [ligand._do_id for ligand in self.ligands]

        # check that protein ID is valid
        if self.protein._do_id is None:
            raise DeepOriginException(
                "Protein has not been uploaded yet. Use .connect() first."
            )

        # check that ligand IDs are valid
        valid_ligand_ids = [ligand._do_id for ligand in self.ligands]

        if None in valid_ligand_ids:
            raise DeepOriginException(
                "Some ligands have not been uploaded yet. Use .connect() first."
            )

        if not set(ligand_ids).issubset(valid_ligand_ids):
            raise DeepOriginException(
                f"Some ligand IDs re not valid. Valid ligand IDs are: {valid_ligand_ids}"
            )

        database_columns = (
            self._ligands_db.cols + self._proteins_db.cols + self._abfe_db.cols
        )

        # only run on ligands that have not been run yet
        # first check that there are no existing runs
        df = api.get_dataframe(DB_ABFE)
        df[df["Protein"] == self.protein._do_id]

        df = df[(df[COL_LIGAND].isin(ligand_ids)) & (~pd.isna(df[COL_RESULT]))]

        already_run_ligands = set(df[COL_LIGAND])
        ligand_ids = set(ligand_ids) - already_run_ligands

        for ligand_id in ligand_ids:
            _start_abfe_run_and_log(
                protein_id=self.protein._do_id,
                ligand_id=ligand_id,
                database_columns=database_columns,
                params=self.params_abfe_end_to_end,
            )


@beartype
def _start_abfe_run_and_log(
    *,
    protein_id: str,
    ligand_id: str,
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

    tool_key = "deeporigin.abfe-end-to-end"

    # make a new row
    response = api.make_database_rows(DB_ABFE, n_rows=1)
    row_id = response.rows[0].hid

    # write ligand ID
    api.set_cell_data(
        ligand_id,
        column_id=COL_LIGAND,
        row_id=row_id,
        database_id=DB_ABFE,
    )

    # write protein ID
    api.set_cell_data(
        protein_id,
        column_id=COL_PROTEIN,
        row_id=row_id,
        database_id=DB_ABFE,
    )

    # write step
    api.set_cell_data(
        "End-to-end",
        column_id=COL_STEP,
        row_id=row_id,
        database_id=DB_ABFE,
    )

    # start job
    params["protein"] = {
        "columnId": COL_PROTEIN,
        "rowId": protein_id,
        "databaseId": DB_PROTEINS,
    }

    params["ligand"] = {
        "columnId": COL_LIGAND,
        "rowId": ligand_id,
        "databaseId": DB_LIGANDS,
    }

    outputs = {
        "output_file": {
            "columnId": COL_OUTPUT,
            "rowId": row_id,
            "databaseId": DB_ABFE,
        },
        "abfe_results_summary": {
            "columnId": COL_RESULT,
            "rowId": row_id,
            "databaseId": DB_ABFE,
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
        database_id=DB_ABFE,
    )
