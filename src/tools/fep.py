"""module to work FEP calculations"""

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
from deeporigin.tools import run
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from deeporigin.tools.utils import query_run_status
from IPython.display import HTML, display

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
    """class for FEP simulations"""

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
            directory (str): directory containing ligand and protein files"""

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
        _ensure_columns(
            database=database,
            required_columns=required_columns,
        )
        self._ligands_db = database

        # proteins
        database = _ensure_database(DB_PROTEINS)
        required_columns = [
            dict(name="Protein", type="file"),
        ]
        _ensure_columns(
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
            dict(name="delta_g", type="file"),
        ]
        database = _ensure_columns(
            database=database,
            required_columns=required_columns,
        )
        self._abfe_db = database

    def show_ligands(self):
        """show ligands in a dataframe"""

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
        display(HTML(df.to_html(escape=False)))

    def connect(self):
        """connects the local instantiation of the simulation to Deep Origin"""

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
        """get ABFE results"""

        df = pd.DataFrame(
            api.get_dataframe(
                DB_ABFE,
                return_type="dict",
                use_file_names=False,
            )
        )

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
        file_ids = list(df["ResultFile"].dropna())
        api.download_files(
            file_ids=file_ids,
            use_file_names=False,
            save_to_dir=FEP_DIR,
        )

        # open each file, read the delta_g, write it to
        # the local dataframe
        for idx, row in df.iterrows():
            if not pd.isna(row["ResultFile"]) and pd.isna(row["delta_g"]):
                file_id = row["ResultFile"].replace("_file:", "")
                delta_g = float(
                    pd.read_csv(os.path.join(FEP_DIR, file_id))["Total"].iloc[0]
                )
                df.loc[idx, "delta_g"] = delta_g
                print(delta_g)

        # drop some columns
        df.drop("Validation Status", axis=1, inplace=True)
        df.drop("JobID", axis=1, inplace=True)
        df.drop("OutputFile", axis=1, inplace=True)
        df.drop("ResultFile", axis=1, inplace=True)
        df.drop("ID", axis=1, inplace=True)
        df.drop("Protein", axis=1, inplace=True)

        # map ligand IDs to ligand file names
        # Create a mapping dictionary: _do_id -> file
        mapping = {
            ligand._do_id: os.path.basename(ligand.file) for ligand in self.ligands
        }
        smiles_mapping = {
            ligand._do_id: ligand.smiles_string for ligand in self.ligands
        }

        # Replace the values in the 'Ligand' column with the corresponding file
        df["SMILES"] = df["Ligand"].map(smiles_mapping)
        df["Ligand"] = df["Ligand"].map(mapping)

        return df

    def show_abfe_results(self):
        """show ABFE results in a dataframe"""

        df = self.get_abfe_results()

        # convert SMILES to aligned images
        smiles_list = list(df["SMILES"])
        df.drop("SMILES", axis=1, inplace=True)

        df["Structure"] = chemistry.smiles_list_to_base64_png_list(smiles_list)

        # Use escape=False to allow the <img> tags to render as images
        display(HTML(df.to_html(escape=False)))

    def abfe_end_to_end(
        self,
        *,
        ligand_ids: Optional[str] = None,
    ) -> str:
        """Function to run an end-to-end job"""

        if ligand_ids is None:
            raise NotImplementedError("Not implemented yet")

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

        for ligand_id in ligand_ids:
            start_abfe_run_and_log(
                protein_id=self.protein._do_id,
                ligand_id=ligand_id,
                database_columns=database_columns,
                params=self.params_abfe_end_to_end,
            )


@beartype
def start_abfe_run_and_log(
    *,
    protein_id: str,
    ligand_id: str,
    params: dict,
    database_columns: list,
):
    """starts a single run of ABFE end to end and logs it in the ABFE database

    Args:
        protein_id (str): protein ID
        ligand_id (str): ligand ID
        params (dict): parameters for the ABFE end-to-end job
        database_columns (list): list of database columns dicts

    """

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
        "columnId": "Protein",
        "rowId": protein_id,
        "databaseId": "Proteins",
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
