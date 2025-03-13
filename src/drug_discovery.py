"""Module to support Drug Discovery workflows using Deep Origin"""

import concurrent.futures
import importlib.resources
import json
import os
from dataclasses import dataclass, field
from typing import Literal, Optional

import more_itertools
import pandas as pd
from beartype import beartype
from deeporigin import chemistry as chem
from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from deeporigin.tools.utils import query_run_statuses
from deeporigin.utils.core import PrettyDict, hash_file, hash_strings

# constants
DB_ABFE = "ABFE"
DB_RBFE = "RBFE"
DB_DOCKING = "Docking"
DB_PROTEINS = "Proteins"
DB_LIGANDS = "Ligands"

VALID_TOOLS = Literal[DB_ABFE, DB_RBFE, DB_DOCKING]

COL_DELTA_DELTA_G = "FEP ΔΔG (kcal/mol)"
COL_DELTA_G = "FEP ΔG (kcal/mol)"
COL_JOBID = "JobID"
COL_LIGAND = "Ligand"  # for ABFE
COL_LIGAND1 = "Ligand1"  # for RBFE
COL_LIGAND2 = "Ligand2"  # for RBFE
COL_CSV_OUTPUT = "OutputFile"  # this will be a CSV
COL_PROTEIN = "Protein"
COL_RESULT = "ResultFile"
COL_SMILES_HASH = "SMILESHash"
COL_COMPLEX_HASH = "ComplexHash"
COL_STEP = "Step"


DATA_DIRS = dict()

DATA_DIRS[DB_ABFE] = os.path.join(os.path.expanduser("~"), ".deeporigin", DB_ABFE)
DATA_DIRS[DB_RBFE] = os.path.join(os.path.expanduser("~"), ".deeporigin", DB_RBFE)
DATA_DIRS[DB_DOCKING] = os.path.join(os.path.expanduser("~"), ".deeporigin", DB_DOCKING)


os.makedirs(DATA_DIRS[DB_ABFE], exist_ok=True)
os.makedirs(DATA_DIRS[DB_RBFE], exist_ok=True)
os.makedirs(DATA_DIRS[DB_DOCKING], exist_ok=True)


@beartype
def _load_params(tool: str) -> PrettyDict:
    """load params for various tools, reading from JSON files"""

    with importlib.resources.open_text("deeporigin.json", f"{tool}.json") as f:
        return PrettyDict(json.load(f))


@beartype
def _ensure_dbs() -> dict:
    """Ensure that there are databases for FEP on the data hub.

    Ensures the following DBs:
      - ligand (list of ligand files)
      - protein (list of protein files)
      - ABFE (ligand, protein, step, job_id, output, results, delta_g)
      - RBFE (ligand1, ligand2, protein, step, job_id, output, results, delta_g)
      - docking (protein, job_id, smiles_hash, complex_hash, output, results)
    """
    databases = PrettyDict()

    # Define a helper function to process a single database
    def process_db(db_key: str, db_name, required_columns: list) -> tuple[str, object]:
        """ensure that a db exists, and ensure that columns exist"""

        db = _ensure_database(db_name)
        db = _ensure_columns(database=db, required_columns=required_columns)

        return db_key, db

    # Mapping each database key to its parameters (name and required columns)
    tasks = {
        "ligands": (
            DB_LIGANDS,
            [dict(name=COL_LIGAND, type="file")],
        ),
        "proteins": (
            DB_PROTEINS,
            [dict(name=COL_PROTEIN, type="file")],
        ),
        "abfe": (
            DB_ABFE,
            [
                dict(name=COL_PROTEIN, type="text"),
                dict(name=COL_LIGAND, type="text"),
                dict(name=COL_STEP, type="text"),
                dict(name=COL_JOBID, type="text"),
                dict(name=COL_CSV_OUTPUT, type="file"),
                dict(name=COL_RESULT, type="file"),
                dict(name=COL_DELTA_G, type="float"),
            ],
        ),
        "rbfe": (
            DB_RBFE,
            [
                dict(name=COL_PROTEIN, type="text"),
                dict(name=COL_LIGAND1, type="text"),
                dict(name=COL_LIGAND2, type="text"),
                dict(name=COL_STEP, type="text"),
                dict(name=COL_JOBID, type="text"),
                dict(name=COL_CSV_OUTPUT, type="file"),
                dict(name=COL_RESULT, type="file"),
                dict(name=COL_DELTA_DELTA_G, type="float"),
            ],
        ),
        "docking": (
            DB_DOCKING,
            [
                dict(name=COL_PROTEIN, type="text"),
                dict(name=COL_JOBID, type="text"),
                dict(name=COL_SMILES_HASH, type="text"),
                dict(name=COL_COMPLEX_HASH, type="text"),
                dict(name=COL_CSV_OUTPUT, type="file"),
                dict(name=COL_RESULT, type="file"),
            ],
        ),
    }

    # Use a ThreadPoolExecutor for I/O-bound network calls
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_to_key = {
            executor.submit(process_db, key, db_name, req_cols): key
            for key, (db_name, req_cols) in tasks.items()
        }

        for future in concurrent.futures.as_completed(future_to_key):
            key, db = future.result()
            databases[key] = db

    return databases


@dataclass
class Complex:
    """class to represent a set of a protein and 1 or many ligands"""

    ligands: list[chem.Ligand]
    protein: chem.Protein

    # these params are not user facing
    _db: Optional[dict] = None
    _params: PrettyDict = PrettyDict()

    """stores a hash of all ligands and the protein. This will be computed post initialization"""
    _hash: Optional[str] = None

    """stores job ids for all jobs, organized by tool. """
    _job_ids: dict = field(
        default_factory=lambda: {DB_DOCKING: [], DB_ABFE: [], DB_RBFE: []}
    )

    def __post_init__(self):
        """various post init tasks"""

        # hash protein and ligands
        protein_hash = hash_file(self.protein.file)
        ligands_hash = hash_strings([ligand.smiles_string for ligand in self.ligands])
        self._hash = hash_strings([protein_hash, ligands_hash])

        # load params for all tools
        self._params.abfe_end_to_end = _load_params("abfe_end_to_end")
        self._params.rbfe_end_to_end = _load_params("rbfe_end_to_end")

    @classmethod
    def from_dir(cls, directory: str) -> "Complex":
        """initialize an FEP class given some files in a directory

        Args:
            directory (str): directory containing ligand and protein files.

        Protein file should be in PDB format. Ligand files should be in SDF format. Each SDF file should contain a single molecule. If your SDF files contain more than one molecule, use `deeporigin.chemistry.split_sdf_file` to split them into separate files.


        """

        # figure out all the SDF files
        sdf_files = sorted(
            [
                os.path.join(directory, f)
                for f in os.listdir(directory)
                if f.lower().endswith(".sdf")
            ]
        )
        mols = []
        for sdf_file in sdf_files:
            mols_from_this_sdf_file = chem.read_molecules_in_sdf_file(sdf_file)
            if len(mols_from_this_sdf_file) == 1:
                # there's only one molecule in thie SDF file, so we should track the file it came from.
                mols_from_this_sdf_file[0]["file"] = sdf_file

            mols.extend(mols_from_this_sdf_file)

        ligands = []

        for mol in mols:
            ligands.append(chem.Ligand(**mol))

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
        protein = chem.Protein(protein_file)

        # Create the Complex instance
        complex = cls(
            ligands=ligands,
            protein=protein,
        )

        return complex

    def connect(self) -> None:
        """connect to the databases on Deep Origin"""

        if self._db is None:
            self._db = _ensure_dbs()

        # ensure that ligands are uploaded
        df = api.get_dataframe(DB_LIGANDS)

        for ligand in self.ligands:
            if ligand.file is None:
                # no ligand file, we only have a SMILES string. this means that there is no need to upload or otherwise manage this ligand

                continue

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

        # fetch all relevant jobIDs
        df = pd.DataFrame(api.get_dataframe(DB_DOCKING, return_type="dict"))
        df = df[df["ComplexHash"] == self._hash]
        self._job_ids["Docking"] = df["JobID"].tolist()

    @beartype
    def get_status(self, tool: VALID_TOOLS) -> dict:
        return query_run_statuses(self._job_ids[tool])

    def _repr_pretty_(self, p, cycle):
        """pretty print a Docking object"""

        if cycle:
            p.text("Complex(...)")
        else:
            p.text("Complex(")

            p.text(f"protein={self.protein.name}")
            p.text(f" with {len(self.ligands)} ligands")
            p.text(")")

    def show_ligands(self):
        """show all ligands in complex object"""

        chem.show_ligands(self.ligands)

    def get_csv_results_for(self, tool: VALID_TOOLS):
        """get CSV results for a particular tool and combine them as need be"""

        df = pd.DataFrame(
            api.get_dataframe(
                tool,
                return_type="dict",
                use_file_names=False,
            )
        )

        df = df[self._hash == df[COL_COMPLEX_HASH]]

        file_ids = list(df[COL_CSV_OUTPUT].dropna())
        existing_files = os.listdir(DATA_DIRS[tool])
        existing_files = ["_file:" + file for file in existing_files]
        missing_files = list(set(file_ids) - set(existing_files))
        if len(missing_files) > 0:
            api.download_files(
                file_ids=missing_files,
                use_file_names=False,
                save_to_dir=DATA_DIRS[tool],
            )

        all_dfs = []
        for file_id in file_ids:
            file_loc = os.path.join(DATA_DIRS[tool], file_id.replace("_file:", ""))
            df = pd.read_csv(file_loc)
            all_dfs.append(df)

        df = pd.concat(all_dfs)
        return df

    def get_docking_results(self) -> pd.DataFrame:
        """get docking results from Deep Origin"""

        # to do -- some way to make sure that we handle failed runs, complete runs, etc.
        status = self.get_status("Docking")

        # download the CSV files for this run

        return self.get_csv_results_for("Docking")

    def show_docking_results(self):
        """show results of bulk Docking run in a table, rendering 2D structures of molecules"""

        df = self.get_docking_results()

        from IPython.display import HTML, display

        smiles_list = list(df["SMILES"])
        images = chem.smiles_list_to_base64_png_list(smiles_list)
        df["Structure"] = images
        df.drop("SMILES", axis=1, inplace=True)
        display(HTML(df.to_html(escape=False)))

    def dock(
        self,
        *,
        box_size: tuple[float, float, float],
        pocket_center: tuple[float, float, float],
        batch_size: int = 32,
    ):
        """Run bulk docking on Deep Origin. Ligands will be split into batches based on the batch_size argument, and will run in parallel on Deep Origin clusters.

        Args:
            box_size (tuple[float, float, float]): box size
            pocket_center (tuple[float, float, float]): pocket center
            batch_size (int, optional): batch size. Defaults to 30.

        """

        if self.protein._do_id is None:
            raise DeepOriginException(
                "Protein must be uploaded to Deep Origin before docking."
            )

        smiles_strings = [ligand.smiles_string for ligand in self.ligands]

        chunks = list(more_itertools.chunked(smiles_strings, batch_size))

        params = dict(
            box_size=box_size,
            pocket_center=pocket_center,
        )

        database_columns = self._db.proteins.cols + self._db.docking.cols

        for chunk in chunks:
            params["smiles_list"] = chunk
            job_id = _start_bulk_docking_run_and_log(
                params=params,
                protein_id=self.protein._do_id,
                database_columns=database_columns,
                complex_hash=self._hash,
            )

            self._job_ids["Docking"].append(job_id)

    @beartype
    def abfe_end_to_end(
        self,
        *,
        ligand_ids: Optional[list[str]] = None,
    ):
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
            self._db.ligands.cols + self._db.proteins.cols + self._db.abfe.cols
        )

        # only run on ligands that have not been run yet
        # first check that there are no existing runs
        df = api.get_dataframe(DB_ABFE)
        df = df[df[COL_PROTEIN] == self.protein._do_id]
        df = df[(df[COL_LIGAND].isin(ligand_ids))]

        already_run_ligands = set(df[COL_LIGAND])
        ligand_ids = set(ligand_ids) - already_run_ligands

        for ligand_id in ligand_ids:
            job_id = _start_abfe_run_and_log(
                protein_id=self.protein._do_id,
                ligand_id=ligand_id,
                database_columns=database_columns,
                params=self._params.abfe_end_to_end,
            )

            self._job_ids["ABFE"].append(job_id)

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
        file_ids = list(df[COL_CSV_OUTPUT].dropna())
        existing_files = os.listdir(DATA_DIRS[DB_ABFE])
        existing_files = ["_file:" + file for file in existing_files]
        missing_files = list(set(file_ids) - set(existing_files))
        if len(missing_files) > 0:
            api.download_files(
                file_ids=missing_files,
                use_file_names=False,
                save_to_dir=DATA_DIRS[DB_ABFE],
            )

        # open each file, read the delta_g, write it to
        # the local dataframe
        for idx, row in df.iterrows():
            if not pd.isna(row[COL_CSV_OUTPUT]) and pd.isna(row[COL_DELTA_G]):
                file_id = row[COL_CSV_OUTPUT].replace("_file:", "")
                delta_g = float(
                    pd.read_csv(os.path.join(DATA_DIRS[DB_ABFE], file_id))[
                        "Total"
                    ].iloc[0]
                )
                df.loc[idx, COL_DELTA_G] = delta_g

        # drop some columns
        df.drop("Validation Status", axis=1, inplace=True)
        df.drop("JobID", axis=1, inplace=True)
        df.drop(COL_CSV_OUTPUT, axis=1, inplace=True)
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


# TO DO: the internal network requests here need to be parallelized
@beartype
def _start_bulk_docking_run_and_log(
    *,
    protein_id: str,
    database_columns: list,
    params: dict,
    complex_hash: str,
) -> str:
    """starts a single run of docking end to end and logs it in the Docking database. Internal function. Do not use.

    Args:
        protein_id (str): protein ID
        params (dict): parameters for the ABFE end-to-end job
        database_columns (list): list of database columns dicts
        complex_hash (str): hash of complex of ligands and protein

    """

    # first check if we actually need to run this
    df = pd.DataFrame(api.get_dataframe(DB_DOCKING, return_type="dict"))

    existing_hashes = list(df[COL_SMILES_HASH])

    smiles_hash = hash_strings(params["smiles_list"])

    if smiles_hash in existing_hashes:
        print("This tranche of ligands has already been docked. Skipping...")

        # TODO transmit job ID correctly, and do so only when  the job has succeeded or is running (Failed jobs should be ignored)
        return

    from deeporigin.tools import run

    tool_key = "deeporigin.bulk-docking"

    # make a new row
    response = api.make_database_rows(DB_DOCKING, n_rows=1)
    row_id = response.rows[0].hid

    # start job
    params["pdb_file"] = {
        "columnId": COL_PROTEIN,
        "rowId": protein_id,
        "databaseId": DB_PROTEINS,
    }

    outputs = {
        "data_file": {
            "columnId": COL_CSV_OUTPUT,
            "rowId": row_id,
            "databaseId": DB_DOCKING,
        },
        "results_sdf": {
            "columnId": COL_RESULT,
            "rowId": row_id,
            "databaseId": DB_DOCKING,
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
        database_id=DB_DOCKING,
    )

    # write hash of all ligands
    api.set_cell_data(
        complex_hash,
        column_id=COL_COMPLEX_HASH,
        row_id=row_id,
        database_id=DB_DOCKING,
    )

    # write hash of ligands we're docking
    api.set_cell_data(
        smiles_hash,
        column_id=COL_SMILES_HASH,
        row_id=row_id,
        database_id=DB_DOCKING,
    )

    # write protein ID
    api.set_cell_data(
        protein_id,
        column_id=COL_PROTEIN,
        row_id=row_id,
        database_id=DB_DOCKING,
    )

    return job_id


@beartype
def _start_abfe_run_and_log(
    *,
    protein_id: str,
    ligand_id: str,
    params: dict,
    database_columns: list,
) -> str:
    """starts a single run of ABFE end to end and logs it in the ABFE database. Internal function. Do not use.

    Args:
        protein_id (str): protein ID
        ligand_id (str): ligand ID
        params (dict): parameters for the ABFE end-to-end job
        database_columns (list): list of database columns dicts

    """

    from deeporigin.tools import run

    tool_key = "deeporigin.abfe-end-to-end"

    # TODO -- take this out of this loop and use n_rows>1
    # so that we can make all rows with a single call
    # make a new row
    response = api.make_database_rows(DB_ABFE, n_rows=1)
    row_id = response.rows[0].hid

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
            "columnId": COL_RESULT,
            "rowId": row_id,
            "databaseId": DB_ABFE,
        },
        "abfe_results_summary": {
            "columnId": COL_CSV_OUTPUT,
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

    return job_id
