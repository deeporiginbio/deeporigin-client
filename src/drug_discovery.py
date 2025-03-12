"""Module to support Drug Discovery workflows using Deep Origin"""

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

DATA_DIRS[DB_ABFE] = os.path.join(os.path.expanduser("~"), ".deeporigin", "abfe")
DATA_DIRS[DB_RBFE] = os.path.join(os.path.expanduser("~"), ".deeporigin", "rbfe")
DATA_DIRS[DB_DOCKING] = os.path.join(os.path.expanduser("~"), ".deeporigin", "docking")


os.makedirs(DATA_DIRS[DB_ABFE], exist_ok=True)
os.makedirs(DATA_DIRS[DB_RBFE], exist_ok=True)
os.makedirs(DATA_DIRS[DB_DOCKING], exist_ok=True)


@beartype
def _load_params(step: str) -> PrettyDict:
    """load default values for abfe end to end run"""

    with importlib.resources.open_text("deeporigin.json", f"{step}.json") as f:
        return PrettyDict(json.load(f))


@beartype
def _ensure_dbs() -> dict:
    """ensure that there are databases for FEP on the data hub

    Ensures the following DBs:

    - ligand (list of ligand files)
    - protein (list of protein files)
    - ABFE (ligand, protein, step, job_id, output, results, delta_g)
    - RBFE (ligand1, ligand2, protein, step, job_id, output, results, delta_g)

    """

    databases = PrettyDict()

    # ligands
    database = _ensure_database(DB_LIGANDS)
    required_columns = [
        dict(name=COL_LIGAND, type="file"),
    ]
    database = _ensure_columns(
        database=database,
        required_columns=required_columns,
    )
    databases.ligands = database

    # proteins
    database = _ensure_database(DB_PROTEINS)
    required_columns = [
        dict(name=COL_PROTEIN, type="file"),
    ]
    database = _ensure_columns(
        database=database,
        required_columns=required_columns,
    )
    databases.proteins = database

    # ABFE
    database = _ensure_database(DB_ABFE)
    required_columns = [
        dict(name=COL_PROTEIN, type="text"),
        dict(name=COL_LIGAND, type="text"),
        dict(name=COL_STEP, type="text"),
        dict(name=COL_JOBID, type="text"),
        dict(name=COL_CSV_OUTPUT, type="file"),
        dict(name=COL_RESULT, type="file"),
        dict(name=COL_DELTA_G, type="float"),
    ]
    database = _ensure_columns(
        database=database,
        required_columns=required_columns,
    )
    databases.abfe = database

    # RBFE
    database = _ensure_database(DB_RBFE)
    required_columns = [
        dict(name=COL_PROTEIN, type="text"),
        dict(name=COL_LIGAND1, type="text"),
        dict(name=COL_LIGAND2, type="text"),
        dict(name=COL_STEP, type="text"),
        dict(name=COL_JOBID, type="text"),
        dict(name=COL_CSV_OUTPUT, type="file"),
        dict(name=COL_RESULT, type="file"),
        dict(name=COL_DELTA_DELTA_G, type="float"),
    ]
    database = _ensure_columns(
        database=database,
        required_columns=required_columns,
    )
    databases.rbfe = database

    # docking
    database = _ensure_database(DB_DOCKING)
    required_columns = [
        dict(name=COL_PROTEIN, type="text"),
        dict(name=COL_JOBID, type="text"),
        dict(name=COL_SMILES_HASH, type="text"),
        dict(name=COL_COMPLEX_HASH, type="text"),
        dict(name=COL_CSV_OUTPUT, type="file"),
        dict(name=COL_RESULT, type="file"),
    ]
    database = _ensure_columns(
        database=database,
        required_columns=required_columns,
    )
    databases.docking = database

    return databases


@dataclass
class Complex:
    """class to represent a set of a protein and 1 or many ligands"""

    ligands: list[chem.Ligand]
    protein: chem.Protein

    # these params are not user facing
    _db: Optional[dict] = None
    _params: Optional[dict] = None

    """stores a hash of all ligands and the protein. This will be computed post initialization"""
    _hash: Optional[str] = None

    """stores job ids for all jobs, organized by tool. """
    _job_ids: dict = field(
        default_factory=lambda: {"docking": [], "abfe": [], "rbfe": []}
    )

    def __post_init__(self):
        """compute hash of this object post initialization"""

        protein_hash = hash_file(self.protein.file)
        ligands_hash = hash_strings([ligand.smiles_string for ligand in self.ligands])
        self._hash = hash_strings([protein_hash, ligands_hash])

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
            mols.extend(chem.read_molecules_in_sdf_file(sdf_file))

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
def _start_bulk_docking_run_and_log(
    *,
    protein_id: str,
    database_columns: list,
    params: dict,
    complex_hash: str,
):
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

    return job_id
