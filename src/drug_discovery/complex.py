"""Module to support Drug Discovery workflows using Deep Origin"""

import concurrent.futures
import json
import os
from dataclasses import dataclass, field
from typing import Optional, get_args

import pandas as pd
from beartype import beartype
from deeporigin.data_hub import api
from deeporigin.drug_discovery import chemistry as chem
from deeporigin.drug_discovery import utils
from deeporigin.drug_discovery.abfe import ABFE
from deeporigin.drug_discovery.docking import Docking
from deeporigin.drug_discovery.rbfe import RBFE
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from deeporigin.tools.utils import query_run_statuses
from deeporigin.utils.core import PrettyDict, hash_file, hash_strings

DATA_DIRS = dict()

DATA_DIRS[utils.DB_ABFE] = os.path.join(
    os.path.expanduser("~"), ".deeporigin", utils.DB_ABFE
)
DATA_DIRS[utils.DB_RBFE] = os.path.join(
    os.path.expanduser("~"), ".deeporigin", utils.DB_RBFE
)
DATA_DIRS[utils.DB_DOCKING] = os.path.join(
    os.path.expanduser("~"), ".deeporigin", utils.DB_DOCKING
)


os.makedirs(DATA_DIRS[utils.DB_ABFE], exist_ok=True)
os.makedirs(DATA_DIRS[utils.DB_RBFE], exist_ok=True)
os.makedirs(DATA_DIRS[utils.DB_DOCKING], exist_ok=True)


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
            utils.DB_LIGANDS,
            [dict(name=utils.COL_LIGAND, type="file")],
        ),
        "proteins": (
            utils.DB_PROTEINS,
            [dict(name=utils.COL_PROTEIN, type="file")],
        ),
        "abfe": (
            utils.DB_ABFE,
            [
                dict(name=utils.COL_PROTEIN, type="text"),
                dict(name=utils.COL_LIGAND1, type="text"),
                dict(name=utils.COL_COMPLEX_HASH, type="text"),
                dict(name=utils.COL_STEP, type="text"),
                dict(name=utils.COL_JOBID, type="text"),
                dict(name=utils.COL_CSV_OUTPUT, type="file"),
                dict(name=utils.COL_RESULT, type="file"),
                dict(name=utils.COL_DELTA_G, type="float"),
            ],
        ),
        "rbfe": (
            utils.DB_RBFE,
            [
                dict(name=utils.COL_PROTEIN, type="text"),
                dict(name=utils.COL_LIGAND1, type="text"),
                dict(name=utils.COL_LIGAND2, type="text"),
                dict(name=utils.COL_COMPLEX_HASH, type="text"),
                dict(name=utils.COL_STEP, type="text"),
                dict(name=utils.COL_JOBID, type="text"),
                dict(name=utils.COL_CSV_OUTPUT, type="file"),
                dict(name=utils.COL_RESULT, type="file"),
                dict(name=utils.COL_DELTA_DELTA_G, type="float"),
            ],
        ),
        "docking": (
            utils.DB_DOCKING,
            [
                dict(name=utils.COL_PROTEIN, type="text"),
                dict(name=utils.COL_JOBID, type="text"),
                dict(name=utils.COL_SMILES_HASH, type="text"),
                dict(name=utils.COL_COMPLEX_HASH, type="text"),
                dict(name=utils.COL_CSV_OUTPUT, type="file"),
                dict(name=utils.COL_RESULT, type="file"),
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

    # Using a private attribute for ligands with the property decorator below
    _ligands: list[chem.Ligand] = field(repr=False)
    protein: chem.Protein

    # these params are not user facing
    _db: Optional[dict] = None

    """stores a hash of all ligands and the protein. This will be computed post initialization"""
    _hash: Optional[str] = None

    """stores job ids for all jobs, organized by tool. """
    _job_ids: dict = field(
        default_factory=lambda: {
            utils.DB_DOCKING: [],
            utils.DB_ABFE: [],
            utils.DB_RBFE: [],
        }
    )

    def __post_init__(self):
        """various post init tasks"""

        # assign references to the complex in the
        # various child classes
        self.docking = Docking(parent=self)
        self.abfe = ABFE(parent=self)
        self.rbfe = RBFE(parent=self)

        # Initialize the _hash
        self._compute_hash()

    def _compute_hash(self):
        """Compute and update the hash based on protein and ligands"""
        protein_hash = hash_file(self.protein.file)
        ligands_hash = hash_strings([ligand.smiles_string for ligand in self.ligands])
        self._hash = hash_strings([protein_hash, ligands_hash])

    @property
    def ligands(self) -> list[chem.Ligand]:
        """Get the current ligands"""
        return self._ligands

    @ligands.setter
    def ligands(self, new_ligands: list[chem.Ligand]):
        """Set new ligands and recompute the hash"""
        self._ligands = new_ligands
        self._compute_hash()

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
        instance = cls(
            _ligands=ligands,
            protein=protein,
        )

        return instance

    def connect(self) -> None:
        """Connect instance of Complex to the databases on Deep Origin. This method uploads ligand and protein files if needed, and retrieves job IDs of tasks for all tools, if they exist.

        Before running any tool, it is required to call this method."""

        if self._db is None:
            self._db = _ensure_dbs()

        # ensure that ligands are uploaded
        df = api.get_dataframe(utils.DB_LIGANDS)

        for ligand in self.ligands:
            if ligand.file is None:
                # no ligand file, we only have a SMILES string. this means that there is no need to upload or otherwise manage this ligand

                continue

            ligand_file = os.path.basename(ligand.file)
            matching_indices = df.index[df[utils.COL_LIGAND] == ligand_file].tolist()
            if len(matching_indices) == 0:
                print(f"Uploading {ligand.file}...")
                response = api.upload_file_to_new_database_row(
                    database_id=utils.DB_LIGANDS,
                    column_id=utils.COL_LIGAND,
                    file_path=str(ligand.file),
                )

                ligand._do_id = response.rows[0].hid
            else:
                ligand._do_id = matching_indices[0]

        # ensure that protein is uploaded
        df = api.get_dataframe(utils.DB_PROTEINS)
        protein_file = os.path.basename(self.protein.file)
        matching_indices = df.index[df[utils.COL_PROTEIN] == protein_file].tolist()
        if len(matching_indices) == 0:
            print(f"Uploading {self.protein.file}...")
            response = api.upload_file_to_new_database_row(
                database_id=utils.DB_PROTEINS,
                column_id=utils.COL_PROTEIN,
                file_path=str(self.protein.file),
            )

            self.protein._do_id = response.rows[0].hid
        else:
            self.protein._do_id = matching_indices[0]

        # fetch all relevant jobIDs
        # TODO -- parallelize this
        for tool in list(get_args(utils.VALID_TOOLS)):
            df = pd.DataFrame(
                api.get_dataframe(
                    tool,
                    return_type="dict",
                )
            )
            df = df[df["ComplexHash"] == self._hash]
            self._job_ids[tool] = df[utils.COL_JOBID].tolist()

    @beartype
    def get_status_for(self, tool: utils.VALID_TOOLS) -> dict:
        """Return status for jobs corresponding to a particular tool

        Args:
            tool: one of "Docking", "ABFE", "RBFE"
        """
        return query_run_statuses(self._job_ids[tool])

    def get_status(self):
        """Returns status for all runs for all tools"""

        data = dict()
        for tool in list(get_args(utils.VALID_TOOLS)):
            data[tool] = self.get_status_for(tool)
        print(json.dumps(data, indent=2))

    def _repr_pretty_(self, p, cycle):
        """pretty print a Docking object"""

        if cycle:
            p.text("Complex(...)")
        else:
            p.text("Complex(")

            p.text(f"protein={self.protein.name}")
            p.text(f" with {len(self.ligands)} ligands")
            p.text(")")

    @beartype
    def show_ligands(self, *, view: str = "2d", limit: Optional[int] = None):
        """Display ligands in the complex object.

        Args:
            view: Visualization type, either "2d" (default) or "3d".
                 - "2d": Shows ligands in a table with 2D structure renderings
                 - "3d": Shows 3D molecular structures using SDF files
            limit: Optional; maximum number of ligands to display.
                  If None, all ligands will be shown.
        """

        if view == "3d":
            files = [ligand.file for ligand in self.ligands]

            if limit is not None:
                files = files[:limit]

            chem.show_molecules_in_sdf_files(files)
        else:
            if limit is not None:
                chem.show_ligands(self.ligands[:limit])
            else:
                chem.show_ligands(self.ligands)

    @beartype
    def get_result_files_for(
        self,
        *,
        tool: utils.VALID_TOOLS,
        ligand_ids: Optional[list[str]] = None,
    ):
        """Retrieve result files for a specific computational tool used with this complex.

        This method fetches data from the specified tool's database, filters for results
        associated with this complex, downloads any missing result files to the local
        storage directory, and returns paths to all result files.

        Args:
            tool: One of "Docking", "ABFE", "RBFE" - specifies which computational
                 method's results to retrieve
            ligand_ids: Optional; list of ligand IDs to filter results. If None,
                       results for all ligands in the complex will be retrieved.


        """

        df = pd.DataFrame(
            api.get_dataframe(
                tool,
                return_type="dict",
                use_file_names=False,
            )
        )

        df = df[self._hash == df[utils.COL_COMPLEX_HASH]]

        if ligand_ids is not None:
            df = df[df[utils.COL_LIGAND1].isin(ligand_ids)]

        # this makes sure that we only retain rows
        # where there is a valid output being generated
        # rows where there is no output (for failed or
        # running jobs) are simply dropped
        df = df.dropna(subset=[utils.COL_RESULT])

        file_ids = list(df[utils.COL_RESULT])

        existing_files = os.listdir(DATA_DIRS[tool])
        existing_files = ["_file:" + file for file in existing_files]
        missing_files = list(set(file_ids) - set(existing_files))

        if len(missing_files) > 0:
            print("Downloading result files. This can take a while...")
            api.download_files(
                file_ids=missing_files,
                use_file_names=False,
                save_to_dir=DATA_DIRS[tool],
            )
            print("Done.")

        return [
            os.path.join(DATA_DIRS[tool], file_id.replace("_file:", ""))
            for file_id in file_ids
        ]

    def get_csv_results_for(self, tool: utils.VALID_TOOLS):
        """Generic method to get CSV results for a particular tool and combine them as need be

        Args:
            tool: One of "Docking", "ABFE", "RBFE"
        """

        df = pd.DataFrame(
            api.get_dataframe(
                tool,
                return_type="dict",
                use_file_names=False,
            )
        )

        df = df[self._hash == df[utils.COL_COMPLEX_HASH]]

        # this makes sure that we only retain rows
        # where there is a valid output being generated
        # rows where there is no output (for failed or
        # running jobs) are simply dropped
        df = df.dropna(subset=[utils.COL_CSV_OUTPUT])
        file_ids = list(df[utils.COL_CSV_OUTPUT])

        existing_files = os.listdir(DATA_DIRS[tool])
        existing_files = ["_file:" + file for file in existing_files]
        missing_files = list(set(file_ids) - set(existing_files))
        if len(missing_files) > 0:
            api.download_files(
                file_ids=missing_files,
                use_file_names=False,
                save_to_dir=DATA_DIRS[tool],
            )

        if utils.COL_LIGAND1 in df.columns:
            ligand1_ids = df[utils.COL_LIGAND1].tolist()
        else:
            ligand1_ids = [None for _ in file_ids]

        if utils.COL_LIGAND2 in df.columns:
            ligand2_ids = df[utils.COL_LIGAND2].tolist()
        else:
            ligand2_ids = [None for _ in file_ids]

        all_dfs = []
        for file_id, ligand1_id, ligand2_id in zip(file_ids, ligand1_ids, ligand2_ids):
            file_loc = os.path.join(DATA_DIRS[tool], file_id.replace("_file:", ""))
            df = pd.read_csv(file_loc)

            # drop some columns
            drop_columns = ["Ligand1", "Ligand2", "Protein"]
            for col in drop_columns:
                df.drop(
                    col,
                    axis=1,
                    inplace=True,
                    errors="ignore",
                )

            smiles_mapping = {
                ligand._do_id: ligand.smiles_string for ligand in self.ligands
            }

            if tool == "RBFE":
                df["Ligand1"] = ligand1_id
                df["Ligand2"] = ligand2_id

                df["SMILES1"] = df[utils.COL_LIGAND1].map(smiles_mapping)
                df["SMILES2"] = df[utils.COL_LIGAND2].map(smiles_mapping)
            elif tool == "ABFE":
                df["Ligand"] = ligand1_id

                df["SMILES"] = df["Ligand"].map(smiles_mapping)

            all_dfs.append(df)

        if len(all_dfs) == 0:
            # no data, so nothing to do
            return pd.DataFrame()

        df = pd.concat(all_dfs)

        return df
