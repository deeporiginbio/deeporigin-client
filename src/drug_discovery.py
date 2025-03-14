"""Module to support Drug Discovery workflows using Deep Origin"""

import concurrent.futures
import importlib.resources
import json
import os
from dataclasses import dataclass, field
from importlib.resources import path
from typing import Literal, Optional, get_args

import more_itertools
import pandas as pd
from beartype import beartype
from deeporigin import chemistry as chem
from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from deeporigin.tools.utils import query_run_statuses
from deeporigin.utils.core import PrettyDict, hash_file, hash_strings

Number = float | int

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
COL_LIGAND1 = "Ligand1"  # for ABFE/RBFE
COL_LIGAND2 = "Ligand2"  # for RBFE
COL_CSV_OUTPUT = "OutputFile"  # this will be a CSV
COL_PROTEIN = "Protein"
COL_LIGAND = "Ligand"
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


# example data
with path("deeporigin.data.brd", "brd.pdb") as file_path:
    EXAMPLE_DATA_DIR = file_path.parent


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
                dict(name=COL_LIGAND1, type="text"),
                dict(name=COL_COMPLEX_HASH, type="text"),
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
                dict(name=COL_COMPLEX_HASH, type="text"),
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
        instance = cls(
            ligands=ligands,
            protein=protein,
        )

        return instance

    def connect(self) -> None:
        """Connect instance of Complex to the databases on Deep Origin. This method uploads ligand and protein files if needed, and retrieves job IDs of tasks for all tools, if they exist.

        Before running any tool, it is required to call this method."""

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
        # TODO -- parallelize this
        for tool in list(get_args(VALID_TOOLS)):
            df = pd.DataFrame(
                api.get_dataframe(
                    tool,
                    return_type="dict",
                )
            )
            df = df[df["ComplexHash"] == self._hash]
            self._job_ids[tool] = df[COL_JOBID].tolist()

    @beartype
    def get_status_for(self, tool: VALID_TOOLS) -> dict:
        """Return status for jobs corresponding to a particular tool

        Args:
            tool: one of "Docking", "ABFE", "RBFE"
        """
        return query_run_statuses(self._job_ids[tool])

    def get_status(self):
        """Returns status for all runs for all tools"""

        data = dict()
        for tool in list(get_args(VALID_TOOLS)):
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

    def show_ligands(self):
        """Show all ligands in complex object in a table, rendering ligands as 2D structures"""

        chem.show_ligands(self.ligands)

    def get_csv_results_for(self, tool: VALID_TOOLS):
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

        df = df[self._hash == df[COL_COMPLEX_HASH]]

        # this makes sure that we only retain rows
        # where there is a valid output being generated
        # rows where there is no output (for failed or
        # running jobs) are simply dropped
        df = df.dropna(subset=[COL_CSV_OUTPUT])
        file_ids = list(df[COL_CSV_OUTPUT])

        existing_files = os.listdir(DATA_DIRS[tool])
        existing_files = ["_file:" + file for file in existing_files]
        missing_files = list(set(file_ids) - set(existing_files))
        if len(missing_files) > 0:
            api.download_files(
                file_ids=missing_files,
                use_file_names=False,
                save_to_dir=DATA_DIRS[tool],
            )

        if COL_LIGAND1 in df.columns:
            ligand1_ids = df[COL_LIGAND1].tolist()
        else:
            ligand1_ids = [None for _ in file_ids]

        if COL_LIGAND2 in df.columns:
            ligand2_ids = df[COL_LIGAND2].tolist()
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

                df["SMILES1"] = df[COL_LIGAND1].map(smiles_mapping)
                df["SMILES2"] = df[COL_LIGAND2].map(smiles_mapping)
            elif tool == "ABFE":
                df["Ligand"] = ligand1_id

                df["SMILES"] = df["Ligand"].map(smiles_mapping)

            all_dfs.append(df)

        if len(all_dfs) == 0:
            # no data, so nothing to do
            return pd.DataFrame()

        df = pd.concat(all_dfs)

        return df

    def get_docking_results(self) -> pd.DataFrame:
        """Get docking results from Deep Origin"""

        # to do -- some way to make sure that we handle failed runs, complete runs, etc.
        # status = self.get_status("Docking")

        df1 = self.get_csv_results_for("Docking")

        df2 = chem.ligands_to_dataframe(self.ligands)
        df2["SMILES"] = df2["Ligand"]
        df2.drop(columns=["Ligand"], inplace=True)

        df = pd.merge(df1, df2, on="SMILES", how="inner")
        return df

    def show_docking_results(self):
        """show results of bulk Docking run in a table, rendering 2D structures of molecules"""

        df = self.get_docking_results()

        from IPython.display import HTML, display

        smiles_list = list(df["SMILES"])
        images = chem.smiles_list_to_base64_png_list(smiles_list)
        df["Structure"] = images
        df.drop("SMILES", axis=1, inplace=True)
        display(HTML(df.to_html(escape=False)))

    @beartype
    def dock(
        self,
        *,
        box_size: tuple[Number, Number, Number],
        pocket_center: tuple[Number, Number, Number],
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

        df = pd.DataFrame(
            api.get_dataframe(
                DB_DOCKING,
                return_type="dict",
                use_file_names=False,
            )
        )

        df = df[self._hash == df[COL_COMPLEX_HASH]]

        # remove failed runs
        statuses = query_run_statuses(df["JobID"].tolist())
        df["Status"] = df["JobID"].replace(statuses)
        df = df["Failed" != df["Status"]]

        smiles_strings = [ligand.smiles_string for ligand in self.ligands]

        chunks = list(more_itertools.chunked(smiles_strings, batch_size))

        params = dict(
            box_size=box_size,
            pocket_center=pocket_center,
        )

        database_columns = self._db.proteins.cols + self._db.docking.cols

        for chunk in chunks:
            params["smiles_list"] = chunk
            smiles_hash = hash_strings(params["smiles_list"])

            if smiles_hash in df[COL_SMILES_HASH].tolist():
                print("Skipping this tranche because this has already been run...")
                continue

            job_id = _start_tool_run(
                params=params,
                protein_id=self.protein._do_id,
                database_columns=database_columns,
                complex_hash=self._hash,
                tool="Docking",
            )

            self._job_ids["Docking"].append(job_id)

    @beartype
    def run_abfe_end_to_end(
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
        df = df[(df[COL_LIGAND1].isin(ligand_ids))]

        already_run_ligands = set(df[COL_LIGAND1])
        ligand_ids = set(ligand_ids) - already_run_ligands

        for ligand_id in ligand_ids:
            job_id = _start_tool_run(
                protein_id=self.protein._do_id,
                ligand1_id=ligand_id,
                database_columns=database_columns,
                params=self._params.abfe_end_to_end,
                tool=DB_ABFE,
                complex_hash=self._hash,
            )

            self._job_ids[DB_ABFE].append(job_id)

    @beartype
    def run_rbfe_end_to_end(
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

        if ligand1_id not in valid_ligand_ids or ligand2_id not in valid_ligand_ids:
            raise DeepOriginException(
                f"Some ligand IDs re not valid. Valid ligand IDs are: {valid_ligand_ids}"
            )

        database_columns = (
            self._db.ligands.cols + self._db.proteins.cols + self._db.rbfe.cols
        )

        # only run on ligands that have not been run yet
        # first check that there are no existing runs
        df = api.get_dataframe(DB_RBFE)
        df = df[df[COL_PROTEIN] == self.protein._do_id]
        df = df[(df[COL_LIGAND1] == ligand1_id) & (df[COL_LIGAND2] == ligand2_id)]

        job_id = _start_tool_run(
            protein_id=self.protein._do_id,
            ligand1_id=ligand1_id,
            ligand2_id=ligand2_id,
            database_columns=database_columns,
            params=self._params.rbfe_end_to_end,
            tool=DB_RBFE,
            complex_hash=self._hash,
        )

        self._job_ids[DB_RBFE].append(job_id)

    def get_abfe_results(self):
        """get ABFE results and return in a dataframe.

        This method returns a dataframe showing the results of ABFE runs associated with this simulation session. The ligand file name and ΔG are shown, together with user-supplied properties"""

        df1 = self.get_csv_results_for(DB_ABFE)

        if len(df1) == 0:
            print("No ABFE results to display.")
            return

        df1["ID"] = df1["Ligand"]
        df1.drop(columns=["Ligand", "SMILES"], inplace=True)

        df2 = chem.ligands_to_dataframe(self.ligands)
        df2["SMILES"] = df2["Ligand"]
        df2.drop(columns=["Ligand"], inplace=True)

        df = pd.merge(df1, df2, on="ID", how="inner")

        return df

    def show_abfe_results(self):
        """Show ABFE results in a dataframe.

        This method returns a dataframe showing the results of ABFE runs associated with this simulation session. The ligand file name, 2-D structure, and ΔG are shown."""

        df = self.get_abfe_results()

        # convert SMILES to aligned images
        smiles_list = list(df["SMILES"])
        df.drop("SMILES", axis=1, inplace=True)

        df["Structure"] = chem.smiles_list_to_base64_png_list(smiles_list)

        # Use escape=False to allow the <img> tags to render as images
        from IPython.display import HTML, display

        display(HTML(df.to_html(escape=False)))

    def get_rbfe_results(self):
        """Fetch RBFE results and return in a dataframe.

        This method returns a dataframe showing the results of RBFE runs associated with this simulation session. The ligand file name, SMILES string and ΔΔG are shown."""

        df = self.get_csv_results_for(DB_RBFE)

        if len(df) == 0:
            print("No RBFE results to display.")
            return pd.DataFrame()

        return df

    def show_rbfe_results(self):
        """Show RBFE results in a dataframe.

        This method returns a dataframe showing the results of RBFE runs associated with this simulation session. The ligand file name, 2-D structure, and ΔΔG are shown."""

        df = self.get_rbfe_results()

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
def _start_tool_run(
    *,
    params: dict,
    database_columns: list,
    tool: VALID_TOOLS,
    protein_id: str,
    complex_hash: str,
    ligand1_id: Optional[str] = None,
    ligand2_id: Optional[str] = None,
) -> str:
    """starts a single run of ABFE end to end and logs it in the ABFE database. Internal function. Do not use.

    Args:
        protein_id (str): protein ID
        ligand_id (str): ligand ID
        params (dict): parameters for the ABFE end-to-end job
        database_columns (list): list of database columns dicts

    """

    # input validation
    if tool == "ABFE" and ligand1_id is None:
        raise ValueError("ligand1_id is required for ABFE")

    if tool == "RBFE" and (ligand1_id is None or ligand2_id is None):
        raise ValueError("ligand1_id and ligand2_id is required for RBFE")

    # tool key mapper
    tool_key_mapper = dict(
        ABFE="deeporigin.abfe-end-to-end",
        RBFE="deeporigin.rbfe-end-to-end",
        Docking="deeporigin.bulk-docking",
    )

    from deeporigin.tools import run

    # make a new row
    response = api.make_database_rows(tool, n_rows=1)
    row_id = response.rows[0].hid

    # a protein is needed for ABFE, RBFE, and docking
    params["protein"] = {
        "columnId": COL_PROTEIN,
        "rowId": protein_id,
        "databaseId": DB_PROTEINS,
    }

    # input ligand files
    if tool == "RBFE":
        params["ligand1"] = {
            "columnId": COL_LIGAND,
            "rowId": ligand1_id,
            "databaseId": DB_LIGANDS,
        }

        params["ligand2"] = {
            "columnId": COL_LIGAND,
            "rowId": ligand2_id,
            "databaseId": DB_LIGANDS,
        }
    elif tool == "ABFE":
        params["ligand"] = {
            "columnId": COL_LIGAND,
            "rowId": ligand1_id,
            "databaseId": DB_LIGANDS,
        }

    # output files
    if tool == "RBFE":
        outputs = {
            "output_file": {
                "columnId": COL_RESULT,
                "rowId": row_id,
                "databaseId": DB_RBFE,
            },
            "rbfe_results_summary": {
                "columnId": COL_CSV_OUTPUT,
                "rowId": row_id,
                "databaseId": DB_RBFE,
            },
        }
    elif tool == "ABFE":
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
    elif tool == "Docking":
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
        tool_key=tool_key_mapper[tool],
        cols=database_columns,
    )

    # write job ID
    api.set_cell_data(
        job_id,
        column_id=COL_JOBID,
        row_id=row_id,
        database_id=tool,
    )

    # write ligand1 ID
    if ligand1_id is not None:
        api.set_cell_data(
            ligand1_id,
            column_id=COL_LIGAND1,
            row_id=row_id,
            database_id=tool,
        )

    # write ligand2 ID
    if ligand2_id is not None:
        api.set_cell_data(
            ligand2_id,
            column_id=COL_LIGAND2,
            row_id=row_id,
            database_id=tool,
        )

    # write protein ID
    api.set_cell_data(
        protein_id,
        column_id=COL_PROTEIN,
        row_id=row_id,
        database_id=tool,
    )

    # write complex hash
    api.set_cell_data(
        complex_hash,
        column_id=COL_COMPLEX_HASH,
        row_id=row_id,
        database_id=tool,
    )

    # write SMILEShash
    if tool == "Docking":
        smiles_hash = smiles_hash = hash_strings(params["smiles_list"])
        api.set_cell_data(
            smiles_hash,
            column_id=COL_SMILES_HASH,
            row_id=row_id,
            database_id=tool,
        )

    # write step
    if tool in ["ABFE", "RBFE"]:
        api.set_cell_data(
            "End-to-end",
            column_id=COL_STEP,
            row_id=row_id,
            database_id=tool,
        )

    return job_id
