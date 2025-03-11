"""this module contains tools for docking on Deep Origin"""

import os
from dataclasses import dataclass
from typing import Optional

import more_itertools
import pandas as pd
from beartype import beartype
from deeporigin import chemistry
from deeporigin.data_hub import api
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from deeporigin.utils.core import hash_strings

DB_PROTEINS = "Proteins"
DB_DOCKING = "Docking"

COL_PROTEIN = "Protein"
COL_SMILES_HASH = "SMILES_hash"
COL_JOBID = "JobID"
COL_OUTPUT = "OutputFile"
COL_RESULT = "ResultFile"
COL_ALL_SMILES_HASH = "All_SMILES_hash"


DOCKING_CSV_DIR = os.path.join(os.path.expanduser("~"), ".deeporigin", "docking", "csv")
DOCKING_SDF_DIR = os.path.join(os.path.expanduser("~"), ".deeporigin", "docking", "sdf")
os.makedirs(DOCKING_CSV_DIR, exist_ok=True)
os.makedirs(DOCKING_SDF_DIR, exist_ok=True)


@dataclass
class Docking:
    """class to run docking on Deep Origin"""

    protein: chemistry.Protein
    smiles_strings: list[str]

    _proteins_db: Optional[dict] = None
    _docking_db: Optional[dict] = None

    _all_smiles_hash: Optional[str] = None

    def _repr_pretty_(self, p, cycle):
        """pretty print a Docking object"""

        if cycle:
            p.text("Docking(...)")
        else:
            p.text("Docking(")

            p.text(f"protein={self.protein.name}")
            p.text(f" with {len(self.smiles_strings)} ligands")
            p.text(")")

    @beartype
    def _ensure_dbs(self):
        """ensure that there are databases for Docking on the data hub

        Ensures the following DBs:

        - protein (list of protein files)
        - Docking

        """

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

        # docking
        database = _ensure_database(DB_DOCKING)
        required_columns = [
            dict(name=COL_PROTEIN, type="text"),
            dict(name=COL_JOBID, type="text"),
            dict(name=COL_SMILES_HASH, type="text"),
            dict(name=COL_ALL_SMILES_HASH, type="text"),
            dict(name=COL_OUTPUT, type="file"),
            dict(name=COL_RESULT, type="file"),
        ]
        database = _ensure_columns(
            database=database,
            required_columns=required_columns,
        )
        self._docking_db = database

    @classmethod
    def from_dir(cls, directory: str) -> "Docking":
        """create a docking object from a directory of PDB and SDF files

        Args:
            directory (str): directory containing PDB and SDF files
        """

        sdf_files = sorted(
            [
                os.path.join(directory, f)
                for f in os.listdir(directory)
                if f.lower().endswith(".sdf")
            ]
        )

        smiles_strings = []
        for file in sdf_files:
            smiles_strings.extend(chemistry.sdf_to_smiles(file))

        smiles_strings = sorted(set(smiles_strings))

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

        # Create the docking instance
        docking = cls(
            smiles_strings=smiles_strings,
            protein=protein,
            _all_smiles_hash=hash_strings(smiles_strings),
        )

        return docking

    def get_results(self) -> pd.DataFrame:
        """get results from Deep Origin and return a pandas dataframe with results"""

        df = pd.DataFrame(
            api.get_dataframe(
                DB_DOCKING,
                return_type="dict",
                use_file_names=False,
            )
        )

        df = df[self._all_smiles_hash == df["All_SMILES_hash"]]

        # download the CSV files for this run
        file_ids = list(df[COL_OUTPUT].dropna())
        existing_files = os.listdir(DOCKING_CSV_DIR)
        existing_files = ["_file:" + file for file in existing_files]
        missing_files = list(set(file_ids) - set(existing_files))
        if len(missing_files) > 0:
            api.download_files(
                file_ids=missing_files,
                use_file_names=False,
                save_to_dir=DOCKING_CSV_DIR,
            )

        all_dfs = []
        for file_id in file_ids:
            file_loc = os.path.join(DOCKING_CSV_DIR, file_id.replace("_file:", ""))
            df = pd.read_csv(file_loc)
            all_dfs.append(df)

        df = pd.concat(all_dfs)

        return df

    def show_results(self):
        """show results of bulk Docking run in a table, rendering 2D structures of molecules"""

        df = self.get_results()

        from IPython.display import HTML, display

        smiles_list = list(df["SMILES"])
        images = chemistry.smiles_list_to_base64_png_list(smiles_list)
        df["Structure"] = images
        df.drop("SMILES", axis=1, inplace=True)
        display(HTML(df.to_html(escape=False)))

    def connect(self):
        """Connects the local instantiation of the simulation to Deep Origin.

        If contained ligand or protein files do not exist on Deep Origin, they will be uploaded. The connect method also connects to existing runs, if any.

        """

        self._ensure_dbs()

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

    def dock(
        self,
        *,
        box_size: tuple[float, float, float],
        pocket_center: tuple[float, float, float],
        batch_size: int = 30,
    ):
        """Run bulk docking on Deep Origin. Ligands will be split into batches based on the batch_size argument, and will run in parallel on Deep Origin clusters.

        Args:
            box_size (tuple[float, float, float]): box size
            pocket_center (tuple[float, float, float]): pocket center
            batch_size (int, optional): batch size. Defaults to 30.

        """

        chunks = list(more_itertools.chunked(self.smiles_strings, batch_size))

        params = dict(
            box_size=box_size,
            pocket_center=pocket_center,
            smiles_list=self.smiles_strings,
        )

        database_columns = self._proteins_db.cols + self._docking_db.cols

        for chunk in chunks:
            params["smiles_list"] = chunk
            _start_bulk_docking_run_and_log(
                params=params,
                protein_id=self.protein._do_id,
                database_columns=database_columns,
                all_smiles_hash=self._all_smiles_hash,
            )


@beartype
def _start_bulk_docking_run_and_log(
    *,
    protein_id: str,
    database_columns: list,
    params: dict,
    all_smiles_hash: str,
):
    """starts a single run of docking end to end and logs it in the Docking database. Internal function. Do not use.

    Args:
        protein_id (str): protein ID
        params (dict): parameters for the ABFE end-to-end job
        database_columns (list): list of database columns dicts
        all_smiles_hash (str): hash of all smiles

    """

    # first check if we actually need to run this
    df = pd.DataFrame(api.get_dataframe(DB_DOCKING, return_type="dict"))

    existing_hashes = list(df[COL_SMILES_HASH])

    smiles_hash = hash_strings(params["smiles_list"])

    if smiles_hash in existing_hashes:
        print("This tranche of ligands has already been docked. Skipping...")
        return

    from deeporigin.tools import run

    tool_key = "deeporigin.bulk-docking"

    # make a new row
    response = api.make_database_rows(DB_DOCKING, n_rows=1)
    row_id = response.rows[0].hid

    # write hash of all ligands
    api.set_cell_data(
        all_smiles_hash,
        column_id=COL_ALL_SMILES_HASH,
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
            "columnId": COL_OUTPUT,
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
