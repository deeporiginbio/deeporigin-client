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


@dataclass
class Docking:
    protein: chemistry.Protein
    smiles_strings: list[str]

    _proteins_db: Optional[dict] = None
    _docking_db: Optional[dict] = None

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
        )

        return docking

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

    def dock(self, *, batch_size: int = 30):
        """dock on Deep Origin"""

        chunks = list(more_itertools.chunked(self.smiles_strings, batch_size))

        params = dict(
            box_size=[
                15,
                15,
                15,
            ],
            pocket_center=[
                13,
                -6,
                22,
            ],
            smiles_list=self.smiles_strings,
        )

        database_columns = self._proteins_db.cols + self._docking_db.cols

        for chunk in chunks:
            params["smiles_list"] = chunk
            _start_bulk_docking_run_and_log(
                params=params,
                protein_id=self.protein._do_id,
                database_columns=database_columns,
            )


@beartype
def _start_bulk_docking_run_and_log(
    *,
    protein_id: str,
    database_columns: list,
    params: dict,
):
    """starts a single run of ABFE end to end and logs it in the ABFE database. Internal function. Do not use.

    Args:
        protein_id (str): protein ID
        ligand_id (str): ligand ID
        params (dict): parameters for the ABFE end-to-end job
        database_columns (list): list of database columns dicts

    """

    # first check if we actually need to run this
    df = pd.DataFrame(api.get_dataframe(DB_DOCKING, return_type="dict"))

    existing_hashes = list(df[COL_SMILES_HASH])

    smiles_hash = hash_strings(params["smiles_list"])

    if smiles_hash in existing_hashes:
        print("This tranche of ligands has already been docked. Skipping...")
        return

    from deeporigin.tools import run

    tool_key = "sgs-test.bulk-docking"

    # make a new row
    response = api.make_database_rows(DB_DOCKING, n_rows=1)
    row_id = response.rows[0].hid

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
