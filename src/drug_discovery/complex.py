"""Module to support Drug Discovery workflows using Deep Origin"""

from dataclasses import dataclass, field
import os
from typing import Optional

from beartype import beartype
import pandas as pd

from deeporigin.data_hub import api
from deeporigin.drug_discovery import chemistry as chem
from deeporigin.drug_discovery import utils
from deeporigin.drug_discovery.abfe import ABFE
from deeporigin.drug_discovery.docking import Docking
from deeporigin.drug_discovery.rbfe import RBFE
from deeporigin.drug_discovery.structures import Ligand, Protein
from deeporigin.files import FilesClient

files_client = FilesClient()


@dataclass
class Complex:
    """class to represent a set of a protein and 1 or many ligands"""

    # Using a private attribute for ligands with the property decorator below
    protein: Protein
    _ligands: list[Ligand] = field(default_factory=list, repr=False)

    # these params are not user facing
    _db: Optional[dict] = None

    """stores a hash of all ligands and the protein. This will be computed post initialization"""
    _hash: Optional[str] = None

    def __init__(
        self,
        protein: Protein,
        *,
        ligands: list[Ligand] | None = None,
        **kwargs,
    ):
        """Initialize a Complex with either ligands or _ligands parameter"""
        self._ligands = ligands if ligands is not None else []
        self.protein = protein

        self.__post_init__()

    def __post_init__(self):
        """various post init tasks"""

        # assign references to the complex in the
        # various child classes
        self.docking = Docking(parent=self)
        self.abfe = ABFE(parent=self)
        self.rbfe = RBFE(parent=self)

    @property
    def ligands(self) -> list[Ligand]:
        """Get the current ligands"""
        return self._ligands

    @ligands.setter
    def ligands(self, new_ligands: list[Ligand]):
        """Set new ligands"""
        self._ligands = new_ligands

    @classmethod
    def from_dir(cls, directory: str) -> "Complex":
        """Initialize a Complex from a directory containing protein and ligand files.

        Args:
            directory (str): Directory containing ligand and protein files.

        The directory should contain:
        - Exactly one PDB file for the protein
        - One or more SDF files for the ligands. Each SDF file can contain one or more molecules.

        Returns:
            Complex: A new Complex instance initialized from the files in the directory.

        Raises:
            ValueError: If no PDB file is found or if multiple PDB files are found.
        """
        # Find all SDF files in the directory
        sdf_files = sorted(
            [
                os.path.join(directory, f)
                for f in os.listdir(directory)
                if f.lower().endswith(".sdf")
            ]
        )

        # Load all ligands from SDF files
        ligands = []
        for sdf_file in sdf_files:
            result = Ligand.from_sdf(sdf_file)
            if isinstance(result, list):
                ligands.extend(result)
            else:
                ligands.append(result)

        # Find PDB file
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
        protein = Protein.from_file(protein_file)

        # Create the Complex instance
        instance = cls(
            protein=protein,
            ligands=ligands,
        )

        return instance

    def _sync_protein_and_ligands(self) -> None:
        """Ensure that the protein and ligands are uploaded to Deep Origin

        Internal method. Do not use."""

        # get a list of all files in the entities directory
        remote_files = files_client.list_folder("entities", recursive=True)
        remote_files = list(remote_files.keys())

        files_to_upload = {}

        protein_path = "entities/proteins/" + os.path.basename(self.protein.file_path)
        if protein_path not in remote_files:
            files_to_upload[str(self.protein.file_path)] = protein_path

        for ligand in self.ligands:
            ligand_path = "entities/ligands/" + os.path.basename(ligand.file_path)
            if ligand_path not in remote_files:
                files_to_upload[str(ligand.file_path)] = ligand_path

        files_client.upload_files(files_to_upload)

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
            files = [ligand.file_path for ligand in self.ligands]

            if limit is not None:
                files = files[:limit]

            chem.show_molecules_in_sdf_files(files)
        else:
            from deeporigin.drug_discovery.structures.ligand import (
                show_ligands as _show_ligands,
            )

            if limit is not None:
                return _show_ligands(self.ligands[:limit])
            else:
                return _show_ligands(self.ligands)

    @beartype
    def get_result_files_for(
        self,
        *,
        tool,
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

        df = df[df[utils.COL_PROTEIN] == self.protein._do_id]

        if ligand_ids is not None:
            df = df[df[utils.COL_LIGAND1].isin(ligand_ids)]

        # this makes sure that we only retain rows
        # where there is a valid output being generated
        # rows where there is no output (for failed or
        # running jobs) are simply dropped
        df = df.dropna(subset=[utils.COL_RESULT])

        file_ids = list(df[utils.COL_RESULT])

        existing_files = os.listdir(utils.DATA_DIRS[tool])
        existing_files = ["_file:" + file for file in existing_files]
        missing_files = list(set(file_ids) - set(existing_files))

        if len(missing_files) > 0:
            print("Downloading result files. This can take a while...")
            api.download_files(
                file_ids=missing_files,
                use_file_names=False,
                save_to_dir=utils.DATA_DIRS[tool],
            )
            print("Done.")

        return [
            os.path.join(utils.DATA_DIRS[tool], file_id.replace("_file:", ""))
            for file_id in file_ids
        ]

    def get_csv_results_for(self, tool):
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

        ligand_ids = [ligand._do_id for ligand in self.ligands]

        df = df[df[utils.COL_PROTEIN] == self.protein._do_id]
        df = df[(df[utils.COL_LIGAND1].isin(ligand_ids))]

        # this makes sure that we only retain rows
        # where there is a valid output being generated
        # rows where there is no output (for failed or
        # running jobs) are simply dropped
        df = df.dropna(subset=[utils.COL_CSV_OUTPUT])
        file_ids = list(df[utils.COL_CSV_OUTPUT])

        existing_files = os.listdir(utils.DATA_DIRS[tool])
        existing_files = ["_file:" + file for file in existing_files]
        missing_files = list(set(file_ids) - set(existing_files))
        if len(missing_files) > 0:
            api.download_files(
                file_ids=missing_files,
                use_file_names=False,
                save_to_dir=utils.DATA_DIRS[tool],
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
        for file_id, ligand1_id, ligand2_id in zip(
            file_ids,
            ligand1_ids,
            ligand2_ids,
            strict=False,
        ):
            file_loc = os.path.join(
                utils.DATA_DIRS[tool], file_id.replace("_file:", "")
            )
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

            smiles_mapping = {ligand._do_id: ligand.smiles for ligand in self.ligands}

            if tool == "RBFE":
                df["Ligand1"] = ligand1_id
                df["Ligand2"] = ligand2_id

                df["SMILES1"] = df[utils.COL_LIGAND1].map(smiles_mapping)
                df["SMILES2"] = df[utils.COL_LIGAND2].map(smiles_mapping)
            elif tool == "ABFE":
                df["Ligand1"] = ligand1_id

                df["SMILES"] = df["Ligand1"].map(smiles_mapping)

            all_dfs.append(df)

        if len(all_dfs) == 0:
            # no data, so nothing to do
            return pd.DataFrame()

        df = pd.concat(all_dfs)

        return df
