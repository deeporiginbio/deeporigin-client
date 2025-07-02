"""Module to support Drug Discovery workflows using Deep Origin"""

from dataclasses import dataclass, field
import os
from typing import Optional

from beartype import beartype

from deeporigin.drug_discovery.abfe import ABFE
from deeporigin.drug_discovery.docking import Docking
from deeporigin.drug_discovery.rbfe import RBFE
from deeporigin.drug_discovery.structures import Ligand, LigandSet, Protein
from deeporigin.exceptions import DeepOriginException
from deeporigin.files import FilesClient
from deeporigin.platform.utils import PlatformClients


@dataclass
class Complex:
    """class to represent a set of a protein and 1 or many ligands"""

    protein: Protein

    # Using a private attribute for ligands with the property decorator below
    ligands: LigandSet = field(default_factory=LigandSet)
    _platform_clients: Optional[PlatformClients] = None
    _prepared_systems: dict[str, str] = field(default_factory=dict, repr=False)

    def __post_init__(self):
        """various post init tasks"""

        # assign references to the complex in the
        # various child classes
        self.docking = Docking(parent=self)
        self.abfe = ABFE(parent=self)
        self.rbfe = RBFE(parent=self)

        self._prepared_systems = {}

    @classmethod
    def from_dir(
        cls,
        directory: str,
        *,
        _platform_clients: Optional[PlatformClients] = None,
    ) -> "Complex":
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
            raise DeepOriginException(
                f"Expected exactly one PDB file in the directory, but found {len(pdb_files)}: {pdb_files}",
                title="Complex.from_dir expects a single PDB file",
            ) from None
        protein_file = pdb_files[0]
        protein = Protein.from_file(protein_file)

        # Create the Complex instance
        instance = cls(
            protein=protein,
            ligands=LigandSet(ligands=ligands),
            _platform_clients=_platform_clients,
        )

        return instance

    @beartype
    def prepare(
        self,
        ligand: Ligand,
        *,
        padding: float = 1.0,
        keep_waters: bool = False,
        is_lig_protonated: bool = True,
        is_protein_protonated: bool = True,
        use_cache: bool = True,
    ) -> None:
        """run system preparation on the protein and one ligand from the Complex

        Args:
            ligand (Ligand): The ligand to prepare.
            padding (float, optional): Padding to add around the system. Defaults to 1.0.
            keep_waters (bool, optional): Whether to keep water molecules. Defaults to False.
            is_lig_protonated (bool, optional): Whether the ligand is already protonated. Defaults to True.
            is_protein_protonated (bool, optional): Whether the protein is already protonated. Defaults to True.
        """
        from deeporigin.functions.sysprep import sysprep

        if ligand.file_path is None:
            ligand_path = ligand.to_sdf()
        else:
            ligand_path = ligand.file_path

        # run sysprep on the ligand
        complex_path = sysprep(
            protein_path=self.protein.file_path,
            padding=padding,
            ligand_path=ligand_path,
            keep_waters=keep_waters,
            is_lig_protonated=is_lig_protonated,
            is_protein_protonated=is_protein_protonated,
            use_cache=use_cache,
        )

        # set this complex path as the prepared system
        self._prepared_systems[ligand.name] = complex_path

        # show it
        Protein.from_file(complex_path).show()

    def _sync_protein_and_ligands(self) -> None:
        """Ensure that the protein and ligands are uploaded to Deep Origin

        Internal method. Do not use."""

        # the reason we are uploading here manually, instead of using ligand.upload()
        # and protein.upload() is so that we can make one call to upload_files, instead
        # of several

        if self._platform_clients is None:
            files_client = FilesClient()
        else:
            files_client = self._platform_clients.FilesApi

        # get a list of all files in the entities directory
        remote_files = files_client.list_folder("entities", recursive=True)
        remote_files = list(remote_files.keys())

        files_to_upload = {}

        protein_path = self.protein._remote_path_base + os.path.basename(
            self.protein.file_path
        )
        self.protein._remote_path = protein_path
        if protein_path not in remote_files:
            files_to_upload[str(self.protein.file_path)] = protein_path

        for ligand in self.ligands:
            if ligand.file_path is None:
                # this ligand isn't being backed by a file, so we can't upload it
                continue
            ligand_path = ligand._remote_path_base + os.path.basename(ligand.file_path)
            ligand._remote_path = ligand_path
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
