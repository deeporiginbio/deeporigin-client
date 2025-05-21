"""Module to support Drug Discovery workflows using Deep Origin"""

from dataclasses import dataclass, field
import os
from typing import Optional

from beartype import beartype
import do_sdk_platform

from deeporigin.drug_discovery import chemistry as chem
from deeporigin.drug_discovery.abfe import ABFE
from deeporigin.drug_discovery.docking import Docking
from deeporigin.drug_discovery.rbfe import RBFE
from deeporigin.drug_discovery.structures import Ligand, Protein
from deeporigin.files import FilesClient


@dataclass
class Complex:
    """class to represent a set of a protein and 1 or many ligands"""

    # Using a private attribute for ligands with the property decorator below
    protein: Protein
    _ligands: list[Ligand] = field(default_factory=list, repr=False)

    # don't set these directly. they are set in the constructor, if needed
    _files_client: Optional[FilesClient] = None
    _tools_client: Optional[do_sdk_platform.api.tools_api.ToolsApi] = None
    _organizations_client: Optional[
        do_sdk_platform.api.organizations_api.OrganizationsApi
    ] = None

    _organization_id: Optional[str] = None

    def __init__(
        self,
        protein: Protein,
        *,
        ligands: list[Ligand] | None = None,
        organization_id: Optional[str] = None,
        token: Optional[str] = None,
        base_url: Optional[str] = None,
    ):
        """Initialize a Complex with either ligands or _ligands parameter"""
        self._ligands = ligands if ligands is not None else []
        self.protein = protein

        if token is not None and base_url is not None and organization_id is not None:
            # we are being passed a token, base_url, and organization_id
            # this allows us to create clients for platform and files API,
            # and these clients will be used for all the operations in the Complex

            from deeporigin.platform.utils import _get_api_client

            self._files_client = FilesClient(
                token=token,
                base_url=base_url,
                organization_id=organization_id,
            )

            api_endpoint = base_url + "/api/"

            self._tools_client = _get_api_client(
                api_name="ToolsApi",
                token=token,
                api_endpoint=api_endpoint,
            )

            self._organizations_client = _get_api_client(
                api_name="OrganizationsApi",
                token=token,
                api_endpoint=api_endpoint,
            )

            self._organization_id = organization_id

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
    def from_dir(
        cls,
        directory: str,
        *,
        organization_id: Optional[str] = None,
        token: Optional[str] = None,
        base_url: Optional[str] = None,
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
            raise ValueError(
                f"Expected exactly one PDB file in the directory, but found {len(pdb_files)}."
            )
        protein_file = pdb_files[0]
        protein = Protein.from_file(protein_file)

        # Create the Complex instance
        instance = cls(
            protein=protein,
            ligands=ligands,
            organization_id=organization_id,
            token=token,
            base_url=base_url,
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
    ) -> None:
        """run system prepartion on the Complex

        Args:
            ligand (Ligand): The ligand to prepare.
            padding (float, optional): Padding to add around the system. Defaults to 1.0.
            keep_waters (bool, optional): Whether to keep water molecules. Defaults to False.
            is_lig_protonated (bool, optional): Whether the ligand is already protonated. Defaults to True.
            is_protein_protonated (bool, optional): Whether the protein is already protonated. Defaults to True.
        """
        from deeporigin.functions.sysprep import sysprep

        # run sysprep on the ligand
        complex_path = sysprep(
            protein_path=self.protein.file_path,
            ligand_path=ligand.file_path,
            padding=padding,
            keep_waters=keep_waters,
            is_lig_protonated=is_lig_protonated,
            is_protein_protonated=is_protein_protonated,
        )

        # show it
        Protein.from_file(complex_path).show()

    def _sync_protein_and_ligands(self) -> None:
        """Ensure that the protein and ligands are uploaded to Deep Origin

        Internal method. Do not use."""

        # get a list of all files in the entities directory
        from deeporigin.files import FilesClient

        if self._files_client is None:
            files_client = FilesClient()
        else:
            files_client = self._files_client

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

        files_client.upload_files(files_to_upload, client=self._files_client)

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
