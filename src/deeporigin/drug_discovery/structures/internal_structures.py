"""
MoleculeModule

Description: This module encapsulates the Molecule class and related utilities for managing and manipulating
small molecule structures in computational biology and cheminformatics workflows. It provides functionalities
to initialize molecules from various sources, handle molecular properties, perform format conversions, and
visualize molecular structures. The primary features and methods include:

- **FileFormat Enum**: Defines supported file formats (`MOL`, `MOL2`, `PDB`, `PDBQT`, `XYZ`, `SDF`) for
  molecular data, ensuring consistent handling across different input types.

- **Molecule Class**: Represents a small molecule, managing its RDKit molecule object (`mol_rdk`), chemical
  properties (e.g., formula, SMILES), coordinate generation, and conformer management. Provides methods for
  molecule manipulation such as adding hydrogens, embedding coordinates, and assigning bond orders from SMILES.

- **Initialization Methods**:
  - **from_smiles_or_name**: Creates a Molecule instance from a SMILES string or a compound name using PubChem API.
  - **mol_from_file**: Loads a Molecule from a specified file based on its format, supporting error handling for
    unsupported formats or parsing failures.
  - **mol_from_smiles**: Generates a Molecule from a SMILES string with optional sanitization.
  - **mol_from_block**: Converts block content (e.g., MOL, PDB) into a Molecule instance by writing to a temporary
    file and parsing it.

- **Property Management**: Calculates and stores molecular properties such as the number of heavy atoms (`hac`),
  molecular formula, and atomic species. Facilitates property assignment and retrieval for further analysis.

- **Conformer Management**: Handles the generation and selection of molecular conformers, allowing for the
  computation of 2D or 3D coordinates and the assignment of unique conformer IDs.

- **Utility Functions**:
  - **convert_to_sdf**: Converts molecular block content to SDF format, ensuring compatibility with various tools.

- **Visualization**: Integrates with the `MoleculeViewer` to render interactive visualizations of molecules within
  Jupyter Notebooks, facilitating intuitive analysis and presentation of molecular structures.

Dependencies: Utilizes libraries such as RDKit for cheminformatics, PubChemPy for API interactions, NumPy and Pandas
for data manipulation, Termcolor for colored terminal output, and integrates with the `MoleculeViewer` from the
`deeporigin_molstar` package for visualization purposes.

Usage Example:
```python
from molecule_module import Molecule, mol_from_smiles, mol_from_file

# Initialize a Molecule instance from a SMILES string
molecule = Molecule.from_smiles_or_name(smiles="CCO", name="Ethanol", add_coords=True)

# Initialize a Molecule instance from a file
molecule = mol_from_file("pdb", "/path/to/molecule.pdb")

# Calculate molecular properties
print(f"Formula: {molecule.formula}")
print(f"SMILES: {molecule.smiles}")

# Add hydrogens and embed coordinates
molecule.add_hydrogens(add_coords=True)
molecule.embed(add_hs=True, seed=42)

# Assign bond orders from a template SMILES string
molecule.assign_bond_order_from_smiles(smiles="CCO")

# Visualize the molecule
molecule.visualize()
```
"""

import base64
from enum import Enum
import random
import tempfile
import warnings

from beartype import beartype
from pubchempy import get_compounds
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, SaltRemover, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

# from ..utilities.conversions import convert_file


warnings.filterwarnings("ignore", category=UserWarning, module="rdkit")
RDLogger.DisableLog("rdApp.*")


class FileFormat(Enum):
    """
    Enumeration of supported file formats for molecular data.

    This enum class defines the various file formats that can be used to represent molecular structures,
    including MOL, MOL2, PDB, PDBQT, XYZ, and SDF formats.
    """

    MOL = "mol"
    MOL2 = "mol2"

    PDB = "pdb"
    PDBQT = "pdbqt"

    XYZ = "xyz"
    SDF = "sdf"


RDKIT_SUPPORTED_INPUT_TYPES = [
    FileFormat.MOL.value,
    FileFormat.MOL2.value,
    FileFormat.PDB.value,
    FileFormat.XYZ.value,
    FileFormat.SDF.value,
]


class Molecule:
    """
    A class representing a small molecule with various properties and manipulation capabilities.

    This class provides a comprehensive interface for working with molecular structures, including
    initialization from different sources, property calculation, coordinate generation, and visualization.

    Attributes:
        m (rdkit.Chem.Mol): The RDKit molecule object
        n (int): Number of atoms in the molecule
        formula (str): Molecular formula
        smiles (str): Canonical SMILES representation
        name (str): Name of the molecule
        contains_boron (bool): Whether the molecule contains boron atoms
    """

    def __init__(self, mol_rdk, name=None, add_coords=False, seed=None):
        """
        Initialize a Molecule instance.

        Args:
            mol_rdk (rdkit.Chem.Mol): RDKit molecule object
            name (str, optional): Name of the molecule
            add_coords (bool, optional): Whether to generate 3D coordinates
            seed (int, optional): Random seed for coordinate generation
        """
        self.m = self.process_mol(mol_rdk)
        self.n = mol_rdk.GetNumAtoms()
        self.formula = rdMolDescriptors.CalcMolFormula(self.m)
        self.smiles = Chem.MolToSmiles(Chem.RemoveHs(self.m), canonical=True)

        self.name = name
        self.contains_boron = any(atom.GetSymbol() == "B" for atom in self.m.GetAtoms())

        if self.name is not None:
            self.m.SetProp("_Name", name)

        if not self.m.GetConformers():
            AllChem.Compute2DCoords(self.m)

        self.set_conformer_id()
        if add_coords:
            seed = random.randint(0, 2**16 - 1) if seed is None else seed
            self.embed(add_hs=False, seed=seed)

        self.m.SetProp("initial_smiles", Chem.MolToSmiles(Chem.RemoveHs(mol_rdk)))

    def process_mol(self, mol):
        """
        Clean the ligand molecule by removing hydrogens and sanitizing the structure.

        Args:
            mol (rdkit.Chem.Mol): Input molecule to process

        Returns:
            rdkit.Chem.Mol: Processed molecule

        Raises:
            ValueError: If salt removal or kekulization fails
        """
        remover = SaltRemover.SaltRemover()

        stripped_mol = remover.StripMol(mol)
        if stripped_mol is None:
            raise ValueError("Salt removal failed.")

        try:
            Chem.Kekulize(stripped_mol, clearAromaticFlags=False)
        except Chem.KekulizeException as e:
            raise ValueError("Kekulization failed.") from e

        return stripped_mol

    def molblock(self):
        """
        Generate a MOL block representation of the molecule.

        Returns:
            str: MOL block string
        """
        return Chem.MolToMolBlock(self.m)

    def species(self):
        """
        Get the atomic symbols of all atoms in the molecule.

        Returns:
            list: List of atomic symbols
        """
        return [a.GetSymbol() for a in self.m.GetAtoms()]

    def coords(self, i: int = 0):
        """
        Get the coordinates of atoms in a specific conformer.

        Args:
            i (int): Conformer index

        Returns:
            numpy.ndarray: Array of atomic coordinates
        """
        conf = self.conformer(i)
        return conf.GetPositions()

    @classmethod
    def from_smiles_or_name(cls, smiles=None, name=None, add_coords=False, seed=None):
        """
        Create a Molecule instance from a SMILES string or compound name.

        Args:
            smiles (str, optional): SMILES string
            name (str, optional): Compound name
            add_coords (bool, optional): Whether to generate 3D coordinates
            seed (int, optional): Random seed for coordinate generation

        Returns:
            Molecule: New Molecule instance

        Raises:
            ValueError: If no compound is found for the given name
            AssertionError: If neither smiles nor name is provided
        """
        if smiles is None and name is not None:
            compounds = get_compounds(name, "name")
            if not compounds:
                raise ValueError(f"No compound found for identifier: {name}")

            compound = compounds[0]
            smiles = compound.isomeric_smiles

        assert smiles is not None

        mol_rdk = Chem.MolFromSmiles(smiles)

        return cls(mol_rdk, name=name, add_coords=add_coords, seed=seed)

    def copy(self):
        """
        Create a deep copy of the molecule.

        Returns:
            Molecule: Copy of the current molecule
        """
        return Molecule(Chem.Mol(self.m), name=self.name)

    def _draw(self):
        """
        Generate an HTML image representation of the molecule.

        Returns:
            str: HTML img tag with base64-encoded PNG
        """
        try:
            mol = self.copy()
            m = Chem.RemoveHs(mol.m)
            AllChem.Compute2DCoords(m)

            drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
            drawer.DrawMolecule(m)
            drawer.FinishDrawing()

            b64_encoded_png = base64.b64encode(drawer.GetDrawingText())
            html_img = (
                '<img src="data:image/png;base64,'
                + b64_encoded_png.decode("utf-8")
                + '">'
            )

            return html_img
        except Exception:
            return ""

    def draw(self):
        """
        Generate a 2D representation of the molecule.

        Returns:
            rdkit.Chem.Mol: Molecule with 2D coordinates
        """
        try:
            mol = self.copy()
            AllChem.Compute2DCoords(mol.m)

            return Chem.RemoveHs(mol.m)
        except Exception:
            return self.m

    def conformer(self, i=0):
        """
        Get a specific conformer of the molecule.

        Args:
            i (int): Conformer index

        Returns:
            rdkit.Chem.Conformer: The requested conformer
        """
        return self.m.GetConformer(i)

    def conformer_id(self):
        """
        Get the ID of the current conformer.

        Returns:
            int: Conformer ID
        """
        return self.m.GetConformer().GetId()

    def set_conformer_id(self, i=0):
        """
        Set the ID of the current conformer.

        Args:
            i (int): New conformer ID
        """
        self.m.GetConformer().SetId(i)

    def embed(self, add_hs=True, seed=-1):
        """
        Generate 3D coordinates for the molecule.

        Args:
            add_hs (bool): Whether to add hydrogens
            seed (int): Random seed for coordinate generation
        """
        if add_hs:
            self.add_hydrogens()

        AllChem.EmbedMolecule(self.m, randomSeed=seed)
        self.set_conformer_id(0)

    def add_hydrogens(self, add_coords=True):
        """
        Add hydrogens to the molecule.

        Args:
            add_coords (bool): Whether to generate coordinates for added hydrogens
        """
        self.m = Chem.AddHs(self.m, addCoords=add_coords)

    def assign_bond_order_from_smiles(self, smiles):
        """
        Assign bond orders from a template SMILES string.

        Args:
            smiles (str): Template SMILES string
        """
        template = Chem.MolFromSmiles(smiles)
        try:
            self.m = AllChem.AssignBondOrdersFromTemplate(template, self.m)
        except Exception:
            self.m = Chem.RemoveHs(self.m)
            self.m = AllChem.AssignBondOrdersFromTemplate(template, self.m)

        self.smiles = smiles


@beartype
def mol_from_file(
    file_type: str,
    file_path: str,
    sanitize: bool = True,
    remove_hs: bool = False,
) -> Molecule:
    """
    Create a Molecule instance from a file.

    Args:
        file_type (str): Type of the input file (must be in RDKIT_SUPPORTED_INPUT_TYPES)
        file_path (str): Path to the input file
        sanitize (bool): Whether to sanitize the molecule
        remove_hs (bool): Whether to remove hydrogens

    Returns:
        Molecule: New Molecule instance

    Raises:
        ValueError: If the file format is invalid or parsing fails
        NotImplementedError: If the file type is not supported
    """
    if file_type in RDKIT_SUPPORTED_INPUT_TYPES:
        mol_rdk = None

        if file_type == FileFormat.MOL.value:
            mol_rdk = Chem.MolFromMolFile(file_path, sanitize, remove_hs)
        elif file_type == FileFormat.MOL2.value:
            mol_rdk = Chem.MolFromMol2File(file_path, sanitize, remove_hs)
        elif file_type == FileFormat.PDB.value:
            mol_rdk = Chem.MolFromPDBFile(file_path, sanitize, remove_hs)
        elif file_type == FileFormat.XYZ.value:
            mol_rdk = Chem.MolFromXYZFile(file_path)
            if sanitize:
                Chem.SanitizeMol(mol_rdk)
            if remove_hs:
                mol_rdk = Chem.RemoveHs(mol_rdk)
        elif file_type == FileFormat.SDF.value:
            mol_rdk = next(iter(Chem.SDMolSupplier(file_path, sanitize, remove_hs)))

        if mol_rdk is None:
            raise ValueError(
                "Invalid file format or file path or failed to sanitize the molecule"
            )

        return Molecule(mol_rdk)
    else:
        raise NotImplementedError(
            f"Conversions of file types not supported yet. File format {file_type} is not supported."
        )


@beartype
def mol_from_smiles(smiles: str, sanitize: bool = True) -> Molecule:
    """
    Create a Molecule instance from a SMILES string.

    Args:
        smiles (str): SMILES string
        sanitize (bool): Whether to sanitize the molecule

    Returns:
        Molecule: New Molecule instance

    Raises:
        ValueError: If the SMILES string is invalid
    """
    mol_rdk = Chem.MolFromSmiles(smiles, sanitize=sanitize)
    if mol_rdk is None:
        raise ValueError("Invalid SMILES string")

    return Molecule(mol_rdk, add_coords=False)


@beartype
def mol_from_block(
    block_type: str,
    block: str,
    sanitize: bool = True,
    remove_hs: bool = False,
) -> Molecule:
    """
    Create a Molecule instance from a block of text.

    Args:
        block_type (str): Type of the input block
        block (str): Text block containing molecular data
        sanitize (bool): Whether to sanitize the molecule
        remove_hs (bool): Whether to remove hydrogens

    Returns:
        Molecule: New Molecule instance
    """
    with tempfile.TemporaryFile(mode="w+") as temp_file:
        temp_file.write(block)
        temp_file.seek(0)  # Reset file pointer to beginning
        return mol_from_file(
            block_type,
            temp_file.name,
            sanitize=sanitize,
            remove_hs=remove_hs,
        )
