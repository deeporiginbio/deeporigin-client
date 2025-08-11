"""
MoleculeModule


"""

import base64
import random
from typing import Optional
import warnings

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, SaltRemover, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

warnings.filterwarnings("ignore", category=UserWarning, module="rdkit")
RDLogger.DisableLog("rdApp.*")


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
    def from_smiles_or_name(
        cls,
        smiles: Optional[str] = None,
        name: Optional[str] = None,
        add_coords: bool = False,
        seed: Optional[int] = None,
    ):
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
            import requests

            try:
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/smiles/JSON"
                response = requests.get(url)
                data = response.json()
                smiles = data["PropertyTable"]["Properties"][0]["SMILES"]
            except Exception:
                raise ValueError(
                    f"Error resolving SMILES string of {name}. Pubchempy did not resolve SMILES"
                ) from None

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
