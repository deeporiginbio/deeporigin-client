import os

import pytest

from deeporigin.drug_discovery.structures import Ligand
from deeporigin.exceptions import DeepOriginException

# Import shared test fixtures
from tests.utils_ligands import bad_ligands, ligands

base_path = os.path.join(os.path.dirname(__file__), "fixtures")


@pytest.mark.parametrize(
    "smiles,name,expected_atoms,equivalent_smiles",
    [
        ("C", "Methane", 1, None),  # Methane
        ("CC", "Ethane", 2, None),  # Ethane
        ("CCO", "Ethanol", 3, None),  # Ethanol
        ("c1ccccc1", "Benzene", 6, ["C1=CC=CC=C1"]),  # Benzene (aromatic notation)
        ("C1=CC=CC=C1", "Benzene_alt", 6, ["c1ccccc1"]),  # Benzene (Kekule notation)
    ],
)
def test_ligand_from_smiles(smiles, name, expected_atoms, equivalent_smiles):
    """Test that we can create a Ligand from a SMILES string using the from_smiles classmethod"""
    from rdkit import Chem

    # Create a ligand using the from_smiles method
    ligand = Ligand.from_smiles(
        smiles=smiles,
        name=name,
        save_to_file=False,
    )

    # Verify the ligand has either the exact SMILES string or an equivalent one
    if equivalent_smiles:
        assert ligand.smiles in [smiles] + equivalent_smiles, (
            f"SMILES {ligand.smiles} not equivalent to {smiles} or any of {equivalent_smiles}"
        )
    else:
        assert ligand.smiles == smiles

    # Verify the name was set correctly
    assert ligand.name == name

    # Verify that the file field is None
    assert ligand.file_path is None

    # Verify that the molecule was properly initialized
    assert ligand.mol is not None
    assert ligand.mol.GetNumAtoms() == expected_atoms

    # Verify that the molecule represents the same chemical structure
    input_mol = Chem.MolFromSmiles(smiles)
    assert Chem.MolToSmiles(input_mol) == Chem.MolToSmiles(ligand.mol)


def test_ligand_from_smiles_invalid():
    """Test that invalid SMILES raises DeepOriginException"""
    with pytest.raises(DeepOriginException, match=r"Cannot create"):
        Ligand.from_smiles(smiles="InvalidSMILES")


@pytest.mark.parametrize(
    "identifier,expected_atoms",
    [
        ("ATP", 31),  # Adenosine triphosphate
        ("ADP", 27),  # Adenosine diphosphate
        ("Oxotremorine", 15),  # Muscarinic acetylcholine receptor agonist
        ("Serotonin", 13),  # 5-hydroxytryptamine (5-HT)
    ],
)
def test_ligand_from_identifier(identifier, expected_atoms):
    """Test that we can create a Ligand from common biochemical identifiers using the from_identifier classmethod"""

    # Create a ligand using the from_identifier method
    ligand = Ligand.from_identifier(identifier=identifier)

    # Verify the name was set correctly
    assert ligand.name == identifier

    # Verify that the file field is None
    assert ligand.file_path is None

    # Verify that the molecule was properly initialized
    assert ligand.mol is not None
    assert ligand.mol.GetNumAtoms() == expected_atoms

    # Verify that the molecule has valid 3D coordinates
    assert ligand.mol.GetNumConformers() > 0
    coords = ligand.mol.GetConformer().GetPositions()
    assert coords.shape[0] == expected_atoms


def test_ligand_from_identifier_invalid():
    """Test that invalid identifier raises appropriate exception"""
    invalid_id = "InvalidMolecule123"
    with pytest.raises(ValueError) as exc_info:
        Ligand.from_identifier(identifier=invalid_id)
    assert str(exc_info.value).startswith(
        f"Error resolving SMILES string of {invalid_id}"
    )


def test_ligand_from_rdkit_mol():
    """Test that we can create a Ligand from an RDKit Mol object using the from_rdkit_mol classmethod"""
    from rdkit import Chem

    # Create test RDKit molecules
    mols = [
        Chem.MolFromSmiles("C"),  # Methane
        Chem.MolFromSmiles("CC"),  # Ethane
        Chem.MolFromSmiles("CCO"),  # Ethanol
        Chem.MolFromSmiles("c1ccccc1"),  # Benzene
    ]

    for mol in mols:
        # Create a ligand using the from_rdkit_mol method
        ligand = Ligand.from_rdkit_mol(mol, name="TestLigand")

        # Verify the ligand has the correct SMILES string
        assert ligand.smiles == Chem.MolToSmiles(mol)

        # Verify the name was set correctly
        assert ligand.name == "TestLigand"

        # Verify that the file field is None
        assert ligand.file_path is None

        # Verify that the molecule was properly initialized
        assert ligand.mol is not None
        assert ligand.mol.GetNumAtoms() == mol.GetNumAtoms()


def test_ligand_from_sdf():
    """Test that we can create a Ligand from an SDF file using the from_sdf classmethod"""
    # Use the brd-7.sdf file which contains exactly one ligand
    # Find the ligand entry for brd-7.sdf from the imported ligands variable
    brd7_ligand = next(ligand for ligand in ligands if "brd-7.sdf" in ligand["file"])
    sdf_file = brd7_ligand["file"]

    # Create a ligand using the from_sdf method
    ligand = Ligand.from_sdf(sdf_file)

    # Verify the ligand was created successfully
    assert isinstance(ligand, Ligand)
    assert ligand.mol is not None
    assert ligand.mol.GetNumAtoms() > 0

    # Verify that the file_path was set correctly
    assert ligand.file_path == sdf_file

    # Verify that the ligand has a name
    assert ligand.name is not None
    assert ligand.name != "Unknown_Ligand"

    # Verify that the ligand has SMILES
    assert ligand.smiles is not None

    # Verify that the ligand has properties (SDF files typically contain properties)
    assert isinstance(ligand.properties, dict)


def test_ligand_base64():
    brd7_ligand = next(ligand for ligand in ligands if "brd-7.sdf" in ligand["file"])
    sdf_file = brd7_ligand["file"]

    ligand = Ligand.from_sdf(sdf_file)

    b64 = ligand.to_base64()
    new_ligand = Ligand.from_base64(b64)

    assert new_ligand.smiles == ligand.smiles


@pytest.mark.parametrize("ligand", bad_ligands)
def test_ligand_errors(ligand):
    with pytest.raises(Exception):  # noqa: B017
        Ligand(
            file_path=ligand["file"],
            smiles=ligand["smiles"],
        )


@pytest.mark.parametrize("ligand", ligands)
def test_ligand(ligand):
    """Test that we can create Ligand instances from various sources"""
    n_ligands = ligand["n_ligands"]

    if n_ligands >= 1:
        return

    result = Ligand.from_sdf(ligand["file"])

    assert isinstance(result, Ligand)
    assert result.mol is not None
    assert result.mol.GetNumAtoms() > 0
    assert (
        result.file_path == ligand["file"]
    )  # Single ligand case should have file_path


def test_ligand_from_sdf_multiple_raises():
    """Test that Ligand.from_sdf raises DeepOriginException for multi-molecule SDF files."""
    with pytest.raises(
        DeepOriginException, match="must contain exactly one molecule, but found 8"
    ):
        Ligand.from_sdf(os.path.join(base_path, "ligands-brd-all.sdf"))
