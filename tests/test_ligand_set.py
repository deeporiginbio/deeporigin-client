import os

import pytest
from rdkit import Chem

from deeporigin.drug_discovery import DATA_DIR
from deeporigin.drug_discovery.structures.ligand import Ligand, LigandSet
from deeporigin.exceptions import DeepOriginException

# Import shared test fixtures


SDF_TEST_CASES = [
    (DATA_DIR / "ligands" / "ligands-brd-all.sdf", 8),
    (DATA_DIR / "ligands" / "42-ligands.sdf", 42),
]

BRD_SMILES = {
    "C/C=C/Cn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
    "C=CCCn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
    "C=CCn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
    "CCCCn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
    "CCCn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
    "CCn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
    "CN(C)C(=O)c1cccc(-c2cn(C)c(=O)c3[nH]ccc23)c1",
    "COCCn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
}


@pytest.mark.parametrize("filename,expected_count", SDF_TEST_CASES)
def test_ligand_set_from_sdf_file(filename, expected_count):
    """Test that we can make many ligands from a single SDF file with many molecules"""
    ligands = LigandSet.from_sdf(filename)
    assert len(ligands.ligands) == expected_count, f"Expected {expected_count} ligands"
    for ligand in ligands.ligands:
        assert isinstance(ligand, Ligand), "Expected a Ligand object"


def test_ligand_set_from_csv():
    """Test that we can create Ligands from a CSV file using the from_csv classmethod"""

    # Get the path to the test CSV file
    csv_path = DATA_DIR / "ligands" / "ligands.csv"

    # Create ligands from the CSV file
    ligands = LigandSet.from_csv(str(csv_path), smiles_column="SMILES")

    # Verify we got the expected number of ligands
    assert len(ligands.ligands) == 30  # Total number of valid SMILES in the file

    # Check a few properties of the first ligand
    first_ligand = ligands.ligands[0]
    assert isinstance(first_ligand, Ligand)
    assert first_ligand.mol is not None
    assert first_ligand.mol.m.GetNumAtoms() > 0
    assert first_ligand.file_path is None

    # Verify properties were correctly loaded
    props = ["score", "binding_energy", "pose_score"]
    df = ligands.to_dataframe()
    for prop in props:
        assert prop in df.columns
        assert df[prop].notna().all()

    # Test with invalid SMILES column
    with pytest.raises(
        DeepOriginException, match="Column 'invalid' not found in CSV file"
    ):
        LigandSet.from_csv(str(csv_path), smiles_column="invalid")

    # Test with non-existent file
    with pytest.raises(FileNotFoundError):
        LigandSet.from_csv("nonexistent.csv")


@pytest.mark.parametrize("filename,expected_count", SDF_TEST_CASES)
def test_sdf_roundtrip(filename, expected_count):
    """Test that we can roundtrip a LigandSet to an SDF file and back for all SDF_TEST_CASES"""

    ligands = LigandSet.from_sdf(filename)
    assert len(ligands.ligands) == expected_count
    sdf_path = ligands.to_sdf()
    assert os.path.exists(sdf_path)

    new_ligands = LigandSet.from_sdf(sdf_path)
    assert len(new_ligands.ligands) == len(ligands.ligands)
    assert set(ligands.to_smiles()) == set(new_ligands.to_smiles()), (
        "SMILES strings should be the same"
    )

    os.unlink(sdf_path)


def test_to_smiles():
    """Test that we can convert a LigandSet to SMILES strings"""

    ligands = LigandSet.from_sdf(DATA_DIR / "ligands" / "ligands-brd-all.sdf")

    assert set(ligands.to_smiles()) == BRD_SMILES, "SMILES strings should be the same"


def test_from_smiles():
    """Test that we can create a LigandSet from a list of SMILES strings."""

    ligands = LigandSet.from_smiles(BRD_SMILES)
    assert isinstance(ligands, LigandSet)
    assert len(ligands) == len(BRD_SMILES)

    # Check that all SMILES are present (order-insensitive)
    assert set(ligands.to_smiles()) == set(BRD_SMILES)
    for ligand in ligands:
        assert isinstance(ligand, Ligand)


def test_minimize():
    """Test that we can minimize a LigandSet"""

    ligands = LigandSet.from_smiles(BRD_SMILES)
    ligands.minimize()


def test_show():
    """Test that we can show a LigandSet"""

    ligands = LigandSet.from_smiles(BRD_SMILES)
    ligands.show()


def test_from_dir():
    """Test that we can create a LigandSet from a directory"""

    ligands = LigandSet.from_dir(DATA_DIR / "brd")
    assert len(ligands) == 8


def test_mcs():
    """Test that LigandSet.mcs() returns a valid SMARTS string for a simple set of ligands."""

    ligands = LigandSet.from_smiles(BRD_SMILES)
    smarts = ligands.mcs()
    assert isinstance(smarts, str)
    assert smarts != ""
    # Check that the SMARTS string can be parsed by RDKit
    mcs_mol = Chem.MolFromSmarts(smarts)
    assert mcs_mol is not None

    # The MCS should have at least 6 atoms
    assert mcs_mol.GetNumAtoms() >= 6
