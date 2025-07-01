import pytest

from deeporigin.drug_discovery import DATA_DIR
from deeporigin.drug_discovery.structures.ligand import Ligand, LigandSet
from deeporigin.exceptions import DeepOriginException

# Import shared test fixtures


SDF_TEST_CASES = [
    (DATA_DIR / "ligands" / "ligands-brd-all.sdf", 8),
    (DATA_DIR / "ligands" / "42-ligands.sdf", 42),
]


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
