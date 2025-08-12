import os

import pytest

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
    assert first_ligand.mol.GetNumAtoms() > 0
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


def test_embed():
    """Test that we can minimize a LigandSet"""

    ligands = LigandSet.from_smiles(BRD_SMILES)
    ligands.embed()


def test_show():
    """Test that we can show a LigandSet"""

    ligands = LigandSet.from_smiles(BRD_SMILES)
    ligands.show()


def test_from_dir():
    """Test that we can create a LigandSet from a directory"""

    ligands = LigandSet.from_dir(DATA_DIR / "brd")
    assert len(ligands) == 8

    for ligand in ligands:
        assert ligand.file_path is not None
        assert os.path.exists(ligand.file_path)


def test_mcs():
    """Test that we can generate the MCS for a set of ligands"""

    from deeporigin.drug_discovery import BRD_DATA_DIR, LigandSet

    ligands = LigandSet.from_dir(BRD_DATA_DIR)
    ligands.mcs()


def test_compute_constraints():
    """Test that we can align a ligandset to a reference ligand"""

    from deeporigin.drug_discovery import BRD_DATA_DIR, LigandSet

    ligands = LigandSet.from_dir(BRD_DATA_DIR)
    ligands.compute_constraints(reference=ligands.ligands[0])


def test_random_sample():
    """Test the random_sample method of LigandSet"""

    # Create a test LigandSet
    test_smiles = ["CCO", "CCCO", "CCCC", "CCCCC", "CCCCCC"]
    ligands = LigandSet.from_smiles(test_smiles)

    # Test basic sampling
    sample = ligands.random_sample(3)
    assert isinstance(sample, LigandSet)
    assert len(sample) == 3
    assert len(sample.ligands) == 3

    # Test that original is unchanged
    assert len(ligands) == 5
    assert len(ligands.ligands) == 5

    # Test that sampled ligands are from original set
    for ligand in sample.ligands:
        assert ligand in ligands.ligands

    # Test edge cases
    sample_all = ligands.random_sample(5)
    assert len(sample_all) == 5
    assert set(sample_all.to_smiles()) == set(ligands.to_smiles())

    sample_one = ligands.random_sample(1)
    assert len(sample_one) == 1
    assert sample_one.ligands[0] in ligands.ligands


def test_random_sample_validation():
    """Test validation in random_sample method"""

    test_smiles = ["CCO", "CCCO", "CCCC"]
    ligands = LigandSet.from_smiles(test_smiles)

    # Test invalid n values
    with pytest.raises(ValueError, match="n must be at least 1"):
        ligands.random_sample(0)

    with pytest.raises(ValueError, match="n must be at least 1"):
        ligands.random_sample(-1)

    with pytest.raises(
        ValueError, match="Cannot sample 5 ligands from a set of 3 ligands"
    ):
        ligands.random_sample(5)

    with pytest.raises(
        ValueError, match="Cannot sample 10 ligands from a set of 3 ligands"
    ):
        ligands.random_sample(10)


def test_random_sample_deterministic():
    """Test that random_sample returns different results on multiple calls"""

    test_smiles = ["CCO", "CCCO", "CCCC", "CCCCC", "CCCCCC", "CCCCCCC"]
    ligands = LigandSet.from_smiles(test_smiles)

    # Sample multiple times and check we get different results
    samples = []
    for _ in range(5):
        sample = ligands.random_sample(3)
        samples.append(sample)

    # Check that at least some samples are different (this is probabilistic but should work)
    sample_smiles = [tuple(sorted(sample.to_smiles())) for sample in samples]
    unique_samples = set(sample_smiles)

    # With 6 ligands, sampling 3 should give us multiple unique combinations
    # This test might occasionally fail due to randomness, but it's very unlikely
    assert len(unique_samples) > 1, "Random sampling should produce different results"


# Test LigandSet functionality
def test_ligandset_operations():
    """Test basic LigandSet operations"""
    from deeporigin.drug_discovery.structures.ligand import LigandSet

    # Create test ligands
    ligand1 = Ligand.from_smiles("CCO", name="Ethanol")
    ligand2 = Ligand.from_smiles("CCCO", name="Propanol")

    # Test LigandSet creation
    ligandset = LigandSet(ligands=[ligand1, ligand2])
    assert len(ligandset) == 2

    # Test iteration
    for ligand in ligandset:
        assert isinstance(ligand, Ligand)

    # Test indexing
    assert ligandset[0] == ligand1
    assert ligandset[1] == ligand2

    # Test containment
    assert ligand1 in ligandset
    assert ligand2 in ligandset


def test_ligandset_addition():
    """Test LigandSet addition operations"""
    from deeporigin.drug_discovery.structures.ligand import LigandSet

    ligand1 = Ligand.from_smiles("CCO", name="Ethanol")
    ligand2 = Ligand.from_smiles("CCCO", name="Propanol")

    set1 = LigandSet(ligands=[ligand1])
    set2 = LigandSet(ligands=[ligand2])

    # Test LigandSet + LigandSet
    combined = set1 + set2
    assert len(combined) == 2

    # Test LigandSet + Ligand
    combined = set1 + ligand2
    assert len(combined) == 2

    # Test Ligand + LigandSet
    combined = ligand2 + set1
    assert len(combined) == 2


def test_ligandset_from_smiles():
    """Test LigandSet creation from SMILES"""
    from deeporigin.drug_discovery.structures.ligand import LigandSet

    smiles_list = ["CCO", "CCCO", "CCCCO"]
    ligandset = LigandSet.from_smiles(smiles_list)

    assert len(ligandset) == 3
    assert all(isinstance(ligand, Ligand) for ligand in ligandset)
    assert ligandset[0].smiles == "CCO"
    assert ligandset[1].smiles == "CCCO"
    assert ligandset[2].smiles == "CCCCO"


def test_ligandset_to_dataframe():
    """Test LigandSet to DataFrame conversion"""
    from deeporigin.drug_discovery.structures.ligand import LigandSet

    ligand1 = Ligand.from_smiles("CCO", name="Ethanol")
    ligand2 = Ligand.from_smiles("CCCO", name="Propanol")

    ligandset = LigandSet(ligands=[ligand1, ligand2])

    # Add properties
    ligand1.set_property("logP", 0.32)
    ligand2.set_property("logP", 0.88)

    df = ligandset.to_dataframe()
    assert len(df) == 2
    assert "Ligand" in df.columns
    assert "logP" in df.columns


def test_ligandset_indexing_and_slicing():
    """Test LigandSet indexing and slicing behavior"""
    from deeporigin.drug_discovery.structures.ligand import LigandSet

    # Create test ligands
    test_smiles = ["CCO", "CCCO", "CCCC", "CCCCC"]
    ligandset = LigandSet.from_smiles(test_smiles)

    # Test single indexing (should return Ligand)
    single_ligand = ligandset[0]
    assert isinstance(single_ligand, Ligand)
    assert single_ligand.smiles == "CCO"

    # Test slicing (should return LigandSet)
    subset = ligandset[1:3]
    assert isinstance(subset, LigandSet)
    assert len(subset) == 2
    assert subset[0].smiles == "CCCO"
    assert subset[1].smiles == "CCCC"

    # Test slice from beginning
    start_slice = ligandset[:2]
    assert isinstance(start_slice, LigandSet)
    assert len(start_slice) == 2

    # Test slice to end
    end_slice = ligandset[2:]
    assert isinstance(end_slice, LigandSet)
    assert len(end_slice) == 2

    # Test that original is unchanged
    assert len(ligandset) == 4
