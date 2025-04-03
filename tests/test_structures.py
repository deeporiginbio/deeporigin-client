from deeporigin.drug_discovery.structures import Ligand, ligands_to_dataframe
from deeporigin.drug_discovery import chemistry
import pytest
from deeporigin.exceptions import DeepOriginException
import os

# Import shared test fixtures
from tests.utils_ligands import ligands, bad_ligands

base_path = os.path.join(os.path.dirname(__file__), "fixtures")


def test_ligands_from_sdf_file():
    """test that we can make many ligands from a single SDF file with many molecules"""

    mols = chemistry.read_molecules_in_sdf_file(
        os.path.join(base_path, "ligands-brd-all.sdf")
    )

    for mol in mols:
        Ligand(**mol)


def test_ligand_from_smiles():
    """Test that we can create a Ligand from a SMILES string using the from_smiles classmethod"""

    # Example SMILES strings - using a couple simple molecules
    smiles_strings = ["C", "CC", "CCO", "c1ccccc1"]

    for smiles in smiles_strings:
        # Create a ligand using the from_smiles method
        ligand = Ligand.from_smiles(smiles)

        # Verify the ligand has the correct SMILES string
        assert ligand.smiles_string == smiles

        # Verify that the file field is None
        assert ligand.file is None

        # Verify that properties are also None initially
        assert ligand.properties is None


@pytest.mark.parametrize("ligand", bad_ligands)
def test_ligand_errors(ligand):
    with pytest.raises(DeepOriginException):
        Ligand(
            file=ligand["file"],
            smiles_string=ligand["smiles_string"],
        )


@pytest.mark.parametrize("ligand", ligands)
def test_ligand(
    ligand,
):
    n_ligands = ligand["n_ligands"]

    if n_ligands > 1:
        with pytest.raises(ValueError, match="Too many molecules."):
            ligand = Ligand(ligand["file"])
    else:
        ligand = Ligand(ligand["file"])


def test_ligand_from_csv(tmp_path):
    """Test that we can create Ligands from a CSV file using the from_csv classmethod"""
    import pandas as pd

    # Create a temporary CSV file with test data
    csv_path = tmp_path / "test_ligands.csv"

    # Create test data with SMILES and properties
    data = {
        "SMILES": [
            "C",
            "CC",
            "CCO",
            "c1ccccc1",
            "",
        ],  # Include an empty SMILES to test skipping
        "Name": ["Methane", "Ethane", "Ethanol", "Benzene", "Invalid"],
        "MolecularWeight": [16.04, 30.07, 46.07, 78.11, 0],
        "LogP": [-0.77, -0.18, -0.31, 2.13, 0],
    }

    # Write to CSV
    pd.DataFrame(data).to_csv(csv_path, index=False)

    # Test with all property columns
    ligands = Ligand.from_csv(
        file=csv_path,
        smiles_column="SMILES",
        properties_columns=["Name", "MolecularWeight", "LogP"],
    )

    # Should have 4 valid ligands (empty SMILES should be skipped)
    assert len(ligands) == 4

    # Check each ligand
    for i, ligand in enumerate(ligands):
        # Skip the empty row which should have been filtered out
        if i >= 4:
            continue

        # Check SMILES
        assert ligand.smiles_string == data["SMILES"][i]

        # Check properties
        assert ligand.properties is not None
        assert ligand.properties["Name"] == data["Name"][i]
        assert ligand.properties["MolecularWeight"] == data["MolecularWeight"][i]
        assert ligand.properties["LogP"] == data["LogP"][i]

    # Test with a non-existent column (should skip with warning)
    ligands = Ligand.from_csv(
        file=csv_path,
        smiles_column="SMILES",
        properties_columns=["Name", "NonExistentColumn"],
    )

    # Should still have 4 valid ligands
    assert len(ligands) == 4

    # All ligands should have Name property but not NonExistentColumn
    for ligand in ligands:
        assert "Name" in ligand.properties
        assert "NonExistentColumn" not in ligand.properties

    # Test with no property columns
    ligands = Ligand.from_csv(file=csv_path, smiles_column="SMILES")

    # Should still have 4 valid ligands
    assert len(ligands) == 4

    # All ligands should have None properties
    for ligand in ligands:
        assert ligand.properties is None

    # Test with invalid SMILES column name
    with pytest.raises(
        ValueError, match="SMILES column 'InvalidColumn' not found in CSV file"
    ):
        Ligand.from_csv(file=csv_path, smiles_column="InvalidColumn")
