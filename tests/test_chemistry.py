"""tests functions in the chemistry module"""

import os

import pytest
from deeporigin.drug_discovery import chemistry
from deeporigin.exceptions import DeepOriginException

base_path = os.path.join(os.path.dirname(__file__), "fixtures")

ligands = [
    {
        "file": os.path.join(base_path, "42-ligands.sdf"),
        "n_ligands": 42,
        "name_by_property": "Compound",
    },
    {
        "file": os.path.join(base_path, "ligands-brd-all.sdf"),
        "n_ligands": 8,
        "name_by_property": "_Name",
    },
    {
        "file": os.path.join(base_path, "brd-7.sdf"),
        "n_ligands": 1,
        "name_by_property": "_Name",
    },
]


bad_ligands = [
    {
        "file": "missing.sdf",
        "smiles_string": None,
    },
    {
        "file": None,
        "smiles_string": None,
    },
]


@pytest.mark.parametrize("ligand", ligands)
def test_count_molecules_in_sdf_file(
    tmp_path,
    ligand,
):
    assert chemistry.count_molecules_in_sdf_file(ligand["file"]) == ligand["n_ligands"]


@pytest.mark.parametrize("ligand", ligands)
def test_split_sdf_file(
    tmp_path,
    ligand,
):
    """
    Test that split_sdf_using_names correctly splits the ligands SDF file
    into separate SDF files, and that the output is cleaned up (by pytest)
    after the test completes.
    """

    # Create an output directory within the pytest temp directory
    output_dir = tmp_path / "split_ligands"
    output_dir.mkdir(exist_ok=True)

    # Call the function to be tested
    sdf_file_names = chemistry.split_sdf_file(
        input_sdf_path=str(ligand["file"]),
        output_prefix="testLig",
        output_dir=str(output_dir),
        name_by_property=ligand["name_by_property"],
    )

    # Check that at least one output file was created
    sdf_files = list(output_dir.glob("*.sdf"))
    assert len(sdf_files) > 0, "No SDF files were created by the splitting function."

    assert len(sdf_file_names) == ligand["n_ligands"], (
        "The number of SDF files is incorrect."
    )

    for sdf_file in sdf_files:
        n_mol = chemistry.count_molecules_in_sdf_file(sdf_file)
        assert n_mol == 1, "The SDF file contains more than one molecule."


@pytest.mark.parametrize("ligand", ligands)
def test_ligand(
    ligand,
):
    n_ligands = ligand["n_ligands"]

    if n_ligands > 1:
        with pytest.raises(ValueError, match="Too many molecules."):
            ligand = chemistry.Ligand(ligand["file"])
    else:
        ligand = chemistry.Ligand(ligand["file"])


@pytest.mark.parametrize("ligand", bad_ligands)
def test_ligand_errors(ligand):
    with pytest.raises(DeepOriginException):
        chemistry.Ligand(
            file=ligand["file"],
            smiles_string=ligand["smiles_string"],
        )


def test_ligands_from_sdf_file():
    """test that we can make many ligands from a single SDF file with many molecules"""

    mols = chemistry.read_molecules_in_sdf_file(
        os.path.join(base_path, "ligands-brd-all.sdf")
    )

    for mol in mols:
        chemistry.Ligand(**mol)


def test_ligand_from_smiles():
    """Test that we can create a Ligand from a SMILES string using the from_smiles classmethod"""

    # Example SMILES strings - using a couple simple molecules
    smiles_strings = ["C", "CC", "CCO", "c1ccccc1"]

    for smiles in smiles_strings:
        # Create a ligand using the from_smiles method
        ligand = chemistry.Ligand.from_smiles(smiles)

        # Verify the ligand has the correct SMILES string
        assert ligand.smiles_string == smiles

        # Verify that the file field is None
        assert ligand.file is None

        # Verify that properties are also None initially
        assert ligand.properties is None


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
    ligands = chemistry.Ligand.from_csv(
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
    ligands = chemistry.Ligand.from_csv(
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
    ligands = chemistry.Ligand.from_csv(file=csv_path, smiles_column="SMILES")

    # Should still have 4 valid ligands
    assert len(ligands) == 4

    # All ligands should have None properties
    for ligand in ligands:
        assert ligand.properties is None

    # Test with invalid SMILES column name
    with pytest.raises(
        ValueError, match="SMILES column 'InvalidColumn' not found in CSV file"
    ):
        chemistry.Ligand.from_csv(file=csv_path, smiles_column="InvalidColumn")


def test_download_protein():
    pdb_id = "1EBY"

    pdb_file = chemistry.download_protein(pdb_id)

    assert os.path.exists(pdb_file), "The downloaded PDB file does not exist."

    os.remove(pdb_file)


@pytest.mark.parametrize("ligand", ligands)
def test_filter_sdf_by_smiles(ligand):
    # Define the whitelist of canonical SMILES to filter by.
    # Replace these placeholder values with actual canonical SMILES strings.

    if ligand["n_ligands"] == 1:
        return

    smiles = chemistry.sdf_to_smiles(ligand["file"])

    smiles = smiles[:5]

    output_sdf = "temp.sdf"

    # Call the filtering function.
    chemistry.filter_sdf_by_smiles(
        input_sdf_file=ligand["file"],
        output_sdf_file=output_sdf,
        keep_only_smiles=smiles,
    )

    # Read back the filtered SDF file.
    actual_smiles = chemistry.sdf_to_smiles(output_sdf)

    actual_smiles = [chemistry.canonicalize_smiles(smi) for smi in actual_smiles]

    smiles = [chemistry.canonicalize_smiles(smi) for smi in smiles]

    assert set(actual_smiles) == set(smiles)

    os.remove(output_sdf)
