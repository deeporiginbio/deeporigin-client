"""tests functions in the chemistry module"""

import os

import pytest

from deeporigin.drug_discovery import chemistry

# Import shared test fixtures
from tests.utils_ligands import ligands

base_path = os.path.join(os.path.dirname(__file__), "fixtures")


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
