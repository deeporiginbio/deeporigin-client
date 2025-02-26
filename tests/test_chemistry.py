"""tests functions in the chemistry module"""

import os

import pytest
from deeporigin import chemistry

ligand_files = ["6xue-paper-ligands-docked.sdf", "ligands-brd-all.sdf", "brd-7.sdf"]
n_ligands = [44, 8, 1]


ligands = [
    {
        "file": os.path.join(os.path.dirname(__file__), "fixtures", file),
        "n_ligands": num,
    }
    for file, num in zip(ligand_files, n_ligands)
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
    )

    # Check that at least one output file was created
    sdf_files = list(output_dir.glob("*.sdf"))
    assert len(sdf_files) > 0, "No SDF files were created by the splitting function."

    assert len(sdf_file_names) == ligand["n_ligands"]

    # Optionally, you can validate each file (e.g., reading them back with RDKit)
    from rdkit import Chem

    for sdf_file in sdf_files:
        mols = list(Chem.SDMolSupplier(str(sdf_file), removeHs=False))
        assert len(mols) == 1, (
            f"Output file {sdf_file} does not contain exactly 1 molecule."
        )
