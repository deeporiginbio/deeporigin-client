"""tests functions in the chemistry module"""

import os

import pytest
from deeporigin import chemistry
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
