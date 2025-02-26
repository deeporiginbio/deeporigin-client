"""tests functions in the chemistry module"""

import os

from deeporigin import chemistry

SDF_FILE = os.path.join(os.path.dirname(__file__), "fixtures", "ligands-brd-all.sdf")


def test_sdf_to_smiles():
    smiles_strings = chemistry.sdf_to_smiles(SDF_FILE)

    assert set(smiles_strings) == {
        "C=CCCn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
        "CN(C)C(=O)c1cccc(-c2cn(C)c(=O)c3[nH]ccc23)c1",
        "C=CCn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
        "COCCn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
        "CCn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
        "CCCCn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
        "CCCn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
        "C/C=C/Cn1cc(-c2cccc(C(=O)N(C)C)c2)c2cc[nH]c2c1=O",
    }


def test_split_sdf_file(tmp_path):
    """
    Test that split_sdf_using_names correctly splits the ligands SDF file
    into separate SDF files, and that the output is cleaned up (by pytest)
    after the test completes.
    """

    # Create an output directory within the pytest temp directory
    output_dir = tmp_path / "split_ligands"
    output_dir.mkdir(exist_ok=True)

    # Call the function to be tested
    chemistry.split_sdf_file(
        input_sdf_path=str(SDF_FILE),
        output_prefix="testLig",
        output_dir=str(output_dir),
    )

    # Check that at least one output file was created
    sdf_files = list(output_dir.glob("*.sdf"))
    assert len(sdf_files) > 0, "No SDF files were created by the splitting function."

    # Optionally, you can validate each file (e.g., reading them back with RDKit)
    from rdkit import Chem

    for sdf_file in sdf_files:
        mols = list(Chem.SDMolSupplier(str(sdf_file), removeHs=False))
        assert len(mols) == 1, (
            f"Output file {sdf_file} does not contain exactly 1 molecule."
        )
