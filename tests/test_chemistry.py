"""tests functions in the chemistry module"""

import os

from deeporigin import chemistry

SDF_FILE = os.path.join(os.path.dirname(__file__), "fixtures", "ligands-brd-all.sdf")


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
