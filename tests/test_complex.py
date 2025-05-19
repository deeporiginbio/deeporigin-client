import os

from deeporigin.drug_discovery import EXAMPLE_DATA_DIR, Complex, Ligand, Protein


def test_from_dir():
    """Test creating a Complex from the example data directory"""
    # Create Complex from directory
    complex_obj = Complex.from_dir(EXAMPLE_DATA_DIR)

    # Verify the complex was created correctly
    assert isinstance(complex_obj, Complex)
    assert isinstance(complex_obj.protein, Protein)

    # Verify we have the expected number of ligands
    # The example directory contains multiple SDF files with multiple molecules each
    assert len(complex_obj.ligands) > 0

    # Verify all ligands are valid
    for ligand in complex_obj.ligands:
        assert isinstance(ligand, Ligand)
        assert ligand.smiles is not None
        assert ligand.name is not None

    # Verify the protein is valid
    assert complex_obj.protein.name == "brd"


def test_construct_complex():
    """Test constructing a Complex by providing protein and ligands directly"""
    # Create protein from PDB file
    pdb_path = os.path.join(EXAMPLE_DATA_DIR, "brd.pdb")
    protein = Protein.from_file(pdb_path)

    # Create ligands from SDF files
    ligands = []
    sdf_files = ["brd-2.sdf", "brd-3.sdf"]  # Testing with a subset of ligands
    for sdf_file in sdf_files:
        sdf_path = os.path.join(EXAMPLE_DATA_DIR, sdf_file)
        result = Ligand.from_sdf(sdf_path)
        if isinstance(result, list):
            ligands.extend(result)
        else:
            ligands.append(result)

    # Create Complex
    complex_obj = Complex(protein=protein, ligands=ligands)

    # Verify the complex was created correctly
    assert isinstance(complex_obj, Complex)
    assert isinstance(complex_obj.protein, Protein)
    assert complex_obj.protein.name == "brd"
    assert len(complex_obj.ligands) > 0

    # Verify all ligands are valid
    for ligand in complex_obj.ligands:
        assert isinstance(ligand, Ligand)
        assert ligand.smiles is not None
        assert ligand.name is not None
