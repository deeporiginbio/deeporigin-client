import os
from pathlib import Path

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

    # Verify we can still compute the hash
    assert complex_obj._hash is not None


def test_hash_changes_with_ligands():
    """Test that the Complex hash changes when ligands are modified"""
    # Create initial complex with protein only
    pdb_path = os.path.join(EXAMPLE_DATA_DIR, "brd.pdb")
    protein = Protein.from_file(pdb_path)
    complex_obj = Complex(protein=protein, ligands=[])

    # Store initial hash
    initial_hash = complex_obj._hash
    previous_hash = initial_hash

    # Add ligands one by one and verify hash changes each time
    sdf_files = sorted([f for f in os.listdir(EXAMPLE_DATA_DIR) if f.endswith(".sdf")])
    for sdf_file in sdf_files[:3]:  # Test with first 3 ligands
        sdf_path = os.path.join(EXAMPLE_DATA_DIR, sdf_file)
        new_ligand = Ligand.from_sdf(sdf_path)
        if isinstance(new_ligand, list):
            new_ligands = new_ligand
        else:
            new_ligands = [new_ligand]

        # Add new ligands to existing ones
        complex_obj.ligands = complex_obj.ligands + new_ligands

        # Verify hash has changed
        current_hash = complex_obj._hash
        assert current_hash != previous_hash, (
            f"Hash didn't change after adding {sdf_file}"
        )
        previous_hash = current_hash

    # Verify final hash is different from initial
    assert complex_obj._hash != initial_hash, (
        "Final hash should be different from initial hash"
    )
