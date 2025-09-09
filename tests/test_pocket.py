"""
Tests for the Pocket class.
"""

from deeporigin.drug_discovery import Pocket, Protein


def test_pocket_from_ligand():
    protein = Protein.from_name("1EBY")

    ligand = protein.extract_ligand()

    pocket = Pocket.from_ligand(ligand)

    assert pocket.file_path is not None, "Pocket file path should not be None"
