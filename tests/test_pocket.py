"""
Tests for the Pocket class.
"""

import numpy as np

from deeporigin.drug_discovery import Pocket, Protein


def test_pocket_from_ligand():
    protein = Protein.from_name("1EBY")

    ligand = protein.extract_ligand()

    pocket = Pocket.from_ligand(ligand)

    assert pocket.file_path is not None, "Pocket file path should not be None"


def test_pocket_get_center():
    """Test getting the center of a pocket."""
    protein = Protein.from_name("1EBY")
    ligand = protein.extract_ligand()
    pocket = Pocket.from_ligand(ligand)

    center = pocket.get_center()
    assert center is not None, "Center should be calculated"

    assert pocket.get_center().shape == (3,), "Pocket center shape is wrong"


def test_pocket_update_coordinates():
    """Test updating pocket coordinates."""
    protein = Protein.from_name("1EBY")
    ligand = protein.extract_ligand()
    pocket = Pocket.from_ligand(ligand)

    original_coords = pocket.coordinates.copy()
    new_coords = original_coords + 1.0  # Shift all coordinates by 1

    pocket.update_coordinates(new_coords)

    assert np.array_equal(pocket.coordinates, new_coords), (
        "Coordinates should be updated"
    )
