import os

import pytest

from deeporigin.drug_discovery import BRD_DATA_DIR, Protein, Pocket
from deeporigin.exceptions import DeepOriginException


def test_from_residue_num():
    """Test creating a pocket from a residue number"""

    # Read in a protein
    pdb_path = os.path.join(BRD_DATA_DIR, "brd.pdb")
    protein = Protein.from_file(pdb_path)

    # Create custom pocket
    custom_pocket = Pocket.from_residue_number(protein, residue_number=77, cutoff=5)

    assert len(custom_pocket) == 1, "Incorrect number of pockets generated."
