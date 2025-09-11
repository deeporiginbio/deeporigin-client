import os

import numpy as np

from deeporigin.drug_discovery import BRD_DATA_DIR, Pocket, Protein


def test_from_residue_num():
    """Test creating a pocket from a residue number"""

    # Read in a protein
    pdb_path = os.path.join(BRD_DATA_DIR, "brd.pdb")
    protein = Protein.from_file(pdb_path)

    # Create custom pocket
    custom_pocket = Pocket.from_residue_number(protein, residue_number=77, cutoff=5)

    assert isinstance(
        custom_pocket.get_center(), np.ndarray
    ) and custom_pocket.get_center().shape == (3,)
