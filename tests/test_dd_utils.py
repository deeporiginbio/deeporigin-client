import os

from deeporigin.drug_discovery import BRD_DATA_DIR
from deeporigin.drug_discovery.external_tools.utils import count_atoms_in_pdb_file


def test_count_atoms_in_pdb_file():
    pdb_file_path = os.path.join(BRD_DATA_DIR, "brd.pdb")
    assert count_atoms_in_pdb_file(pdb_file_path) == 2503
