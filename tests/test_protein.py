import os

from deeporigin.drug_discovery import Protein


def test_from_pdb_id():
    conotoxin = Protein.from_pdb_id("2JUQ")

    os.remove(conotoxin.file_path)

    conotoxin = Protein.from_pdb_id("2JUQ")
