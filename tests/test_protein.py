import os

import pytest

from deeporigin.drug_discovery import Protein
from deeporigin.exceptions import DeepOriginException


def test_from_pdb_id():
    conotoxin = Protein.from_pdb_id("2JUQ")

    os.remove(conotoxin.file_path)

    conotoxin = Protein.from_pdb_id("2JUQ")


def test_from_pdb_id_with_invalid_id():
    with pytest.raises(DeepOriginException, match=r".*Failed to create Protein.*"):
        Protein.from_pdb_id("foobar")
