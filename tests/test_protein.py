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


def test_sequence():
    protein = Protein.from_pdb_id("1EBY")
    sequences = protein.sequence
    assert len(sequences) == 2

    for seq in sequences:
        assert (
            seq
            == "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
        )


def test_find_missing_residues():
    protein = Protein.from_pdb_id("5QSP")
    missing = protein.find_missing_residues()
    # The expected output is based on the documentation example
    expected = {
        "A": [(511, 514), (547, 550), (679, 682), (841, 855)],
        "B": [(509, 516), (546, 551), (679, 684), (840, 854)],
    }
    assert missing == expected


def test_pdb_id():
    protein = Protein.from_pdb_id("1EBY")
    assert protein.pdb_id == "1EBY"
