import os

import pytest

from deeporigin.drug_discovery import BRD_DATA_DIR, Protein
from deeporigin.exceptions import DeepOriginException


def test_from_file():
    protein = Protein.from_file(BRD_DATA_DIR / "brd.pdb")

    assert (
        str(protein.sequence[0])
        == "STNPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDMGTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTE"
    )


def test_from_name():
    protein = Protein.from_name("conotoxin")
    assert protein.pdb_id == "1P1P"

    assert str(protein.sequence[0]) == "GCCGSYPNAACHPCSCKDR"


def test_from_pdb_id():
    conotoxin = Protein.from_pdb_id("2JUQ")

    os.remove(conotoxin.file_path)

    _ = Protein.from_pdb_id("2JUQ")


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


def test_extract_ligand():
    protein = Protein.from_pdb_id("1EBY")
    ligand = protein.extract_ligand()

    assert (
        ligand.smiles
        == "OC(N[C@H]1C2CCCCC2C[C@H]1O)[C@H](OCC1CCCCC1)[C@H](O)[C@@H](O)[C@@H](OCC1CCCCC1)[C@@H](O)N[C@H]1C2CCCCC2C[C@H]1O"
    )


def test_protein_base64():
    """Test that we can convert a Protein to base64 and back"""
    # Create a protein using from_pdb_id
    protein = Protein.from_pdb_id("1EBY")

    # Convert to base64
    b64 = protein.to_base64()

    # Convert back from base64
    new_protein = Protein.from_base64(b64)

    # Verify the structures have the same number of atoms
    assert len(new_protein.structure) == len(protein.structure)

    # Verify the structures have the same coordinates (within numerical precision)
    import numpy as np

    np.testing.assert_array_almost_equal(
        new_protein.structure.coord,
        protein.structure.coord,
        decimal=3,
    )
