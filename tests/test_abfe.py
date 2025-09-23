"""tests for abfe"""

import pytest

from deeporigin.drug_discovery import BRD_DATA_DIR, Complex, Ligand, Protein
from deeporigin.exceptions import DeepOriginException


def test_abfe_charged_ligand():
    """test that abfe raises an error if a charged ligand is provided"""

    ligand = Ligand.from_smiles("C[N+]1=CCCC1")
    protein = Protein.from_file(BRD_DATA_DIR / "brd.pdb")
    sim = Complex(protein=protein, ligands=ligand)

    with pytest.raises(
        DeepOriginException,
        match="ABFE does not currently support charged ligands",
    ):
        sim.abfe.run()


def test_abfe_prepared_system():
    """test that abfe raises an error if a prepared system is not provided"""
    ligand = Ligand.from_smiles("CCO")
    protein = Protein.from_file(BRD_DATA_DIR / "brd.pdb")
    sim = Complex(protein=protein, ligands=ligand)
    with pytest.raises(DeepOriginException, match="Please prepare the system using"):
        sim.abfe.run()
