import os

import pytest

from deeporigin.drug_discovery.structures import Ligand
from deeporigin.exceptions import DeepOriginException

# Import shared test fixtures
from tests.utils_ligands import bad_ligands, ligands

base_path = os.path.join(os.path.dirname(__file__), "fixtures")


def test_ligand_from_smiles():
    """Test that we can create a Ligand from a SMILES string using the from_smiles classmethod"""
    from rdkit import Chem

    # Test cases with different types of molecules
    test_cases = [
        {
            "smiles": "C",  # Methane
            "name": "Methane",
            "expected_atoms": 1,
        },
        {
            "smiles": "CC",  # Ethane
            "name": "Ethane",
            "expected_atoms": 2,
        },
        {
            "smiles": "CCO",  # Ethanol
            "name": "Ethanol",
            "expected_atoms": 3,
        },
        {
            "smiles": "c1ccccc1",  # Benzene (aromatic notation)
            "name": "Benzene",
            "expected_atoms": 6,
            "equivalent_smiles": ["C1=CC=CC=C1"],  # Benzene (Kekule notation)
        },
        {
            "smiles": "C1=CC=CC=C1",  # Benzene (Kekule notation)
            "name": "Benzene_alt",
            "expected_atoms": 6,
            "equivalent_smiles": ["c1ccccc1"],  # Benzene (aromatic notation)
        },
    ]

    for case in test_cases:
        # Create a ligand using the from_smiles method
        ligand = Ligand.from_smiles(
            smiles=case["smiles"],
            name=case["name"],
            save_to_file=False,
        )

        # Verify the ligand has either the exact SMILES string or an equivalent one
        if "equivalent_smiles" in case:
            assert ligand.smiles in [case["smiles"]] + case["equivalent_smiles"], (
                f"SMILES {ligand.smiles} not equivalent to {case['smiles']} or any of {case['equivalent_smiles']}"
            )
        else:
            assert ligand.smiles == case["smiles"]

        # Verify the name was set correctly
        assert ligand.name == case["name"]

        # Verify that the file field is None
        assert ligand.file_path is None

        # Verify that the molecule was properly initialized
        assert ligand.mol is not None
        assert ligand.mol.m.GetNumAtoms() == case["expected_atoms"]

        # Verify that the molecule represents the same chemical structure
        input_mol = Chem.MolFromSmiles(case["smiles"])
        assert Chem.MolToSmiles(input_mol) == Chem.MolToSmiles(ligand.mol.m)

    # Test with invalid SMILES
    with pytest.raises(DeepOriginException, match=r"Cannot create"):
        Ligand.from_smiles(smiles="InvalidSMILES")


# this is marked as such because the pubchem API has become very flaky
@pytest.mark.xfail(strict=False)
def test_ligand_from_identifier():
    """Test that we can create a Ligand from common biochemical identifiers using the from_identifier classmethod"""

    # Test cases with different biologically relevant molecules
    test_cases = [
        {
            "identifier": "ATP",  # Adenosine triphosphate
            "name": "ATP",
            "expected_atoms": 31,  # Heavy atoms only
        },
        {
            "identifier": "ADP",  # Adenosine diphosphate
            "name": "ADP",
            "expected_atoms": 27,  # Heavy atoms only
        },
        {
            "identifier": "Oxotremorine",  # Muscarinic acetylcholine receptor agonist
            "name": "Oxotremorine",
            "expected_atoms": 15,  # Heavy atoms only
        },
        {
            "identifier": "Serotonin",  # 5-hydroxytryptamine (5-HT)
            "name": "Serotonin",
            "expected_atoms": 13,  # Heavy atoms only (N,C,C,c,c,n,c,c,c,c,O,c,c)
        },
    ]

    for case in test_cases:
        # Create a ligand using the from_identifier method
        ligand = Ligand.from_identifier(
            identifier=case["identifier"],
            name=case["name"],
            save_to_file=False,
        )

        # Verify the name was set correctly
        assert ligand.name == case["name"]

        # Verify that the file field is None
        assert ligand.file_path is None

        # Verify that the molecule was properly initialized
        assert ligand.mol is not None
        assert ligand.mol.m.GetNumAtoms() == case["expected_atoms"]

        # Verify that the identifier was stored
        assert ligand.identifier == case["identifier"]

        # Verify that the molecule has valid 3D coordinates
        assert ligand.mol.m.GetNumConformers() > 0
        coords = ligand.mol.m.GetConformer().GetPositions()
        assert coords.shape[0] == case["expected_atoms"]

    # Test with invalid identifier
    invalid_id = "InvalidMolecule123"
    with pytest.raises(DeepOriginException) as exc_info:
        Ligand.from_identifier(identifier=invalid_id)
    assert str(exc_info.value).startswith(
        f"Could not resolve chemical identifier '{invalid_id}'"
    )


def test_ligand_from_rdkit_mol():
    """Test that we can create a Ligand from an RDKit Mol object using the from_rdkit_mol classmethod"""
    from rdkit import Chem

    # Create test RDKit molecules
    mols = [
        Chem.MolFromSmiles("C"),  # Methane
        Chem.MolFromSmiles("CC"),  # Ethane
        Chem.MolFromSmiles("CCO"),  # Ethanol
        Chem.MolFromSmiles("c1ccccc1"),  # Benzene
    ]

    for mol in mols:
        # Create a ligand using the from_rdkit_mol method
        ligand = Ligand.from_rdkit_mol(mol, name="TestLigand")

        # Verify the ligand has the correct SMILES string
        assert ligand.smiles == Chem.MolToSmiles(mol)

        # Verify the name was set correctly
        assert ligand.name == "TestLigand"

        # Verify that the file field is None
        assert ligand.file_path is None

        # Verify that the molecule was properly initialized
        assert ligand.mol is not None
        assert ligand.mol.m.GetNumAtoms() == mol.GetNumAtoms()


@pytest.mark.parametrize("ligand", bad_ligands)
def test_ligand_errors(ligand):
    with pytest.raises(Exception):  # noqa: B017
        Ligand(
            file_path=ligand["file"],
            smiles=ligand["smiles"],
        )


@pytest.mark.parametrize("ligand", ligands)
def test_ligand(ligand):
    """Test that we can create Ligand instances from various sources"""
    n_ligands = ligand["n_ligands"]

    if n_ligands >= 1:
        return

    result = Ligand.from_sdf(ligand["file"])

    assert isinstance(result, Ligand)
    assert result.mol is not None
    assert result.mol.m.GetNumAtoms() > 0
    assert (
        result.file_path == ligand["file"]
    )  # Single ligand case should have file_path


def test_ligand_from_sdf_multiple_raises():
    """Test that Ligand.from_sdf raises DeepOriginException for multi-molecule SDF files."""
    with pytest.raises(
        DeepOriginException, match="must contain exactly one molecule, but found 8"
    ):
        Ligand.from_sdf(os.path.join(base_path, "ligands-brd-all.sdf"))
