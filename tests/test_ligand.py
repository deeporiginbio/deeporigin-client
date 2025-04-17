import os

import pytest

from deeporigin.drug_discovery import chemistry
from deeporigin.drug_discovery.structures import Ligand
from deeporigin.exceptions import DeepOriginException

# Import shared test fixtures
from tests.utils_ligands import bad_ligands, ligands

base_path = os.path.join(os.path.dirname(__file__), "fixtures")


def test_ligands_from_sdf_file():
    """test that we can make many ligands from a single SDF file with many molecules"""

    mols = chemistry.read_molecules_in_sdf_file(
        os.path.join(base_path, "ligands-brd-all.sdf")
    )

    for mol in mols:
        Ligand.from_smiles(mol["smiles"], properties=mol["properties"])


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
    with pytest.raises(DeepOriginException) as exc_info:
        Ligand.from_smiles(smiles="InvalidSMILES")
    assert str(exc_info.value) == "Invalid SMILES string"


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
            "identifier": "oxotremorine",  # Muscarinic acetylcholine receptor agonist
            "name": "Oxotremorine",
            "expected_atoms": 15,  # Heavy atoms only
        },
        {
            "identifier": "serotonin",  # 5-hydroxytryptamine (5-HT)
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
    result = Ligand.from_sdf(ligand["file"])

    if n_ligands > 1:
        assert isinstance(result, list)
        assert len(result) == n_ligands
        for lig in result:
            assert isinstance(lig, Ligand)
            assert lig.mol is not None
            assert lig.mol.m.GetNumAtoms() > 0
            assert lig.file_path is None  # Multi-ligand case should have no file_path
    else:
        assert isinstance(result, Ligand)
        assert result.mol is not None
        assert result.mol.m.GetNumAtoms() > 0
        assert (
            result.file_path == ligand["file"]
        )  # Single ligand case should have file_path


def test_ligand_from_csv():
    """Test that we can create Ligands from a CSV file using the from_csv classmethod"""
    from pathlib import Path

    # Get the path to the test CSV file
    csv_path = Path(base_path) / "data.csv"

    # Create ligands from the CSV file
    ligands = Ligand.from_csv(str(csv_path), smiles_column="SMILES")

    # Verify we got the expected number of ligands
    assert len(ligands) == 30  # Total number of valid SMILES in the file

    # Check a few properties of the first ligand
    first_ligand = ligands[0]
    assert isinstance(first_ligand, Ligand)
    assert first_ligand.mol is not None
    assert first_ligand.mol.m.GetNumAtoms() > 0
    assert first_ligand.file_path is None

    # Verify properties were correctly loaded
    assert "score" in first_ligand.properties
    assert "binding_energy" in first_ligand.properties
    assert "pose_score" in first_ligand.properties

    # Test with invalid SMILES column
    with pytest.raises(ValueError, match="Column 'invalid' not found in CSV file"):
        Ligand.from_csv(str(csv_path), smiles_column="invalid")

    # Test with non-existent file
    with pytest.raises(FileNotFoundError):
        Ligand.from_csv("nonexistent.csv")


@pytest.mark.parametrize(
    "test_case",
    [
        pytest.param(
            {
                "file": "brd-7.sdf",
                "expected_type": Ligand,
                "expected_name": "cmpd 7 (n-propyl)",
                "has_atoms": True,
                "is_multi": False,
                "should_have_file_path": True,
            },
            id="single_molecule_sdf",
        ),
        pytest.param(
            {
                "file": "ligands-brd-all.sdf",
                "expected_type": list,
                "expected_name": None,  # Multiple molecules, names vary
                "has_atoms": True,
                "is_multi": True,
                "should_have_file_path": False,
            },
            id="multi_molecule_sdf",
        ),
        pytest.param(
            {
                "file": "nonexistent.sdf",
                "expected_error": FileNotFoundError,
            },
            id="invalid_file_path",
        ),
        pytest.param(
            {
                "file": "42-ligands.sdf",
                "expected_type": list,
                "expected_name": None,
                "has_atoms": True,
                "is_multi": True,
                "should_have_file_path": False,
            },
            id="large_multi_molecule_sdf",
        ),
        pytest.param(
            {
                "file": "brd-7.sdf",
                "expected_type": Ligand,
                "has_atoms": True,
                "sanitize": False,
                "should_have_file_path": True,
            },
            id="sanitize_false",
        ),
        pytest.param(
            {
                "file": "brd-7.sdf",
                "expected_type": Ligand,
                "has_atoms": True,
                "removeHs": True,
                "check_no_hydrogens": True,
                "should_have_file_path": True,
            },
            id="remove_hydrogens",
        ),
    ],
)
def test_ligand_from_sdf(test_case):
    """Test that we can create Ligand instances from SDF files using the from_sdf classmethod"""
    from pathlib import Path

    # Get the full path to the test file
    file_path = Path(base_path) / test_case["file"]

    # Handle expected error cases
    if "expected_error" in test_case:
        with pytest.raises(test_case["expected_error"]) as exc_info:
            if "error_message" in test_case:
                Ligand.from_sdf(str(file_path))
                assert test_case["error_message"] in str(exc_info.value)
            else:
                Ligand.from_sdf(str(file_path))
        return

    # Handle normal cases
    kwargs = {}
    if "sanitize" in test_case:
        kwargs["sanitize"] = test_case["sanitize"]
    if "removeHs" in test_case:
        kwargs["removeHs"] = test_case["removeHs"]

    result = Ligand.from_sdf(str(file_path), **kwargs)

    # Verify the result type
    assert isinstance(result, test_case["expected_type"])

    # Handle single molecule case
    if not test_case.get("is_multi", False):
        ligand = result
        assert ligand.mol is not None
        assert ligand.mol.m.GetNumAtoms() > 0

        # Check file_path based on whether it should have one
        if test_case.get("should_have_file_path", False):
            assert ligand.file_path == str(file_path)
        else:
            assert ligand.file_path is None

        if "expected_name" in test_case:
            assert ligand.name == test_case["expected_name"]

        if test_case.get("check_no_hydrogens", False):
            assert all(atom.GetAtomicNum() != 1 for atom in ligand.mol.m.GetAtoms())

    # Handle multi-molecule case
    else:
        assert len(result) > 1
        for ligand in result:
            assert isinstance(ligand, Ligand)
            assert ligand.mol is not None
            assert ligand.mol.m.GetNumAtoms() > 0
            # In multi-molecule case, file_path should always be None
            assert ligand.file_path is None
