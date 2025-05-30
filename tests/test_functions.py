"""this module contains tests for functions. These are meant to be run against a live instance"""

import pytest

from tests.utils import config  # noqa: F401


def test_molprops(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")

    from deeporigin.drug_discovery import Ligand

    ligand = Ligand.from_identifier("serotonin")

    props = ligand.admet_properties()

    assert isinstance(props, dict), "Expected a dictionary"
    assert "logP" in props, "Expected logP to be in the properties"
    assert "logD" in props, "Expected logD to be in the properties"
    assert "logS" in props, "Expected logS to be in the properties"


def test_pocket_finder(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")
    from deeporigin.drug_discovery import Protein

    protein = Protein.from_pdb_id("1EBY")
    pockets = protein.find_pockets(pocket_count=1)

    assert len(pockets) == 1, "Incorrect number of pockets"


def test_docking(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")
    from deeporigin.drug_discovery import Ligand, Protein

    protein = Protein.from_pdb_id("1EBY")
    pockets = protein.find_pockets(pocket_count=1)
    pocket = pockets[0]

    ligand = Ligand.from_smiles("CN(C)C(=O)c1cccc(-c2cn(C)c(=O)c3[nH]ccc23)c1")

    poses_sdf = protein.dock(ligand=ligand, pocket=pocket)

    assert isinstance(poses_sdf, str), (
        "Expected a string to be returned by dock function"
    )


def test_sysprep(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")

    from deeporigin.drug_discovery import EXAMPLE_DATA_DIR, Complex
    from deeporigin.functions.sysprep import sysprep

    sim = Complex.from_dir(EXAMPLE_DATA_DIR)

    # this is chosen to be one where it takes >1 min
    sysprep(
        protein_path=sim.protein.file_path,
        ligand_path=sim.ligands[3].file_path,
        is_lig_protonated=True,
    )
