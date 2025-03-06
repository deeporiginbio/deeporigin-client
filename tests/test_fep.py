import os

import pytest
from deeporigin.tools import fep

base_path = os.path.join(os.path.dirname(__file__), "fixtures")

ligands = [
    {
        "file": os.path.join(base_path, "42-ligands.sdf"),
        "n_ligands": 42,
        "name_by_property": "Compound",
    },
    {
        "file": os.path.join(base_path, "ligands-brd-all.sdf"),
        "n_ligands": 8,
        "name_by_property": "_Name",
    },
    {
        "file": os.path.join(base_path, "brd-7.sdf"),
        "n_ligands": 1,
        "name_by_property": "_Name",
    },
]


@pytest.mark.parametrize("ligand", ligands)
def test_ligand(
    ligand,
):
    n_ligands = ligand["n_ligands"]

    if n_ligands > 1:
        with pytest.raises(ValueError, match="Too many molecules."):
            ligand = fep.Ligand(ligand["file"])
    else:
        ligand = fep.Ligand(ligand["file"])
