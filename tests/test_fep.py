import os

import pytest
from deeporigin.tools import fep

ligand_files = ["6xue-paper-ligands-docked.sdf", "ligands-brd-all.sdf", "brd-7.sdf"]
n_ligands = [44, 8, 1]


ligands = [
    {
        "file": os.path.join(os.path.dirname(__file__), "fixtures", file),
        "n_ligands": num,
    }
    for file, num in zip(ligand_files, n_ligands)
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
