"""Shared test fixtures and setup code for drug discovery tests"""

import glob
import os

from deeporigin.drug_discovery import BRD_DATA_DIR

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


single_ligand_files = glob.glob(os.path.join(BRD_DATA_DIR, "*.sdf"))
single_ligand_files.sort()

single_ligand_hashes = [
    "2b689066dc53d6b81f579eacbaf08204a662528b2e0f482b3e12bb3bb0cd40a8",
    "45c64ae0e47040937e117edc6a37e9559e46a2c0017f697b2ae01f722e6b7638",
    "d9738bc8979d5e2cc510ff053d0f5de57013553693b288697b6fc893ed7f3a0d",
    "90f971e9523808dbfcab511157be359ac0143a6d483951aa9e5ec72ef87aa5b8",
    "e7196a8038b24d1fd8a6c4a4dc195ced555b9229b46ab056eef50c94bd96e6b5",
    "65cf3d57b4e4ddcd9258d7ca00b96ceb9626d4dc46d9c2a6063d06a27524ca3b",
    "94e35940aee80601d64cbc3147d914bcda33ed6d15c962965799c0e1b3538440",
    "e6a872953ee3dea7b00c9952ec9d6c646da77b993d48ddaf4e789c3381040192",
]

bad_ligands = [
    {
        "file": "missing.sdf",
        "smiles_string": None,
    },
    {
        "file": None,
        "smiles_string": None,
    },
]
