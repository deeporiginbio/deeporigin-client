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
    "95b16b7c028d5914ba0aaa8c1a123765c453f39199ccb0f85c728a3497b7e5af",
    "39c73ef3d7f99dc4d3b3ce505e2694252cffa8d55f8717613bd0d56391f9e748",
    "12ec11311b8cd62100ef6c64ab22a84da875902e28ced762d47fdd7804cc2ad4",
    "447df636990e04bc4165e266a6b94d9b9843cb779401d48ad8a27f414f38f6f7",
    "2da6dc605b748bd8b16eb833be73151b98cd0c89a97db990c26094eaf03529ee",
    "defcf9db5c113704a7f372d3a5c342c06c5641335cea4e355e65a53f2b1c3118",
    "4897c616dd1dfb06f08170b7f5ca67ce7f6779c2fd03970315c4387e80274225",
    "bac7b4d01c1a7ab102b1c9955a1839730a5099b08eba93807e12f6ab22adfb67",
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
