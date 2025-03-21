# example data
from importlib.resources import path

from deeporigin.drug_discovery import chemistry
from deeporigin.drug_discovery.complex import Complex

with path("deeporigin.data.brd", "brd.pdb") as file_path:
    EXAMPLE_DATA_DIR = file_path.parent
