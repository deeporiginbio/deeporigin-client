# example data
from importlib.resources import path

from deeporigin.drug_discovery.structures import Protein, Ligand, Pocket

__all__ = ["chemistry", "Complex", "Protein", "Ligand", "Pocket"]


# Lazy imports for better performance
def __getattr__(name):
    if name == "Complex":
        from deeporigin.drug_discovery.complex import Complex

        return Complex
    elif name == "chemistry":
        from deeporigin.drug_discovery import chemistry

        return chemistry

    raise AttributeError(f"module {__name__} has no attribute {name}")


with path("deeporigin.data.brd", "brd.pdb") as file_path:
    EXAMPLE_DATA_DIR = file_path.parent
