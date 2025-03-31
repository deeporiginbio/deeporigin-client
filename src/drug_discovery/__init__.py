# example data
from importlib.resources import path

# eager import
from deeporigin.drug_discovery import chemistry

__all__ = ["chemistry", "Complex", "Protein", "Ligand"]


def __getattr__(name):
    if name == "Complex":
        # lazy import
        from deeporigin.drug_discovery.complex import Complex

        return Complex
    elif name == "Protein":
        # lazy import
        from deeporigin.drug_discovery.chemistry import Protein

        return Protein
    elif name == "Ligand":
        # lazy import
        from deeporigin.drug_discovery.chemistry import Ligand

        return Ligand
    raise AttributeError(f"module {__name__} has no attribute {name}")


with path("deeporigin.data.brd", "brd.pdb") as file_path:
    EXAMPLE_DATA_DIR = file_path.parent
