This document describes how to create a `Complex` object, that can be used to run Docking, ABFE and RBFE. 

## From a directory 

```python
from deeporigin.drug_discovery import Complex

# here, we're using the example data directory
sim = Complex.from_dir(EXAMPLE_DATA_DIR)
```

The directory should contain:

- exactly one PDB file 
- `N` SDF files, with one molecule per SDF file

## From `Protein` and `Ligand` objects

A `Complex` object can be also be constructed using `Protein` and `Ligand` objects. 

```python
from deeporigin.drug_discovery import Complex

sim = Complex(protein=protein, ligands=ligands)
```

??? tip "Constructing Ligands"
    To see how to construct Ligands, see [this](./ligands.md)

??? tip "Constructing the Protein"
    To see how to construct the Protein, see [this](./proteins.md)