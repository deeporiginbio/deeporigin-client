This document describes how to create a `Complex` object, that can be used to run Docking, ABFE and RBFE. 

## From a directory 

```python
from deeporigin.drug_discovery import Complex

# here, we're using the example data directory
sim = Complex.from_dir(EXAMPLE_DATA_DIR)
```

The directory should contain:

- exactly one PDB file for the protein
- one or more SDF files for the ligands. Each SDF file can contain one or more molecules.

## From `Protein` and `Ligand` objects

A `Complex` object can be also be constructed using `Protein` and `Ligand` objects. 

```python
from deeporigin.drug_discovery import Complex

sim = Complex(protein=protein, ligands=ligands)
```

## Modifying a Complex

You can modify a Complex object by adding or replacing ligands. The Complex's hash will automatically update when ligands are modified.

```python
# Create a complex with just the protein
complex = Complex(protein=protein, ligands=[])

# Add a single ligand
new_ligand = Ligand.from_sdf("ligand.sdf")
complex.ligands = complex.ligands + [new_ligand]

# Add multiple ligands
more_ligands = Ligand.from_sdf("multiple_ligands.sdf")  # Returns a list if file contains multiple molecules
complex.ligands = complex.ligands + more_ligands

# Replace all ligands
complex.ligands = [new_ligand]  # Replace with a single ligand
```

??? tip "Constructing Ligands"
    To see how to construct Ligands, see [this](./ligands.md)

??? tip "Constructing the Protein"
    To see how to construct the Protein, see [this](./proteins.md)