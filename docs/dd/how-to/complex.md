This document describes how to create a `Complex` object, that can be used to run Docking, ABFE and RBFE. 


## Creating a Complex

### From a directory 

```python
# here, we're using the example data directory
from deeporigin.drug_discovery import Complex, BRD_DATA_DIR

sim = Complex.from_dir(BRD_DATA_DIR)
```

The directory should contain:

- exactly one PDB file for the protein
- one or more SDF files for the ligands. Each SDF file can contain one or more molecules.

### From `Protein` and `Ligand` objects

A `Complex` object can be also be constructed using `Protein` and `Ligand` objects. 

```python
from deeporigin.drug_discovery import Complex, BRD_DATA_DIR, Protein, Ligand

protein = Protein.from_file(BRD_DATA_DIR / "brd.pdb")
ligand = Ligand.from_sdf(BRD_DATA_DIR / "brd-2.sdf")

sim = Complex(protein=protein, ligands=[ligand])
```

### From `LigandSet` objects


A `Complex` object can also be constructed using `Protein` and `LigandSet` objects. 

```python
from deeporigin.drug_discovery import Complex, BRD_DATA_DIR, Protein, LigandSet

protein = Protein.from_file(BRD_DATA_DIR / "brd.pdb")
ligands = LigandSet.from_dir(BRD_DATA_DIR)

sim = Complex(protein=protein, ligands=ligands)
```

## Modifying a Complex

You can modify a Complex object by adding or replacing ligands. 

```python
from deeporigin.drug_discovery import Complex, Protein, Ligand, BRD_DATA_DIR
protein = Protein.from_file(BRD_DATA_DIR / "brd.pdb")
ligand = Ligand.from_sdf(BRD_DATA_DIR / "brd-2.sdf")
ligand3 = Ligand.from_sdf(BRD_DATA_DIR / "brd-3.sdf")

# Create a complex with just the protein
sim = Complex(protein=protein)

# Add a single ligand
sim.ligands = ligand3

# Add another ligand
ligand4 = Ligand.from_sdf(BRD_DATA_DIR / "brd-4.sdf")
sim.ligands = sim.ligands + ligand4

# Replace all ligands
ligand5 = Ligand.from_sdf(BRD_DATA_DIR / "brd-5.sdf")
sim.ligands = ligand5  # Replace with a single ligand
```

??? tip "Constructing Ligands"
    To see how to construct Ligands, see [this](./ligands.md)

??? tip "Constructing the Protein"
    To see how to construct the Protein, see [this](./proteins.md)

## Preparing a complex

To prepare a protein-ligand complex for simulation or further analysis, use the `Complex.prepare()` method. This method runs system preparation on a given ligand in the context of the complex's protein, handling tasks such as protonation, water retention, and box padding.

### Usage

```{.python notest}
sim.prepare(ligand)
```

Typically, you would call this method using a ligand in the complex:

```{.python notest}
sim.prepare(sim.ligands[0])
```

### Arguments

- `ligand` (`Ligand`): The ligand to prepare. This should be an instance of the `Ligand` class associated with your complex.
- `padding` (`float`, optional): Padding (in angstroms) to add around the system. Default is `1.0`.
- `keep_waters` (`bool`, optional): Whether to keep water molecules in the prepared system. Default is `False`.
- `is_lig_protonated` (`bool`, optional): Whether the ligand is already protonated. Default is `True`.
- `is_protein_protonated` (`bool`, optional): Whether the protein is already protonated. Default is `True`.

You can override any of these defaults as needed:

```{.python notest}
sim.prepare(
    ligand,
    padding=2.0,
    keep_waters=True,
    is_lig_protonated=False,
    is_protein_protonated=False,
)
```

### What happens during preparation?

- The method calls the system preparation backend, which processes the protein and ligand files.
- The prepared complex is saved to a cache for efficiency.
- The resulting structure can be visualized or used for downstream workflows.

After preparation, the prepared complex structure will be shown automatically.

