# PDBFixer

[OpenMM](https://github.com/openmm)'s [PDBFixer](https://github.com/openmm/pdbfixer) is a Python library for preparing molecular structures for simulations, particularly when working with files in the [Protein Data Bank](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html) (PDB) format. It can address common issues such as missing heavy atoms, hydrogen atoms, or residues, and can add missing chains or fix terminal capping. It is often used in the preprocessing phase of molecular dynamics simulations.

## File Inputs:


### PDB File

A PDB file containing the structure of the molecule or system to be fixed.


## File Outputs:

### Fixed PDB File

The corrected and prepared molecular structure, which can be saved as a PDB file.

## Parameters


### `addMissingResidues` 

| Type | Default |
| --- | --- |
| Bool | true |

Whether to add missing residues based on the structure and sequence.


### `addMissingAtoms`

| Type | Default |
| --- | --- |
| Bool | true |

### `addMissingHydrogens`

| Type | Default |
| --- | --- |
| Bool | true |

Whether to add hydrogen atoms to the structure, with optional control over pH for titration states.

### `pH` 

| Type | Default |
| --- | --- |
| Float | true |

### `removeHeterogens`

| Type | Default |
| --- | --- |
| Bool | true |

Whether to remove small molecules or ions from the structure, with options to retain water molecules.

### `removeWater`

| Type | Default |
| --- | --- |
| Bool | true |


