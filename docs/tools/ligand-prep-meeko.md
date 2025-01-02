

[Meeko](https://meeko.readthedocs.io/en/release/lig_prep_basic.html) is a Python library that facilitates the preparation of small-molecule ligands for docking simulations, particularly with AutoDock tools. It automates tasks such as protonation, torsion bond assignment, and preparation of input files in the appropriate format.

By ensuring consistent and reliable ligand preparation, Meeko plays a vital role in generating high-quality input files for docking, such as with [AutoDock Vina](vina.md).

## File Inputs


### 1. Ligand Structure File

A molecule file in [SDF format](https://pmc.ncbi.nlm.nih.gov/articles/PMC3163917/). This file provides the 3D atomic coordinates and bond information for the ligand.


## Output file

### 1. PDBQT File

The prepared ligand in the [PDBQT format](https://userguide.mdanalysis.org/2.6.0/formats/reference/pdbqt.html), which includes atomic coordinates, torsions, charges, and connectivity data.


## Parameters

No parameters can be passed to this tool.

