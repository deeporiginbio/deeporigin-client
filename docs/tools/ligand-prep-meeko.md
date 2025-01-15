

[Meeko](https://meeko.readthedocs.io/en/release/lig_prep_basic.html) is a Python library that facilitates the preparation of small-molecule ligands for docking simulations, particularly with AutoDock tools. It automates tasks such as protonation, torsion bond assignment, and preparation of input files in the appropriate format.

By ensuring consistent and reliable ligand preparation, Meeko plays a vital role in generating high-quality input files for docking, such as with [AutoDock Vina](vina.md).

## File Inputs


### Ligand Structure File

A molecule file in [SDF format](https://pmc.ncbi.nlm.nih.gov/articles/PMC3163917/). This file provides the 3D atomic coordinates and bond information for the ligand.


## Parameters

No parameters can be passed to this tool.

## Running on Deep Origin

To run this tool on Deep Origin, follow these steps:

### 1. Create a database to store input and output files

Navigate to DataHub on Deep Origin, and create a database for this tool. Add columns for input and output files. Create a new row and upload your input files to the input column. 

### 2. Start a tool run on Deep Origin

For this tool run, we will use the following parameters:


To start a tool run, use:

```python
from deeporigin.tools import run

job_id = run.ligand_prep(
    database_id="<your-db-name>",
    row_id="<row-name>",
    output_column_name="<output-column-name>",
    ligand_column_name="<receptor-column-name>",
)
```

`run.ligand_prep` returns the ID of the tool run, that can be used to monitor the status of the run and terminate it if needed. 

`run.ligand_prep` prints a message that looks like:

```bash
🧬 Job started with ID: 9f7a3741-e392-45fb-a349-804b7fca07d7
```


## Outputs

### PDBQT File

The prepared ligand in the [PDBQT format](https://userguide.mdanalysis.org/2.6.0/formats/reference/pdbqt.html), which includes atomic coordinates, torsions, charges, and connectivity data.
