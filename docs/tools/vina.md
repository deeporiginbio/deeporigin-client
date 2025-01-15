# AutoDock Vina

[AutoDock Vina](https://vina.scripps.edu/) is a molecular docking tool widely used in computational drug discovery. It predicts the binding mode of a small molecule (ligand) to a protein (receptor) by modeling their interactions. Vina is known for its high performance, accuracy, and user-friendly interface. It calculates the binding affinity and provides docking poses of the ligand in the receptor’s active site.

Vina is commonly used to screen potential drug candidates, study protein-ligand interactions, and explore binding mechanisms. It employs a scoring function to evaluate the strength of binding and an efficient optimization algorithm to search for the best docking configuration


## File Inputs 



### 1. Receptor File

- Format: [PDBQT](https://userguide.mdanalysis.org/2.6.0/formats/reference/pdbqt.html) (Protein Data Bank, with charges and torsions).
- This file represents the target (e.g., a protein or DNA structure).

Upload this file to the Data Hub. Typically, you would upload this file to a database dedicated to runs with this tool, in a column named `receptor`.

### 2. Ligand File

- Format: [PDBQT](https://userguide.mdanalysis.org/2.6.0/formats/reference/pdbqt.html).
- Represents the small molecule (e.g., a drug or compound).
- Prepared by calculating torsions and assigning charges using AutoDock Tools or MGLTools.
- Other input formats (e.g., PDB or MOL2) must be converted to PDBQT using preparation tools before use in AutoDock Vina.

Upload this file to the Data Hub. Typically, you would upload this file to a database dedicated to runs with this tool, in a column named `ligand`.

## Parameters

This section describes the parameters for a tool run, that are passed in the `start_run` function.


### Search Space

Defines the area of the receptor where docking will occur. This is critical for focusing on the active site or binding pocket.

#### `center_x, center_y, center_z`

- Coordinates of the center of the search box, specified in Ångstroms.
- Should be based on the binding site of the receptor (obtained from experimental data or visual inspection).

#### `size_x, size_y, size_z`

- Dimensions of the search box along each axis (in Ångstroms).
- Determines the search region’s size. A larger box covers more area but increases computation time.
- For most cases, sizes between 20-30 Ångstroms per side are typical for a flexible ligand.

### Docking Parameters

#### `energy_range`

- The energy difference (in kcal/mol) between the best pose and the worst acceptable pose.
- Smaller values prioritize only low-energy poses; larger values allow more diverse poses.

#### `exhaustiveness`

- Determines the thoroughness of the search.
- Higher values increase the number of sampling attempts, improving accuracy but requiring more computational time.
- Default is 8; lower values (e.g., 1-4) are faster but less exhaustive.

#### `num_modes`

- The maximum number of docking poses to generate.
- Vina will output up to this many unique poses for analysis.

## Running Vina on Deep Origin

To run AutoDock Vina on Deep Origin, follow these steps:


### 1. Create a database to store input and output files

### 2. Start a tool run on Deep Origin

For this tool run, we will use the following parameters:

```python
search_space = {
    "center_x": 15.190,
    "center_y": 53.903,
    "center_z": 16.917,
    "size_x": 1.1,
    "size_y": 1.1,
    "size_z": 1.1,
}

docking = {
    "energy_range": 0.3,
    "exhaustiveness": 1,
    "num_modes": 9,
}
```

To start a tool run, use:

```python
from deeporigin.tools import run

job_id = run.autodock_vina(
    database_id="<your-db-name>",
    row_id="<row-name>",
    search_space=search_space,
    docking=docking,
    output_column_name="<output-column-name>",
    receptor_column_name="<receptor-column-name>",
    ligand_column_name="<ligand-column-name>",
)
```

`run.autodock_vina` returns the ID of the tool run, that can be used to monitor the status of the run and terminate it if needed. 

`run.autodock_vina` prints a message that looks like:

```bash
🧬 Job started with ID: 9f7a3741-e392-45fb-a349-804b7fca07d7
```



## Outputs

### Docked Poses

- Format: PDBQT.
- Contains the ligand poses docked in the receptor.
- Each pose includes the coordinates, torsional flexibility, and orientation of the ligand.
