This document describes how to use the outputs of Docking as inputs to FEP tools.

## Assumptions

We assume that you have 

1. [:material-page-previous: created a `Complex` object](../tutorial/getting-started.md)
2. [:material-page-previous: run Docking](../tutorial/docking.md)
3. [:material-page-previous: Generated a SDF file with docked poses](../tutorial/docking.md#exporting-a-sdf-with-docked-poses)

## Split Docking outputs SDF file

Docking generates a single SDF file with all the ligands in it. FEP tools (ABFE and RBFE) require SDF files with individual molecules in their own SDF file. To get here, we use a utility function as follows:

```{.python notest}
from deeporigin.drug_discovery import chemistry

sdf_files = chemistry.split_sdf_file(
    input_sdf_path="filtered_ligands.sdf", 
)
```

`sdf_files` contains a list of paths to the generated SDF files. 

We can use this to generate a new list of Ligands and modify the `Complex`


```{.python notest}
from deeporigin.drug_discovery import Ligand

sim.ligands = [Ligand(file) for file in files]
```

This object can be used to run [:material-page-previous: ABFE](../tutorial/abfe.md) or [:material-page-previous:RBFE](../tutorial/rbfe.md). 