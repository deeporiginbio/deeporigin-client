This document describes how to filter the outputs of Docking based on various properties. 

Here we assume that you have constructed a `Complex` object and successfully run [Docking](../tutorial/docking.md). 
Following convention, we assume that the `Complex` object is called `sim`.

## Fetch Docking results

First, we get the results of Docking in a pandas DataFrame using:

```{.python notest}
df = sim.docking.get_results()
```
## Fetch combined SDF

We can generate a single SDF file with all the poses using:

```{.python notest}
sim.docking.get_poses("poses.sdf")
```

## Filter results   

We can now filter this dataframe using any criteria we want. For example, we can only retain ligands that have a `pose_score` greater than `0.9` using:

```{.python notest}
df = df[df["pose_score"] > 0.9]
```

The filtered dataframe now only has ligands that matches the required criterion. 



## Filter Ligands by SMILES strings

We can create a new SDF file with only these ligands using:

```{.python notest}
from deeporigin.drug_discovery import Ligand

ligands = Ligand.from_sdf("..path/to/poses.sdf") # from docking 

smiles_strings = df["SMILES"]

ligands = [ligand for ligand in ligands if ligand.properties["SMILES"] in smiles_strings]
```
