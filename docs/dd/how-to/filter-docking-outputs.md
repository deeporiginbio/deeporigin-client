This document describes how to filter the outputs of Docking based on various properties. 

Here we assume that you have constructed a `Complex` object and successfully run [Docking](../tutorial/docking.md). 
Following convention, we assume that the `Complex` object is called `sim`.

## Fetch Docking results

First, we get the results of Docking in a pandas DataFrame using:

```python
df = sim.docking.get_results()
```
## Fetch combined SDF

We can generate a single SDF file with all the poses using:

```python
sim.docking.get_poses("poses.sdf")
```

## Filter results   

We can now filter this dataframe using any criteria we want. For example, we can only retain ligands that have a `pose_score` greater than `0.9` using:

```python
df = df[df["pose_score"] > 0.9]
```

The filtered dataframe now only has ligands that matches the required criterion. 



## Filter SDF by SMILES strings

We can create a new SDF file with only these ligands using:

```python
from deeporigin.drug_discovery import chemistry

chemistry.filter_sdf_by_smiles(
    input_sdf_file="poses.sdf",
    output_sdf_file="filtered_ligands.sdf",
    keep_only_smiles=df["SMILES"],
)
```

We can visualize these ligands using:

```python
from deeporigin.drug_discovery import chemistry

chemistry.show_molecules_in_sdf_file("filtered_ligands.sdf")
```

A visualization similar to the following will be shown:

<iframe 
    src="./brd-ligands.html" 
    width="100%" 
    height="650" 
    style="border:none;"
    title="Ligands visualization"
></iframe>
