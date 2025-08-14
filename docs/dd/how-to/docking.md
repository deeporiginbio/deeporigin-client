
# Docking Ligands to a Protein

This document describes how to dock ligands to a Protein. 


## Prerequisites

- You have a prepared [Protein](./protein.md)
- You have a [Ligand or LigandSet](./ligands.md)
- You have found pockets in your Protein using the [Pocket Finder](./find-pockets.md)

## Docking a single `Ligand`

A single Ligand can be docked to a Protein using:

```{.python notest}
poses = protein.dock(
    pocket=pocket,
    ligand=ligand,
)
```

where `pocket` is a `Pocket` object generated using the [Pocket Finder Tool](find-pockets.md) :octicons-book-24: . 


### Viewing docked poses

Docked poses for that ligand can be viewed using:

```{.python notest}
protein.show(poses=poses)
```

You will see something similar to the following. Use the arrows to inspect individual poses. 

<iframe 
    src="../../images/1eby-docked-poses.html" 
    width="100%" 
    height="650" 
    style="border:none;"
    title="Docked poses of ligand in 1EBY"
></iframe>

### Viewing pose scores and binding energy

Every pose is assigned a pose score and a binding energy. These can be viewed using:

```{.python notest}
poses
```

To work with a dataframe containing this data, use:

```{.python notest}
df = poses.to_dataframe()
```

### Exporting poses to SDF

Poses can be saved to a SDF file using:


```{.python notest}
poses.to_sdf()
```

## Docking a `LigandSet` 

### Using Batch Jobs

!!! tip "Tutorial"
    Follow [the tutorial](../tutorial/docking.md) on how to dock ligands using Batch Jobs. This is best suited for large jobs with 100+ ligands. 

### Using Functions

Several ligands in a LigandSet can be docked to a Protein using:

```{.python notest}
poses = protein.dock(
    ligands=ligands,
    pocket=pocket,
)
```

`poses` contains all poses for all ligands in the LigandSet. To filter poses to keep only top poses, use:

```{.python notest}
poses = poses.filter_top_poses()
```

These poses can by visualized as before:

```{.python notest}
protein.show(poses=poses)
```


## Constrained Docking

We can use constrained docking to dock a Ligand to a Protein while constraining certain atoms to certain locations.

Typically, these constraints are computed a reference docked pose for another Ligand, using a Maximum Common Substructure (MCS) shared across Ligands. 

