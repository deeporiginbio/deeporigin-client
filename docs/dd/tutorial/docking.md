# Docking

This document describes how to [dock :octicons-link-external-16:](https://en.wikipedia.org/wiki/Docking_(molecular)) a set of ligands to a protein  using Deep Origin tools. 

## Prerequisites

We assume that we have an initialized and configured a `Complex` object:

```python
from deeporigin.drug_discovery import Complex, BRD_DATA_DIR

sim = Complex.from_dir(BRD_DATA_DIR) 
```

For more details on how to get started, see [:material-page-previous: Getting Started ](./getting-started.md).


## Find pockets

First we find pockets in the protein using:

```{.python notest}
pockets = sim.protein.find_pockets(pocket_count=1)
```

We can visualize the pocket using:

```{.python notest}
sim.protein.show(pockets=pockets)
```

You should see something along the lines of:

<iframe 
    src="./pockets.html" 
    width="100%" 
    height="650" 
    style="border:none;"
    title="Protein visualization"
></iframe>

## Starting a docking run

### Using pocket

To dock all ligands in the complex to the protein, using the pocket we found, we can do:


```{.python notest}
pocket = pockets[0] # or choose as needed
sim.docking.run(pocket=pocket)

```

### Using residue ID

!!! info "Coming Soon"
    Ability to dock ligands using a residue ID is coming soon.

### Using pocket center and box size


To dock all ligands to the protein, parallelizing and batching across all ligands, we do the following:


```{.python notest}
job = sim.docking.run(
    box_size=(15, 15, 15),      # match to your protein
    pocket_center=(13, -6, 22), # match to your protein
)
```

??? info "Controlling batch size"

    By default, all ligands are docked in batches of 32 ligands. 

    This can be controlled in two ways. First, you can control the batch size using the `batch_size` parameter.

    ```{.python notest}
    sim.dock(
        batch_size=32,
        ... 
    )
    ```

    You can also specify the number of workers using:

    ```{.python notest}
    sim.dock(
        n_workers=2,
        ...
    )
    ```

    You can specify either the number of workers or the batch size, but not both. 


This queues up tasks on Deep Origin. When it completes, the results of docking can be viewed.

## Viewing status of docking

A job object is returned from `docking.run`. This job object can be inspected to show the status of the job when created. 

The job object can also be used monitor a job as it completes:

```{.python notest}
job.watch()
```

Doing so creates a widget that automatically updates and monitors a job as long as its running. 

To stop watching a job, do:

```{.python notest}
job.stop_watching()
```

## Results

### Viewing results

After completion of docking, we can view results using:

```{.python notest}
sim.docking.show_results()
```  

This shows a table similar to:

![Docking results](../../images/tools/docking-results.png)

### Viewing docked poses

To view the docked poses of all ligands in the complex, use:

```{.python notest}
sim.docking.show_poses()
```

<iframe 
    src="./docked-poses.html" 
    width="100%" 
    height="650" 
    style="border:none;"
    title="Protein visualization"
></iframe>


### Exporting for further analysis

To obtain the raw dataframe for further analysis, use:

```{.python notest}
df = sim.docking.get_results()
```

### Get all docked poses

Use the `get_poses` function to return a [`LigandSet`](../ref/ligandset.md) with all docked poses

```{.python notest}
poses = sim.docking.get_poses()
```

and use a `LigandSet` method to write to SDF:

```{.python notest}
poses.to_sdf("/path/to/docked/ligands.sdf")
```

