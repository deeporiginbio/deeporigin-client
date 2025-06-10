# Docking

This document describes how to [dock :octicons-link-external-16:](https://en.wikipedia.org/wiki/Docking_(molecular)) a set of ligands to a protein  using Deep Origin tools. 

## Prerequisites

We assume that we have an initialized and configured `Complex` object:

```python
from deeporigin.drug_discovery import Complex, EXAMPLE_DATA_DIR
sim = Complex.from_dir(EXAMPLE_DATA_DIR) # or replace with your folder
```

For more details on how to get started, see [:material-page-previous: Getting Started ](./getting-started.md).


## Find pockets

First we find pockets in the protein using:

```py
pockets = sim.protein.find_pockets(pocket_count=1)
```

We can visualize the pocket using:

```py
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


```py
pocket = pockets[0] # or choose as needed
sim.docking.run(pocket=pocket)

```

### Using residue ID

!!! info "Coming Soon"
    Ability to dock ligands using a residue ID is coming soon.

### Using pocket center and box size


To dock all ligands to the protein, parallelizing and batching across all ligands, we do the following:


```python
job = sim.docking.run(
    box_size=(15, 15, 15),      # match to your protein
    pocket_center=(13, -6, 22), # match to your protein
)
```

??? info "Controlling batch size"

    By default, all ligands are docked in batches of 32 ligands. 

    This can be controlled in two ways. First, you can control the batch size using the `batch_size` parameter.

    ```python
    sim.dock(
        batch_size=32,
        ... 
    )
    ```

    You can also specify the number of workers using:

    ```python
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

```py
job.watch()
```

Doing so creates a widget that automatically updates and monitors a job as long as its running. 

To stop watching a job, do:

```py
job.stop_watching()
```

## Results

### Viewing results

After completion of docking, we can view results using:

```python
sim.docking.show_results()
```  

This shows a table similar to:

![Docking results](../../images/tools/docking-results.png)

### Viewing docked poses

To view the docked poses of all ligands in the complex, use:

```python
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

```python
df = sim.docking.get_results()
```

### Exporting a SDF with docked poses

To export a SDF with docked poses, use:

```python
sim.docking.get_poses("/path/to/output.sdf")
```

This generates a SDF file with the docked poses for all ligands in the Complex. 