# Docking

This document describes how to [dock](https://en.wikipedia.org/wiki/Docking_(molecular)) a set of ligands to a protein  using Deep Origin tools. 

## Prerequisites 

Make sure you have [installed](../../install.md), [configured](../../configure.md), and [authenticated](../../how-to/auth.md) with the Deep Origin python client.

!!! tip "Use uv to install `deeporigin`" 
    We strongly recommend using [uv](https://docs.astral.sh/uv/) to install Deep Origin in a separate project, following instructions [here](../../install.md#using-uv-to-set-up-deeporigin-on-your-computer). This gives you all the dependencies you need, together with a stand-alone install of Jupyter Lab that you can use with `deeporigin`




## Required input files 

You will need to have the following input files on your local computer:

1. A set of ligand files, in SDF format.
2. A protein PDB file 


## Running Bulk Docking

!!! tip "Jupyter notebooks"
    It is assumed that you are working in a Jupyter notebook (or similar IPython environment). This makes it easier to run the workflow, and some functions assume that you are in a Jupyter notebook.

First, we import necessary functions and modules from the `deeporigin` package:

```python
from deeporigin.tools import docking
```

We then specify where our input files are:

```python
input_dir = "/path/to/input/dir"
```

### Initialization

Here, we data structures to store input and intermediate files, and upload input files to Deep Origin:

```python
sim = docking.Docking.from_dir(input_dir)
```


### Connecting to Deep Origin

We can connect the instance we created to Deep Origin using

```python
sim.connect()
```


This function creates necessary databases on Deep Origin, uploads ligand and protein files if needed, and updates the status of runs on Deep Origin.





### Running bulk docking in parallel

To dock all ligands to the protein in bulk, parallelizing and batching across all ligands, we do the following:



```python
sim.dock(
    batch_size=30,
    box_size=[15, 15, 15],      # match to your protein
    pocket_center=[13, -6, 22], # match to your protein
)
```

This queues up tasks on Deep Origin. When it completes, the results of docking can be viewed.



### Results

After completion of bulk docking, we can view results using:

```python
sim.show_results()
```  

This shows a table similar to:

![Docking results](../../images/tools/docking-results.png)

To obtain the raw dataframe for further analysis, use:

```python
df = sim.get_results()
```