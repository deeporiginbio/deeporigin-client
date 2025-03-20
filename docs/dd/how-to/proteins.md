This document describes how to work with proteins and use them in Deep Origin tools. 

## The `Protein` class

The [`Protein` class](../ref/chemistry.md#src.chemistry.Protein) is the primary way to work with proteins in Deep Origin.

## Constructing a protein

### From a file

A protein can be constructed from a file:

```python
from deeporigin.chemistry import Protein
protein = Protein("/path/to/pdb")
```

### From a PDB ID

A protein can also be constructed from a PDB ID:


```python
from deeporigin.chemistry import Protein
protein = Protein(pdb_id="1EBY")
```


## Visualizing a protein

??? warning "Browser support"
    These visualizations work best on Google Chrome. We are aware of issues on other browsers, especially Safari on macOS.

A protein object can be visualized using `show`:

```python
protein.show()
```

A visualization such as this will be shown:

<iframe 
    src="./protein.html" 
    width="100%" 
    height="600" 
    style="border:none;"
    title="Protein visualization"
></iframe>

!!! info "Jupyter notebook required"
    Visualizations such as these require this code to be run in a jupyter notebook. We recommend using [these instructions](../../install.md) to install Jupyter.


