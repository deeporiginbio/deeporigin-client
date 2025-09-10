This document describes how to find pockets in a Protein using the Deep Origin Pocket Finder

# Find pockets

First, we create a Protein, for example, using:

```python
from deeporigin.drug_discovery import  Protein, BRD_DATA_DIR
protein = Protein.from_file(BRD_DATA_DIR / "brd.pdb")
```
We can then find pockets in this protein using the Pocket Finder tool:

```{.python notest}
pockets = protein.find_pockets(pocket_count=1)
```

Pockets can be visualized using:

```{.python notest}
protein.show(pockets=pockets)
```

You should see something like:

<iframe 
    src="../../images/pockets.html" 
    width="100%" 
    height="650" 
    style="border:none;"
    title="Protein visualization"
></iframe>