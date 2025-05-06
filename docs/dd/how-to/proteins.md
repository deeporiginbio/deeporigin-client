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

## Modifying proteins

### Removing specific residues

You can remove specific residue names from a protein structure using the `remove_resnames` method:

```python
# Remove water molecules (HOH) and ions (NA, CL)
protein.remove_resnames(exclude_resnames=["HOH", "NA", "CL"])
```

!!! note "Residue name format"
    Residue names in PDB files are always uppercase (e.g., "HOH" for water, "NA" for sodium, "CL" for chloride).

This method modifies the protein structure in place by removing the specified residue names. If no residue names are provided, the protein structure remains unchanged.

### Removing HETATM records

You can remove HETATM records from the protein structure using the `remove_hetatm` method:

```python
# Remove all HETATM records except water and zinc
protein.remove_hetatm(keep_resnames=["HOH"], remove_metals=["ZN"])
```

This method modifies the protein structure in place by removing HETATM records (heteroatoms) from the structure. You can:
- Keep specific residues by providing their names in `keep_resnames`
- Keep specific metals by providing their names in `remove_metals` (these metals will be excluded from removal)
- If no arguments are provided, all HETATM records will be removed

!!! note "Metal names"
    Metal names should be provided in uppercase (e.g., "ZN" for zinc, "FE" for iron).

### Removing water molecules

You can remove all water molecules from a protein structure using the `remove_water` method:

```python
protein.remove_water()
```

This method modifies the protein structure in place by removing all water molecules (residue name "HOH"). Unlike `remove_resnames`, this method does not return a new protein structure but instead modifies the existing one.
