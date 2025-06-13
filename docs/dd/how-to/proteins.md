This document describes how to work with proteins and use them in Deep Origin tools. 

## The `Protein` class

The [`Protein` class](../ref/structures.md#src.drug_discovery.structures.Protein)  is the primary way to work with proteins in Deep Origin.

## Constructing a protein

### From a file

A protein can be constructed from a file:

```python
from deeporigin.drug_discovery import Protein, EXAMPLE_DATA_DIR
protein = Protein.from_file(EXAMPLE_DATA_DIR / "brd.pdb")
```

### From a PDB ID

A protein can also be constructed from a PDB ID:


```python
from deeporigin.drug_discovery import Protein

protein = Protein.from_pdb_id("1EBY")
```


## Visualizing a protein

??? warning "Browser support"
    These visualizations work best on Google Chrome. We are aware of issues on other browsers, especially Safari on macOS.

A protein object can be visualized using `show`:

```python
from deeporigin.drug_discovery import Protein

protein = Protein.from_pdb_id("1EBY")
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

## Modifying and preparing a protein

### Removing specific residues

You can remove specific residue names from a protein structure using the `remove_resnames` method:

```python
from deeporigin.drug_discovery import Protein

protein = Protein.from_pdb_id("1EBY")

# Remove water molecules (HOH) and ions (NA, CL)
protein.remove_resnames(exclude_resnames=["HOH", "NA", "CL"])
```

!!! note "Residue name format"
    Residue names in PDB files are always uppercase (e.g., "HOH" for water, "NA" for sodium, "CL" for chloride).

This method modifies the protein structure in place by removing the specified residue names. If no residue names are provided, the protein structure remains unchanged.

### Removing HETATM records

You can remove HETATM records from the protein structure using the `remove_hetatm` method:

```python

from deeporigin.drug_discovery import Protein

protein = Protein.from_pdb_id("1EBY")

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

from deeporigin.drug_discovery import Protein

protein = Protein.from_pdb_id("1EBY")

protein.remove_water()
```

This method modifies the protein structure in place by removing all water molecules (residue name "HOH"). Unlike `remove_resnames`, this method does not return a new protein structure but instead modifies the existing one.


### Removing Specific Residues

You can remove specific residue types while keeping others:

```python

from deeporigin.drug_discovery import Protein

protein = Protein.from_pdb_id("1EBY")


# Remove specific residue names
protein.remove_resnames(exclude_resnames=["HOH", "SO4"])

# Remove only water molecules
protein_no_water = protein.remove_water()
```

### Chain Selection

You can select specific chains from the protein structure:

```python notest
# Select a single chain
chain_a = protein.select_chain('A')

# Select multiple chains
chains_ab = protein.select_chains(['A', 'B'])
```

## Best Practices

1. Always visualize the structure before and after preparation to ensure the desired changes were made
2. When working with metalloproteins, use `remove_hetatm()` with appropriate `keep_resnames` or `remove_metals` parameters
3. For multi-chain proteins, consider selecting specific chains before preparation
4. Save the prepared structure using `to_pdb()` if you need to use it later:

```python notest
# Save the prepared structure
protein.to_pdb("prepared_protein.pdb")
```

## Common Use Cases

### Preparing a Protein for Docking

```python notest
# Load and prepare protein
protein = Protein.from_pdb_id("1EBY")
protein.remove_water()  # Remove water molecules
protein.remove_hetatm(keep_resnames=['ZN'])  # Keep important cofactors
protein.to_pdb("docking_ready.pdb")
```

### Working with Metalloproteins

```python notest
# Load and prepare metalloprotein
protein = Protein.from_pdb_id("1XYZ")
protein.remove_hetatm(remove_metals=['ZN', 'MG'])  # Keep metal ions
protein.show()
```

### Multi-chain Protein Preparation

```python notest
# Load and prepare multi-chain protein
protein = Protein.from_pdb_id("1ABC")
chains_ab = protein.select_chains(['A', 'B'])  # Select chains A and B
chains_ab.remove_water()  # Remove water from selected chains
chains_ab.show()
```
