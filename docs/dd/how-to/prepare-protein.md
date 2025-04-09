# Preparing Proteins

This guide demonstrates how to prepare protein structures for drug discovery workflows using the `Protein` class.

## Basic Protein Preparation

The most common preparation step is removing heterogeneous atoms (HETATM records) from the protein structure. Here's how to do it:

```python
from deeporigin.drug_discovery import Protein

# Load a protein structure from PDB ID
protein = Protein.from_pdb_id("1EBY")

# Visualize the original structure (includes water and heteroatoms)
protein.show()

# Remove all heterogeneous atoms
protein.remove_hetatm()

# Visualize the cleaned structure
protein.show()
```

## Advanced Preparation Options

### Removing Specific Residues

You can remove specific residue types while keeping others:

```python
# Remove specific residue names
protein.remove_resnames(exclude_resnames=["HOH", "SO4"])

# Remove only water molecules
protein_no_water = protein.remove_water()
```

### Chain Selection

You can select specific chains from the protein structure:

```python
# Select a single chain
chain_a = protein.select_chain('A')

# Select multiple chains
chains_ab = protein.select_chains(['A', 'B'])
```

### Keeping Specific Heteroatoms

When removing heterogeneous atoms, you can specify which ones to keep:

```python
# Remove all HETATM records except for specified residues
protein.remove_hetatm(keep_resnames=['ZN', 'MG'])

# Remove all HETATM records except for specific metals
protein.remove_hetatm(remove_metals=['ZN', 'MG'])
```

## Best Practices

1. Always visualize the structure before and after preparation to ensure the desired changes were made
2. When working with metalloproteins, use `remove_hetatm()` with appropriate `keep_resnames` or `remove_metals` parameters
3. For multi-chain proteins, consider selecting specific chains before preparation
4. Save the prepared structure using `to_pdb()` if you need to use it later:

```python
# Save the prepared structure
protein.to_pdb("prepared_protein.pdb")
```

## Common Use Cases

### Preparing a Protein for Docking

```python
# Load and prepare protein
protein = Protein.from_pdb_id("1EBY")
protein.remove_water()  # Remove water molecules
protein.remove_hetatm(keep_resnames=['ZN'])  # Keep important cofactors
protein.to_pdb("docking_ready.pdb")
```

### Working with Metalloproteins

```python
# Load and prepare metalloprotein
protein = Protein.from_pdb_id("1XYZ")
protein.remove_hetatm(remove_metals=['ZN', 'MG'])  # Keep metal ions
protein.show()
```

### Multi-chain Protein Preparation

```python
# Load and prepare multi-chain protein
protein = Protein.from_pdb_id("1ABC")
chains_ab = protein.select_chains(['A', 'B'])  # Select chains A and B
chains_ab.remove_water()  # Remove water from selected chains
chains_ab.show()
```
