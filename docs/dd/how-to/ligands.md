This document describes how to work with ligands (molecules) and use them in Deep Origin tools. 

There are two classes that help you work with ligands:

- `Ligand`
- `LigandSet`

The [`Ligand` class](../ref/structures.md#src.drug_discovery.structures.Ligand) is the primary way to work with ligands in Deep Origin.

## Constructing a Ligand or LigandSet

### Single Ligand from a SDF file

A single Ligand can be constructed from a file:

```python
from deeporigin.drug_discovery import Ligand, EXAMPLE_DATA_DIR

ligand = Ligand.from_sdf(EXAMPLE_DATA_DIR / "brd-2.sdf")
```

### Many ligands from a SDF file

A LigandSet can be constructed from a SDF File:

```python
from deeporigin.drug_discovery import LigandSet, EXAMPLE_DATA_DIR

ligands = Ligand.from_sdf(EXAMPLE_DATA_DIR / "ligands-brd-all.sdf")
```

### From a SMILES string

A ligand can be constructed from a SMILES string, which is a compact way to represent molecular structures:

```python
from deeporigin.drug_discovery import Ligand


# Basic usage with just a SMILES string
ligand = Ligand.from_smiles(smiles="CCO")  # Ethanol

# With additional parameters
ligand = Ligand.from_smiles(
    smiles="c1ccccc1",  # Benzene
    name="Benzene",     # Optional name for the ligand
    save_to_file=False  # Optional: whether to save the ligand to file
)
```

The `from_smiles` constructor:

- Takes a SMILES string as input
- Optionally accepts a name for the ligand
- Optionally accepts a `save_to_file` parameter to control file persistence
- Automatically validates the SMILES string and creates a proper molecular representation
- Returns a `Ligand` instance that can be used for further operations

!!! note "SMILES Validation"
    The constructor will raise an exception if the provided SMILES string is invalid or cannot be parsed into a valid molecule.


### From a Chemical Identifier

You can create a ligand from common chemical identifiers (like PubChem names, common names, or drug names). This is particularly useful when working with well-known biochemical molecules:

```python
from deeporigin.drug_discovery import Ligand

# Create ligands from common biochemical names
atp = Ligand.from_identifier(
    identifier="ATP",  # Adenosine triphosphate
    name="ATP"
)

serotonin = Ligand.from_identifier(
    identifier="serotonin",  # 5-hydroxytryptamine (5-HT)
    name="Serotonin"
)
```

The `from_identifier` constructor:

- Accepts common chemical names and identifiers
- Automatically resolves the identifier to a molecular structure
- Creates a 3D conformation of the molecule
- Particularly useful for well-known biochemical molecules like:
    - Nucleotides (ATP, ADP, GTP, etc.)
    - Neurotransmitters (serotonin, dopamine, etc.)
    - Drug molecules (by their generic names)
    - Common metabolites and cofactors

!!! note "Identifier Resolution"
    The constructor will attempt to resolve the identifier using chemical databases. If the identifier cannot be resolved, it will raise an exception.

### From an RDKit Mol object

If you're working with RDKit molecules directly, you can create a Ligand from an RDKit Mol object:

```python
from deeporigin.drug_discovery import Ligand
from rdkit import Chem

# Create an RDKit molecule
mol = Chem.MolFromSmiles("CCO")  # Ethanol

# Convert to a Ligand
ligand = Ligand.from_rdkit_mol(
    mol=mol,
    name="Ethanol",  # Optional name for the ligand
    save_to_file=False  # Optional: whether to save the ligand to file
)
```

This is particularly useful when you're working with RDKit's molecular manipulation functions and want to convert the results into a Deep Origin Ligand object for further processing or visualization.


The method will:

- Read the CSV file using pandas
- Extract SMILES strings from the specified column
- Create a Ligand instance for each valid SMILES
- Store all other columns as properties in each Ligand instance
- Skip any rows with empty or invalid SMILES strings

!!! note "Error Handling"
    The method will raise:
    - `FileNotFoundError` if the CSV file does not exist
    - `ValueError` if the specified SMILES column is not found in the CSV file

### From a CSV file

You can also create a `LigandSet` from a CSV file containing SMILES strings and optional properties:

```python
from deeporigin.drug_discovery import LigandSet, EXAMPLE_DATA_DIR

ligands = LigandSet.from_csv(
    file_path=EXAMPLE_DATA_DIR / "ligands.csv",
    smiles_column="SMILES"  # Optional, defaults to "smiles"
)
```

## Visualization

### Visualizing individual ligands

??? warning "Browser support"
    These visualizations work best on Google Chrome. We are aware of issues on other browsers, especially Safari on macOS.

A ligand object can be visualized using `show`:

```python
from deeporigin.drug_discovery import Ligand

ligand = Ligand.from_identifier("serotonin")

ligand.show()
```

If a ligand is backed by a SDF file, a 3D visualization will be shown, similar to:

A visualization such as this will be shown:

<iframe 
    src="./ligand.html" 
    width="100%" 
    height="600" 
    style="border:none;"
    title="ligand visualization"
></iframe>

!!! info "Jupyter notebook required"
    Visualizations such as these require this code to be run in a jupyter notebook. We recommend using [these instructions](../../install.md) to install Jupyter.


If a ligand is not backed by a SDF file, a 2D visualization will be shown:

![](../../images/ligand.png)

### Visualizing LigandSets

## Predicting ADMET Properties

You can predict ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity) properties for a ligand using the `admet_properties` method:

```{.python notest}
# Predict ADMET properties
properties = ligand.admet_properties()
```

The method returns a dictionary containing various ADMET-related predictions:

```python
{
    'smiles': 'Cn1c(=O)n(Cc2ccccc2)c(=O)c2c1nc(SCCO)n2Cc1ccccc1',
    'properties': {
        'logS': {'value': -4.004},  # Aqueous solubility
        'logP': 3.686,             # Partition coefficient
        'logD': 2.528,             # Distribution coefficient
        'hERG': {'probability': 0.264},  # hERG inhibition risk
        'ames': {'probability': 0.213}, # Ames mutagenicity
        'cyp': {                                # Cytochrome P450 inhibition
            'probabilities': {
                'cyp1a2': 0.134,
                'cyp2c9': 0.744,
                'cyp2c19': 0.853,
                'cyp2d6': 0.0252,
                'cyp3a4': 0.4718
            }
        },
        'pains': {                              # PAINS (Pan Assay Interference Compounds)
            'has_pains': None,
            'pains_fragments': []
        }
    }
}
```

The predicted properties are automatically stored in the ligand's properties dictionary and can be accessed later using the `get_property` method:

```{.python notest}
# Access a specific property
logP = ligand.get_property('logP')
```

!!! note "Property Storage"
    All predicted properties are automatically stored in the ligand's properties dictionary and can be accessed at any time using the `get_property` method.

## Exporting ligands

### To SDF files

To write a ligand to a SDF file, use:

```python
from deeporigin.drug_discovery import Ligand

ligand = Ligand.from_identifier("serotonin")
ligand.to_sdf()
```

### To mol files

To write a ligand to a mol file, use:

```python
from deeporigin.drug_discovery import Ligand

ligand = Ligand.from_identifier("serotonin")
ligand.to_mol()
```

### To PDB files

To write a ligand to a PDB file, use:

```python
from deeporigin.drug_discovery import Ligand

ligand = Ligand.from_identifier("serotonin")
ligand.to_pdb()
```


### To Pandas DataFrames

To convert a LigandSet to a Pandas DataFrame, use:

```python
from deeporigin.drug_discovery import LigandSet, EXAMPLE_DATA_DIR

ligands = LigandSet.from_csv(
    file_path=EXAMPLE_DATA_DIR / "ligands.csv",
    smiles_column="SMILES"  # Optional, defaults to "smiles"
)
df = ligands.to_dataframe()
```




This will read the CSV, create a `Ligand` for each valid row, and store them in the set.


This will display a table of ligands with their structures and properties, either in a Jupyter notebook or as a pandas DataFrame, depending on the environment.

### From an SDF file

You can create a `LigandSet` from an SDF file containing one or more molecules:

```python
from deeporigin.drug_discovery import LigandSet

ligand_set = LigandSet.from_sdf(
    file_path="molecules.sdf"
)
```

This will load all valid molecules in the SDF file as `Ligand` objects in the set. If the file contains no valid molecules, an empty set will be returned.

