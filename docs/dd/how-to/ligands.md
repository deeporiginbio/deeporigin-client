This document describes how to work with ligands (molecules) and use them in Deep Origin tools. 

## The `Ligand` class

The [`Ligand` class](../ref/structures.md#src.drug_discovery.structures.Ligand) is the primary way to work with ligands in Deep Origin.

## Constructing a ligand

### From a file

A ligand can be constructed from a file:

```python
from deeporigin.drug_discovery import Ligand, EXAMPLE_DATA_DIR

ligand = Ligand.from_sdf(EXAMPLE_DATA_DIR / "brd-2.sdf")
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

### From a CSV file

You can create multiple ligands at once from a CSV file using the `from_csv` class method. This is useful when you have a dataset of molecules with their SMILES strings and associated properties. The `from_csv` method:

- Requires a path to the CSV file and the name of the column containing SMILES strings
- Automatically extracts all non-SMILES columns as properties
- Skips rows with empty or invalid SMILES
- Returns a list of `Ligand` objects

This approach is ideal for processing large datasets of molecules where each row represents a different compound.

#### Basic Usage

For the simplest case, just specify the file path and which column contains the SMILES strings:

```{.python notest}
from deeporigin.drug_discovery import Ligand

# Basic usage - just extracting SMILES from a column
ligands = Ligand.from_csv(
    file_path="molecules.csv",
    smiles_column="SMILES"  # Optional, defaults to "smiles"
)
```

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

## Visualizing a ligand

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