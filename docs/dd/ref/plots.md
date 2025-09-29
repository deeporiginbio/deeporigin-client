# Plots

The `deeporigin.plots` module provides visualization functions for drug discovery data, including interactive scatter plots with molecular structure visualization.

## Functions

::: src.plots.scatter
    options:
      show_signature_annotations: true
      show_root_heading: true
      show_root_full_path: false
      heading_level: 3

## Examples

### Basic Scatter Plot with Molecule Images

```python
import numpy as np
from deeporigin.plots import scatter

# Sample data
x = np.array([1, 2, 3, 4, 5])
y = np.array([2, 4, 6, 8, 10])
smiles_list = [
    "CCO",           # ethanol
    "CC(=O)O",       # acetic acid
    "c1ccccc1",      # benzene
    "CCN(CC)CC",     # triethylamine
    "CC(C)O"         # isopropanol
]

# Create and display scatter plot
scatter(
    x=x, 
    y=y, 
    smiles_list=smiles_list, 
    x_label="X Coordinate", 
    y_label="Y Coordinate",
    title="Molecule Analysis Plot"
)
```

### Working with Drug Discovery Data

```{.python notest}
import pandas as pd
from deeporigin.plots import scatter

# Load data from a CSV file with SMILES and properties
df = pd.read_csv("ligand_data.csv")

# Create and display scatter plot of molecular weight vs logP
scatter(
    x=df["molecular_weight"].values,
    y=df["logp"].values,
    smiles_list=df["smiles"].tolist(),
    x_label="Molecular Weight (Da)",
    y_label="LogP",
    title="Drug Discovery: Molecular Properties Analysis"
)
```

## Features

- **Interactive Hover**: Hover over any point to see the molecular structure image
- **SMILES Validation**: Automatically filters out invalid SMILES strings
- **High-Quality Images**: Generates 200x200 pixel molecular structure images
- **Responsive Design**: Follows mouse movement for optimal user experience
- **Error Handling**: Gracefully handles invalid SMILES and rendering errors

## Requirements

The plots module requires the following optional dependencies:

- `bokeh`: For interactive plotting
- `rdkit`: For molecular structure rendering

Install them with:

```bash
pip install deeporigin[plots,tools]
```

## Notes

- Invalid SMILES strings are automatically filtered out
- If all SMILES strings are invalid, a `ValueError` is raised
- The function returns a Bokeh figure object that can be further customized
- Molecular images are generated at 200x200 pixels for optimal display
