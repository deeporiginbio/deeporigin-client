# Visualizing ABFE Trajectories

This guide explains how to visualize molecular dynamics trajectories from Absolute Binding Free Energy (ABFE) simulations in Deep Origin.

## Overview

ABFE simulations generate molecular dynamics trajectories that show how ligands interact with proteins over time. Visualizing these trajectories can provide valuable insights into binding mechanisms, protein-ligand interactions, and conformational changes.

Deep Origin provides tools to easily visualize these trajectories using the `show_trajectory` method in the `ABFE` class.

## Prerequisites

- A completed ABFE simulation run
- Valid ligand ID from your simulation
- The Deep Origin Python package properly installed and configured

## Visualizing Trajectories

### Using the `show_trajectory` Method

The `show_trajectory` method allows you to visualize either the molecular dynamics (md) phase or the binding phase of an ABFE run.

Before using this method, you'll need to have a properly set up Complex object with completed ABFE calculations. For a complete walkthrough of setting up your environment, creating a Complex object, and running ABFE calculations, please refer to the [Getting Started](../tutorial/getting-started.md) tutorial.

Once you've completed the steps in the tutorial and have run your ABFE calculations, you can return to this guide to visualize the trajectories.


### Behind the Scenes

When you call `show_trajectory`, the following happens:
1. The method retrieves result files for the specified ABFE run
2. It extracts the trajectory files from a ZIP archive (if not already extracted)
3. It generates a 3D visualization using the `deeporigin_molstar` library
4. The visualization is displayed in your Jupyter notebook

## Examples

### Visualizing Molecular Dynamics Step

```{.python notest}
# Show the molecular dynamics trajectory
sim.abfe.show_trajectory(ligands=sim.ligands[0], step="md")
```

### Visualizing Binding Step

```{.python notest}
# Show the binding trajectory for the default window (window 1)
sim.abfe.show_trajectory(ligand_id="ligand1", step="binding")

# Show the binding trajectory for a specific window
sim.abfe.show_trajectory(ligand_id="ligand1", step="binding", window=5)
```

When visualizing the binding step, you can specify which window to visualize using the `window` parameter. The windows represent different stages of the binding process. If you specify an invalid window number, the method will show you a list of available windows.

<iframe 
    src="../../images/traj.html" 
    width="100%" 
    height="600" 
    style="border:none;"
    title="Protein visualization"
></iframe>


## Troubleshooting

If you encounter issues when visualizing trajectories:

- Ensure your ABFE run completed successfully
- Verify that the ligand ID is correct
- Check that the specified step ("md" or "binding") is available
- Ensure you have sufficient disk space for extracting trajectory files

## Additional Resources

- For more information on running ABFE simulations, see the [ABFE documentation](../tutorial/abfe.md)
- For details on the `Complex` class, refer to the [API reference](../ref/abfe.md)