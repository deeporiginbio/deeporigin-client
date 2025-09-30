This document describes how to use the outputs of Docking as inputs to FEP tools.

## Assumptions

We assume that you have 

1. [:material-page-previous: created a `Complex` object](../tutorial/getting-started.md)
2. [:material-page-previous: run Docking](../tutorial/docking.md)

## Pick poses to use

If you have used docking on a number of ligands, you can find the best poses for each ligand using:

```{.python notest}
poses = sim.docking.get_poses()
poses = poses.filter_top_poses()
```

These poses can then be used as ligands:

```{.python notest}
sim.ligands = poses
```

and ABFE can be run on them as described [:material-page-previous: here](../tutorial/abfe.md)