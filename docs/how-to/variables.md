# Install variables and secrets into a workstation

Variables and secrets specified in the Deep Orgin platform, either
at the organization level or the account level, can be installed into
workstations by running the following command:

```bash
deeporigin variables install
```

This command will install your variables and secrets into your workstation and print a message similar to that below:

```bash
No variables were modified
No variables were added
No variables were deleted
1 variables were unmodified:
  EnvironmentVariable: foo
```

!!! warning "Only on Deep Origin workstations"
    This functionality is only meant for use in Deep Origin workstations.
