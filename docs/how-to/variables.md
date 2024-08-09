# Install variables and secrets 

[Variables and secrets](https://docs.deeporigin.io/docs/os/variables-secrets) captured in the Deep Orgin OS, either
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

!!! info "Only on Deep Origin workstations"
    This function of the CLI only works on Deep Origin workstations. This function will not work on your local computer.
