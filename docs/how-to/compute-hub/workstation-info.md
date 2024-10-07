# Get information about your Deep Origin workstation

Inside a Deep Origin workstation, run the following command to get information about the workstation:

```bash
deeporigin context
```

This will display a table similar to below:

```bash
Workstation ID: average-possum-3x3
User ID: None
Organization ID: likely-aardvark-ewo
Environment: edge
Debug: False
```

!!! warning "Only supported on Deep Origin workstations"
    Note that this command is only supported on Deep Origin workstations. If you run this command outside of a workstation, you will likely see this output:

    ```bash
    Workstation ID: None
    User ID: None
    Organization ID: None
    Environment: None
    Debug: False
    ```