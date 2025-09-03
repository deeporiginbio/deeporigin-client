# Configuration

## On a Deep Origin workstation

On a Deep Origin workstation, no configuration is needed! Within a workstation, the Deep Origin CLI and Python client are automatically configured.

## On your local computer

To run this package outside of a Deep Origin workstation (for example, on your own computer), first you need to configure this package. After installing this package, set your organization key and environment.

=== "python"


    ```python
    from deeporigin import config
    config.set_value("org_key", "org-key")
    config.set_value("env", "prod")  # or "staging" / "edge"
    ```


=== "CLI"

    ```bash
    deeporigin config set org_key [org-key]
    deeporigin config set env prod
    ```

## View configuration

To view the configuration for this package, run:

=== "python"

    ```python
    from deeporigin import config
    config.get_value()

    ```

=== "CLI"

    ```bash
    deeporigin config show
    ```

    This will display a table such as below:

    ```
    ╭─────────────────────────────┬─────────────────────────────────────╮
    │ Variable                    │ Value                               │
    ├─────────────────────────────┼─────────────────────────────────────┤
    │ org_key                     │ deeporigin                          │
    │ env                         │ prod                                │
    ╰─────────────────────────────┴─────────────────────────────────────╯
    ```
