# Configuration

## On a Deep Origin workstation

On a Deep Origin workstation, no configuration is needed! Within a workstation, the Deep Origin CLI and Python client are automatically configured.

## On your local computer

To run this package outside of a Deep Origin workstation (for example, on your own computer), first you need to configure this package. After installing this package, run the following to configure your organization, replacing `org-id` with the ID of the Deep Origin organization that you would like to work with.

=== "python"


    ```python
    from deeporigin import config
    config.set("organization_id", "org-id")
    ```


=== "CLI"

    ```bash
    deeporigin config set organization_id [org-id]
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
    │ organization_id             │ likely-aardvark-ewo                 │
    │ bench_id                    │ average-possum-3x3                  │
    │ env                         │ us-west-1                           │
    │ api_endpoint                │ https://os.prod.deeporigin.io/api   │
    │ nucleus_api_route           │ datahub-api/api/                    │
    │ auth_domain                 │ https://formicbio-dev.us.auth0.com  │
    │ auth_device_code_endpoint   │ oauth/device/code/                  │
    │ auth_token_endpoint         │ oauth/token/                        │
    │ auth_audience               │ https://os.deeporigin.io/api        │
    │ auth_grant_type             │ urn:ietf:params:device_code         │
    │ auth_client_id              │ <secret>                            │
    │ auth_client_secret          │ <secret>                            │
    │ variables_cache_filename    │ ~/variables.yml                     │
    │ feature_flags               │                                     │
    ╰─────────────────────────────┴─────────────────────────────────────╯
    ```
