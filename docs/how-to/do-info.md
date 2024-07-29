

## View Deep Origin Config Info


To view configuration information on Deep Origin, run:



```bash
deeporigin config
```

and you will see a table that looks like this:

```
╭─────────────────────────────┬─────────────────────────────────────╮
│ Variable                    │ Value                               │
├─────────────────────────────┼─────────────────────────────────────┤
│ organization_id             │ likely-aardvark-ewo                 │
│ bench_id                    │ average-possum-3x3                  │
│ env                         │ us-west-1                           │
│ api_endpoint                │ https://os.prod.deeporigin.io/api   │
│ nucleus_api_route           │ nucleus-api/api/                    │
│ graphql_api_route           │ api/graphql/                        │
│ auth_domain                 │ https://formicbio-dev.us.auth0.com  │
│ auth_device_code_endpoint   │ oauth/device/code/                  │
│ auth_token_endpoint         │ oauth/token/                        │
│ auth_audience               │ https://os.deeporigin.io/api        │
│ auth_grant_type             │ urn:ietf:params:device_code         │
│ auth_client_id              │ <secret>                            │
│ auth_client_secret          │ <secret>                            │
│ api_tokens_filename         │ ~/api_tokens                        │
│ variables_cache_filename    │ ~/variables.yml                     │
│ feature_flags               │                                     │
╰─────────────────────────────┴─────────────────────────────────────╯
```

## View Deep Origin Context Information

Context information from Deep Origin can be displayed using the CLI:

```bash
deeporigin context
```

and you will see a table that looks like this:


```bash
Bench ID: average-possum-3x3
User ID: None
Organization ID: likely-aardvark-ewo
Environment: edge
Debug: False
```


!!! warning "Only supported on ComputeBenches"
    Note that this command is only supported on ComputeBenches. If you run this command outside of a ComputeBench, you will likely see something like this:

    ```bash
    Bench ID: None
    User ID: None
    Organization ID: None
    Environment: None
    Debug: False
    ```