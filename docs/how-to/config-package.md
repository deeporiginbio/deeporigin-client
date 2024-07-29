# Configure this package

## View the configuration for this package

To view the configuration for this package, run:

```bash
deeporigin config show
```

This will display a table like below:

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

## Set a configuration variable

To set a variable, run:

```bash
deeporigin config set [variable-name] [variable-value]
```
