# Configuration

## On a Deep Origin workstation

On a Deep Origin workstation, no configuration is needed! Within a workstation, the Deep Origin CLI and Python client are automatically configured.

## On your local computer

To run this package outside of a Deep Origin workstation (for example, on your own computer), first you need to configure this package. After installing this package, run the following to configure your organization, replacing `org-id` with the ID of the Deep Origin organization that you would like to work with.

```bash
deeporigin config set organization_id [org-id]
```
