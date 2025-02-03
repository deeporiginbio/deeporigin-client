# Get information about workstations

To list all workstations on Deep Origin, use:

```python
from deeporigin.platform import api
api.get_workstations()
```

This returns a list of objects, where each object correspond to a workstation. A typical entry looks like:

```json
{
    "attributes": {
        "accessMethods": [
            {
                "icon": "/assets/icons/catalog-items/jupyterlab.svg",
                "id": "jupyterlab",
                "name": "JupyterLab"
            },
            {
                "icon": "/assets/icons/catalog-items/code-server.svg",
                "id": "code-server",
                "name": "VS Code (web)"
            }
        ],
        "accessSettings": {
            "publicKey": "ssh-ed25519 ... ",
            "ssh": true
        },
        "autoStopIdleCPUThreshold": 0,
        "autoStopIdleDuration": 30,
        "blueprint": "deeporigin/deeporigin-python:staging",
        "cloudProvider": {
            "region": "us-west-2",
            "vendor": "aws"
        },
        "clusterId": "3bb775e4-8be6-4936-a6b9",
        "created": "2024-10-05T17:01:06.840Z",
        "description": "dfd",
        "drn": "drn:...",
        "enableAutoStop": true,
        "name": "forthcoming-tyrannosaurus-8fd",
        "nextUserActions": [
            "DELETE"
        ],
        "orgHandle": "deeporigin-com",
        "requestedResources": {
            "cpu": 8,
            "gpu": 0,
            "gpuSize": "NONE",
            "memory": 32,
            "storage": 250
        },
        "state": {
            "error": "",
            "isError": false,
            "stage": "READY",
            "status": "TERMINATED"
        },
        "status": "TERMINATED",
        "summary": "",
        "templateVersion": "v0.1.0",
        "updated": "2024-10-07T12:46:46.511Z",
        "userHandle": "google-apps|srinivas@deeporigin.com",
        "volumeDrns": [
            "..."
        ],
        "wasAutoStopped": false
    },
    "id": "...",
    "type": "ComputeBench"
}
```