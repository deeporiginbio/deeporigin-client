This document describes how to run tools on the Deep Origin platform. 

## Running packaged tools

In general, tools are best run by calling functions in the `tools.run` module. For example, to run [Autodock Vina](../../tools/vina.md), use:

```python
from deeporigin.tools import run
run.autodock_vina(...)
```

## Running any tool

Any tool in the Deep Origin tool registry can be run using a generic tool run function. 

```python
from deeporigin.tools.utils import run_tool
run_tool(data)
```

where `data` is a dictionary representing the DTO (Data Transfer Object) for the tool.

The DTO requires the following keys:


### `toolId`

The toolId is used to identify the tool in the Deep Origin tool registry. For AutoDock Vina, this is `deeporigin/autodock-vina`.


### `inputs`

The `inputs` key is a dictionary of input parameters for the tool. 

### `outputs`

The `outputs` key is a dictionary of output parameters for the tool.


### `clusterId` (Optional)

The `clusterId` (optional) key is the id of the cluster to run the tool on. If this is not provided, it will be run on the default (us-west-2) cluster. 