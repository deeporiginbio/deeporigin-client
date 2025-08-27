# Using platform APIs using Deep Origin Platform Client

This document describes how to use the Deep Origin Platform Client. 

## Background

The typical way an end-user would use the Deep Origin python package would be to simply call functions. These functions call various APIs on the Deep Origin platform, using tokens and config information that is read from disk. This approach offers convenience for users who are taking actions as themselves on the platform, within a single organization.  

## Multi-user, multi-org

To make actions in multiple organizations, or as multiple users, a `client` can be passed to every function. 

First, construct a client using:


```python
from deeporigin.platform import Client

client = Client(token="my-secret-token", org_key="my-org")
```

Then, pass the client to any function. For example, to list tools:

```{.python notest}
from deeporigin.platform import tools_api

tools = tools_api.get_all_tools(client=client)
```