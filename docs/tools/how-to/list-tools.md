This document describes how to discover and list tools on the Deep Origin platform. 

First, available tools are listed on the panel on left, with documentation on each tool and how to run them. 

## List tools using the API

All tools on the Deep Origin platform can by listed using the API.

```python
from deeporigin.platform import tools
all_tools = tools.get_all_tools()
```