
# Deep Origin CLI and Python Client

<div class="grid cards" markdown>

- :fontawesome-brands-python: **Pure python**       
Lightweight implementation, written in pure Python. Works on any system that can run Python.
- :fontawesome-solid-handshake-simple: **Talks to Deep Origin APIs**  
Talks to Deep Origin APIs, including APIs to work with managed data, and ComputeBench.
- :octicons-unlock-24: **Freely available**         
Freely available and open source. Install on your computer!
- :material-download: **Easy install**     
Just `pip install deeporigin`.

</div>

## Overview


The Deep Origin CLI and Python client allow you to programmatically
interact with the [Deep Origin platform](https://os.deeporigin.io/).

```python
from deeporigin.managed_data import api

api.get_row_data("data-1")

 # {'Status': 'Processing',
 #  'Age (years)': 15,
 #  'Gender': 'F',
 #  'Order date': '2024-03-01T00:00:00',
 #  'Received by client': '2023-04-05T00:00:00',
 #  'Sent to client': '2024-03-03T00:00:00',
 #  'Sent by client': '2024-03-05T00:00:00',
 #  'Raw reads': ['_file:UjKrB0QibhgBDAVvGQ3VP'],
 #  'Date completed': '2024-03-08T00:00:00',
 # }
```
