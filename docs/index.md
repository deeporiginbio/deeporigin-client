
# Deep Origin CLI and Python client

<div class="grid cards" markdown>

- :fontawesome-solid-handshake-simple: **Access your DO resources**  
Work with your DO data, variables, secrets, and other resources through a CLI or Python.
- :octicons-unlock-24: **Free and open-source**
Install onto your computer to use your data, variables, secrets, and other resources outside DO.
- :material-download: **Easy to install**
Just run `pip install deeporigin`.
- :fontawesome-brands-python: **Pure Python**
Lightweight, written in pure Python. Works on any system that can run Python.

</div>

## Examples

The Deep Origin CLI and Python client allow you to programmatically
interact with the [Deep Origin OS](https://os.deeporigin.io/).
The example below illustrates how to use the CLI and Python library to
retrieve the database row with ID `data-1`.

=== "Terminal"

    ```bash
    deeporigin data show data-1

    ╭──────────────┬─────────────────────────────────╮
    │ Column      │ Value                        │
    ├──────────────┼─────────────────────────────────┤
    │ Boolean     │ False                        │
    │ Float       │ 112                          │
    │ Select      │ sdsd                         │
    │ Date        │ 2024-06-19 00:00:00          │
    │ File        │ _file:hnU7F62xeW8j0l1kR7YP1  │
    ...
    ╰──────────────┴─────────────────────────────────╯
    ```

=== "Python"

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
