## List rows 


!!! Note "Rows"
    In this context, a row can mean a database, a workspace, or a row.

=== "CLI"

    To list all workspaces, databases and rows on Deep Origin:

    ```bash
    deeporigin data ls
    ```

=== "Python Client"
    
    The python client offers a number of ways to list resources 
    on Deep Origin.

    First, we start off by importing the necessary modules:


    ```python
    from deeporigin.managed_data import _api, api
    ```

    We can list rows using:

    ```python
    _api.list_rows()
    ```

    To find only databases, we add a [`RowType`][src.managed_data.schema.RowType] argument:

    ```python
    _api.list_rows(row_type="databases")
    ```

## List files 

=== "CLI"

    To list files on Deep Origin:

    ```bash
    deeporigin data list-files 
    ```

    This shows a table with file_IDs, names, status, and other
    information about each file. 

    To list unassigned files, pass the `--unassigned` flag:

    ```bash
    deeporigin data list-files --unassigned
    ```

    The CLI can also return JSON output that can be piped
    to a JSON parser:

    ```bash
    deeporigin data list-files --unassigned --json | jq
    ```


=== "Python Client"

    First, we start off by importing the necessary modules:


    ```python
    from deeporigin.managed_data._api import list_files
    ```

    We can list all files on Deep Origin using:

    ```python
    list_files()
    ```

    To find only unassigned files, we can use:

    ```python
    list_files(is_unassigned=True)
    ```

    To find files that are assigned to a specific row:

    ```python
    list_files(assigned_row_ids=["row-1"])
    ```