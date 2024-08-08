This document describes how to create databases and folders in the Deep Origin Data Hub. 


## Create folders

Folders (or workspaces) can be created by specifying a name, and, optionally, a parent. 

=== "CLI"

    If no parent is specified, the folder will be created at the root level.

    ```bash
    deeporigin data new --name <name>  --folder
    ```

    To create a folder within another folder, specify the parent:


    ```bash
    deeporigin data new --name <name>  --folder --parent <parent-id>
    ```

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.create_workspace(name="test-workspace")
    ```

    To create a folder within another folder, specify the parent:


    ```py
    api.create_workspace(
        name="test-workspace-2",
        parent_id="parent-id",
    )
    ```

## Create databases

Databases can be created by specifying a name, and, optionally, a parent. 

=== "CLI"


    If no parent is specified, the database will be created at the root level.

    ```bash
    deeporigin data new --name <name>  --database
    ```

    To create a folder within another folder, specify the parent:


    ```bash
    deeporigin data new --name <name>  --database --parent <parent-id>
    ```

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.create_database(name="test-database")
    ```

    To create a folder within another folder, specify the parent:


    ```py
    api.create_database(
        name="test-database-2",
        parent_id="parent-id",
    )
    ```
