This document describes how to create databases, columns in databases, and folders in the Deep Origin Data Hub. 


## Folders

Folders (or workspaces) can be created by specifying a name, and, optionally, a parent. 

=== "CLI"

    If no parent is specified, the folder will be created at the root level.

    ```bash
    deeporigin data new folder --name <name>  
    ```

    To create a folder within another folder, specify the parent:


    ```bash
    deeporigin data new folder \
        --name <name> \ 
        --parent <parent-id>
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

## Databases

Databases can be created by specifying a name, and, optionally, a parent. 

=== "CLI"


    If no parent is specified, the database will be created at the root level.

    ```bash
    deeporigin data new --name <name>  --database
    ```

    To create a folder within another folder, specify the parent:


    ```bash
    deeporigin data new database \
        --name <name> \
        --parent <parent-id>
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


## Database columns

!!! warning "Work in progress"
    There is limited support for creating database columns from the python client and CLI at this time. Not all features are supported yet. 


Create a new database column in an existing database using:

=== "CLI"


    ```bash
    deeporigin data new column \
        --name <name> \
        --database <database-id> \
        --type <type>
    ```


=== "Python"


    ```py
    from deeporigin.data_hub import api
    api.add_database_column(
        database_id="existing-database-id",
        key="unique-key",
        type="integer",
        name="unique-name",
    )
    ```


This code creates a new column in that database. To configure the type of the column, use the `type` argument. The type must be one of [DataType](../../ref/data-hub/types.md#src.utils.DataType). 
