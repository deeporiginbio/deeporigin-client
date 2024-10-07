This document describes how to create databases, columns in databases, and folders (workspaces) in the Deep Origin data hub.

## Folders

Folders can be created by specifying a name, and, optionally, a parent.

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
    api.create_workspace(name="test-folder")
    ```

    To create a folder within another folder, specify the parent:


    ```py
    api.create_workspace(
        name="test-folder-2",
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
    Currently, this package has limited support for creating database columns. We plan to expand the capabilities of this package.

To create a new database column in an existing database, run:

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

This code creates a new column in the existing database. To configure the type of the column, use the `type` argument. The type must a member of [DataType](../../ref/data-hub/types.md#src.utils.constants.DataType).
