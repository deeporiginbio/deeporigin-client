This document describes how to delete objects in the Deep Origin data hub:

- Database rows
- Database columns
- Databases
- Folders (workspaces)

!!! danger "Exercise caution"
    - Deleting a folder deletes all of the databases in the folder.
    - Deleting a column destroys all of the data in that column, including all of the files assigned to the cells in that column.
    - Deleting a database deletes all of the rows in the database.

    All resources will be deleted without asking for confirmation. 

## Delete database rows

To delete rows in a database, run:

=== "CLI"

    ```bash
    deeporigin data delete --ids <ids>
    ```

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.delete_rows(row_ids=["row-1","database-1","folder-1"])
    ```


## Delete databases

To delete an entire database, run:

=== "CLI"

    ```bash
    deeporigin data delete --ids <ids>
    ```

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.delete_database(database_id=<row-ID>)
    ```


## Delete folders (workspaces)

To delete an entire folder, including all contained databases, run:

=== "CLI"

    ```bash
    deeporigin data delete --ids <ids>
    ```

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.delete_workspace(workspace_id=<folder-ID>)
    ```



## Delete database columns

To delete columns in databases, run the following command, specifying the IDs of the columns to delete and their parent database:

=== "CLI"

    ```bash
    deeporigin data delete \
        --ids <ids> \
        --database <database-id> \
        --columns
    ```

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.delete_database_column(
        column_id="col-id",
        database_id="database-id",
    )
    ```
