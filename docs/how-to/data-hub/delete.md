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

## Delete database rows, databases, and folders

To delete multiple rows, databases and folders, run:

=== "CLI"

    ```bash
    deeporigin data delete --ids <ids>
    ```

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.delete_rows(row_ids=["row-1","database-1","folder-1"])
    ```

    !!! Info "Rows?"
        Rows, folders and databases can all be deleted using the `delete_rows` method. 

## Delete database columns

To delete columns in databases, run:

=== "CLI"

    ```bash
    deeporigin data delete --ids <ids> --columns
    ```

    !!! warning "Column IDs"
        Currently, columns can only be deleted using their IDs, which are distinct from their names. To view the IDs of the columns of a database, run: 

        ```bash
        deeporigin data describe <database-id>
        ```

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.delete_database_column(column_id="col-id")
    ```
