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

To delete a row in a database, run:

=== "CLI"

    ```bash
    deeporigin data delete \
        --row <row-id> \
        --database <database-id>

    ```

=== "Python"



    ```py
    from deeporigin.data_hub import api
    api.delete_rows(row_ids=["row-1", "row-2"], database_id="database-id")
    ```

    !!! tip "Deleting multiple rows"
        Note that the python client allows you to delete multiple rows at once.


## Delete databases

To delete an entire database, run:

=== "CLI"

    ```bash
    deeporigin data delete --database <db-id>
    ```

    !!! tip "Aliases"
        The following alias also works: `deeporigin data delete -d <db-id>`.

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.delete_database(database_id=<db-id>)
    ```


## Delete folders (workspaces)

To delete an entire folder, including all contained databases, run:

=== "CLI"

    ```bash
    deeporigin data delete --folder <id>
    ```

    !!! tip "Aliases"
        The following aliases also work `--workspace`, `-w`, `--ws`, `-f` instead of `--folder`.

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.delete_workspace(workspace_id=<folder-id>)
    ```



## Delete database column

To delete a column in a databas, run the following command, specifying the ID of the column to delete and its parent database:

=== "CLI"

    ```bash
    deeporigin data delete \
        --column <ids> \
        --database <database-id> \
    ```

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.delete_database_column(
        column_id="col-id",
        database_id="database-id",
    )
    ```
