This document describes how to delete 

- database rows
- database columns
- databases and 
- folders 

in the Deep Origin Data Hub. 

!!! warning "Exercise caution"
    - Deleting a workspace deletes all databases inside it. 
    - Deleting a column destroys all data in that column, including all files assigned to cells in that column. 
    - Deleting a database deletes all rows inside it. 

    All resources will be deleted without asking for confirmation. 

## Delete database rows, databases, and folders

=== "CLI"

    Multiple rows, databases and folders can be deleted using:

    ```bash
    deeporigin data delete --ids <ids>
    ```



=== "Python"



    ```py
    from deeporigin.data_hub import api
    api.delete_rows(row_ids=["row-1","database-1","workspace-1"])
    ```

    !!! Info "Rows?"
        Rows, workspaces and databases can all be deleted using the `delete_rows` method. 





## Delete database columns 

Columns in databases can be deleted using:

=== "CLI"

    

    ```bash
    deeporigin data delete --ids <ids> --columns
    ```

    !!! warning "Column IDs"
        Currently, columns can only be deleted using Column IDs, which are not the same as column names. To view the column IDs of a database, use 

        ```bash
        deeporigin data describe <database-id>
        ```



=== "Python"



    ```py
    from deeporigin.data_hub import api
    api.delete_database_column(column_id="col-id")
    ```

