This document describes how to delete databases and folders in the Deep Origin Data Hub. 

## Delete folders or databases or rows

=== "CLI"

    Mulitiple rows, databases and folders can be deleted using:

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


!!! warning "Exercise caution"
    Deleting a workspace deletes all databases inside it. Deleting a database deletes all rows inside it. All resources will be deleted without asking for confirmation. 