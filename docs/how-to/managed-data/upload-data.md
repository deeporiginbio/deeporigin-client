
This page describes how to upload data to Deep Origin Data,
and assign it to a cell in a row.

## Upload and Assign Files


=== "CLI"

    To upload a file to Deep Origin Data, use the `upload` command:

    ```bash
    deeporigin data upload /path/to/test.fasta 
    ```

    This will upload a file to the "Staging Area", and this
    file is not assigned to any database or cell. An example
    response is shown below:

    ```
    ╭──────────────────┬──────────────────────────────╮
    │ Property         │ Value                        │
    ├──────────────────┼──────────────────────────────┤
    │ name             │ test.fasta                   │
    │ contentType      │ data/fasta                   │
    │ contentLength    │ 554588                       │
    │ id               │ _file:36ufKT2Sej22coSEOizpK  │
    │ status           │ ready                        │
    │ uri              │ s3://foo                     │
    │ dateCreated      │ 2024-06-18T14:48:33.501Z     │
    │ dateUpdated      │ 2024-06-18T14:48:33.501Z     │
    │ createdByUserDrn │ levins@deeporigin.com        │
    ╰──────────────────┴──────────────────────────────╯


    ```

    ??? tip "JSON output with `--json`"
        [JSON](https://www.json.org/) output can be requested by adding `--json`, and allows
        you to pipe out to a JSON processor like [jq](https://jqlang.github.io/jq/):

        ```bash
        deeporigin upload /path/to/test.fasta | jq
        ```

        ```json
        {
          "name": "test.fasta",
          "contentType": "data/fasta",
          "contentLength": 554588,
          "id": "_file:36ufKT2Sej22coSEOizpK",
          "status": "ready",
          "uri": "s3://foo",
          "dateCreated": "2024-06-18T14:51:43.876Z",
          "dateUpdated": "2024-06-18T14:51:43.876Z",
          "createdByUserDrn": "levins@deeporigin.com"
        }

        ```

    To upload a file and assign it to an existing cell, use:


    ```bash
    deeporigin data upload /path/to/test.fasta \
        --column <column_id> \
        --database <database_id> \
        --row <row_id>
    ```

    To upload a file and assign it to a column in a new row,
    use: 


    ```bash
    deeporigin data upload /path/to/test.fasta \
        --column <column_id> \
        --database <database+_id> 
    ```

 

=== "Python client"
    

    We can upload a file using 

    ```py
    from deeporigin.managed_data import api
    api.upload_file("/path/to/file.fasta")
    ```

    This uploads files to the "Staging Area" of Deep Origin
    Data. An example is shown below:

    ```json
    {
        "name": "file.fasta",
        "contentType": "data/foo",
        "contentLength": 55454688,
        "id": "_file:6Hdhyc3t8xZ6pmyCrQy1t",
        "status": "ready",
        "uri": "s3://data.<org-name>/_file:6Hdhyc3t8xZ6pmyCrQy1t",
        "dateCreated": "2024-06-18T14:18:37.409Z",
        "dateUpdated": "2024-06-18T14:18:37.409Z",
        "createdByUserDrn": "haldane@deeporigin.com",
    }
    ```


    To assign an uploaded file to a cell in an existing row, use the `assign_files_to_cell` function as follows:

    ```py
        api.assign_files_to_cell(
        file_ids=["_file:6Hdhyc3t8xZ6pmyCrQy1t"],
        database_id="db-dna",
        column_id="base_sequence_file",
        row_id="row-id",
    )

    ```

    To assign an uploaded file to a cell on a new row, use the `assign_files_to_cell` function as follows:

    ```py
        api.assign_files_to_cell(
        file_ids=["_file:6Hdhyc3t8xZ6pmyCrQy1t"],
        database_id="db-dna",
        column_id="base_sequence_file",
    )

    ```


