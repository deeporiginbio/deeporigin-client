
# Upload files

This page describes how to upload files to the Deep Origin
data hub, and assign them to individual cells.

## Upload a file to the data hub

To upload a file to the Deep Origin data hub, run the following commands:

=== "CLI"

    ```bash
    deeporigin data upload /path/to/test.fasta
    ```

    This will upload the file to your data hub, but the
    file will not yet be assigned to any database or cell. An example
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

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.upload_file("/path/to/file.fasta")
    ```

    This will upload the file to your data hub. An example is shown below:

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

## Upload a file to an existing row

To upload a file and assign it to a cell in an existing row, run the following commands:

=== "CLI"

    ```bash
    deeporigin data upload /path/to/test.fasta \
        --column <column_id> \
        --database <database_id> \
        --row <row_id>
    ```

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.upload_file("/path/to/file.fasta")
    ```

    To assign an uploaded file to a cell in an existing row, run the `assign_files_to_cell` function:

    ```py
        api.assign_files_to_cell(
        file_ids=["_file:6Hdhyc3t8xZ6pmyCrQy1t"],
        database_id="db-dna",
        column_id="base_sequence_file",
        row_id="row-id",
    )

    ```

## Upload a file to a new row of a database

To upload a file and assign it to a column in a new row, run the following commands:

=== "CLI"

    ```bash
    deeporigin data upload /path/to/test.fasta \
        --column <column_id> \
        --database <database_id> 
    ```
    
    This omits the `--row` parameter, which
    would create a new row in the database.

=== "Python"

    First, upload a file to your data hub by running:

    ```py
    from deeporigin.data_hub import api
    api.upload_file("/path/to/file.fasta")
    ```

    Second, assign the file to a new row by running:

    ```py
    api.assign_files_to_cell(
        file_ids=["_file:6Hdhyc3t8xZ6pmyCrQy1t"],
        database_id="db-dna",
        column_id="base_sequence_file",
    )

    ```
