!!! danger "Deprecated"
    The Deep Origin DataHub is deprecated. 


This page describes how to list resources on Deep Origin.
Listing folders (workspaces), databases, rows and files can be used
to discover the resources available on Deep Origin, and show
their IDs for further queries.

## List folders

To list all of your Deep Origin folders, run:

=== "CLI"

    ```bash
    deeporigin data list --folders
    ```

    This will display a screen similar to below:

    ```
    ╭─────────────────┬────────────┬───────────────╮
    │ Name            │ Type       │ ID            │
    ├─────────────────┼────────────┼───────────────┤
    │ Secret Project  │ workspace  │ secret        │
    │ QC Efforts      │ workspace  │ qc-efforts    │
    │ Covid Target    │ workspace  │ corona        │
    ╰─────────────────┴────────────┴───────────────╯
    ```

    ??? tip "JSON output with `--json`"
        [JSON](https://www.json.org/) output can be requested by adding `--json`, and allows
        you to pipe out to a JSON processor like [jq](https://jqlang.github.io/jq/):

        ```bash
        deeporigin data list --folders --json | jq
        ```

        ```json
        [
          {
            "id": "_workspace:ncZhbnYXXfg0zWNcGKTnz",
            "parentId": null,
            "hid": "secret",
            "type": "workspace",
            "name": "Secret Project"
          },
          {
            "id": "_workspace:sDTZKGZXOkhGw6XSg2Jla",
            "parentId": null,
            "hid": "qc-efforts",
            "type": "workspace",
            "name": "QC Efforts"
          },
          {
            "id": "_workspace:UWhlb4Wyzh2R7bySapY2m",
            "parentId": null,
            "hid": "corona",
            "type": "workspace",
            "name": "Covid Target"
          }
        ]
        ```

    ??? tip "Combining outputs"
        You can combine multiple output types. For example, to list all folders and databases:

        ```bash
        deeporigin data list --folders --databases
        ```

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.list_rows(row_type="workspace")
    ```

## List databases

To list all of your databases in Deep Origin, run:

=== "CLI"    

    ```bash
    deeporigin data list --databases
    ```

    This will display a screen similar to below:

    ```
    ╭────────────────┬───────────┬────────────────╮
    │ Name           │ Type      │ ID             │
    ├────────────────┼───────────┼────────────────┤
    │ First Database │ database  │ first          │
    │ QC Efforts     │ database  │ qc-efforts     │
    │ Covid Target   │ database  │ corona         │
    ╰────────────────┴───────────┴────────────────╯
    ```

    ??? tip "JSON output with `--json`"

        [JSON](https://www.json.org/) output can be requested by adding `--json`, and allows
        you to pipe out to a JSON processor like [jq](https://jqlang.github.io/jq/):

        ```bash
        deeporigin data list --databases --json | jq
        ```

        ```json
        [
          {
            "id": "_database:ncZhbnYXXfg0zWNcGKTnz",
            "parentId": "_workspace:ncZhbnYXXfg0zWNcGKTnz",
            "hid": "db-first",
            "type": "database",
            "name": "First Database"
          },
          {
            "id": "_database:sDTZKGZXOkhGw6XSg2Jla",
            "parentId": "_workspace:sDTZKGZXOkhGw6XSg2Jla",
            "hid": "db-qc",
            "type": "database",
            "name": "QC Efforts"
          },
          {
            "id": "_database:UWhlb4Wyzh2R7bySapY2m",
            "parentId": "_workspace:UWhlb4Wyzh2R7bySapY2m",
            "hid": "db-covid",
            "type": "database",
            "name": "Covid Target"
          }
        ]
        ```

    ??? tip "Combining outputs"
        You can combine multiple output types. For example, to list all folders and databases:

        ```bash
        deeporigin data list --folders --databases
        ```

=== "Python"

    ```python
    from deeporigin.data_hub import api
    api.list_rows(row_type="database")
    ```

## List rows

To list all of your database rows in Deep Origin:

=== "CLI"

    ```bash
    deeporigin data list --rows
    ```

    This will display a screen similar to below:

    ```
    ╭────────┬────────┬────────╮
    │ Name   │ Type   │ ID     │
    ├────────┼────────┼────────┤
    │        │ row    │ data-1 │
    │        │ row    │ data-2 │
    │        │ row    │ data-3 │
    │        │ row    │ data-4 │
    │        │ row    │ data-5 │
    ╰────────┴────────┴────────╯
    ```

=== "Python"

    ```python
    from deeporigin.data_hub import api
    api.list_rows(row_type="row")
    ```

## List files

To list all of your files in Deep Origin, run:

=== "CLI"

    ```bash
    deeporigin data list --files
    ```

    This will display a screen similar to below:

    ```
    ╭────────────┬──────────┬─────────────────────────────╮
    │ Name       │ Status   │ ID                          │
    ├────────────┼──────────┼─────────────────────────────┤
    │ db-dna.csv │ ready    │ _file:gBAK9tzFC5Cegx4NmSETc │
    │ seqs.gz    │ ready    │ _file:FgVjcv8zzyPho6FME8QFp │
    │ db-rna.csv │ ready    │ _file:hnU7F62xeW8j0l1kR7YP1 │
    ╰────────────┴──────────┴─────────────────────────────╯
    ```


    ??? tip "JSON output with `--json`"
        [JSON](https://www.json.org/) output can be requested by adding `--json`, and allows
        you to pipe out to a JSON processor like [jq](https://jqlang.github.io/jq/):

        ```bash
        deeporigin data list --databases --json | jq
        ```

        ```json
        [
          {
            "file": {
              "id": "_file:gBAK9tzFC5Cegx4NmSETc",
              "uri": "s3://_file:gBAK9tzFC5Cegx4NmSETc",
              "name": "db-dna.csv",
              "status": "ready",
              "contentLength": 234,
              "contentType": "text/csv",
              "dateCreated": "2024-05-08 01:01:48.925",
              "dateUpdated": "2024-05-08 01:01:48.925",
              "createdByUserDrn": "scientist@deeporigin.com"
            }
          },
          {
            "file": {
              "id": "_file:FgVjcv8zzyPho6FME8QFp",
              "uri": "s3://_file:FgVjcv8zzyPho6FME8QFp",
              "name": "seqs.gz",
              "status": "ready",
              "contentLength": 554588,
              "contentType": "zip/gz",
              "dateCreated": "2024-05-08 18:08:09.149",
              "dateUpdated": "2024-05-08 18:08:09.149",
              "createdByUserDrn": "scientist@deeporigin.com"
            },
            "assignments": [
              {
                "rowId": "_row:WORR9xeGvG6mSg0yyDRlk"
              }
            ]
          },
          {
            "file": {
              "id": "_file:hnU7F62xeW8j0l1kR7YP1",
              "uri": "s3://_file:hnU7F62xeW8j0l1kR7YP1",
              "name": "db-rna.csv",
              "status": "ready",
              "contentLength": 234,
              "contentType": "text/csv",
              "dateCreated": "2024-05-08 18:07:57.655",
              "dateUpdated": "2024-05-08 18:07:57.655",
              "createdByUserDrn": "scientist@deeporigin.com"
            },
            "assignments": [
              {
                "rowId": "_row:mlNnmNkfktz7GT5qpjyrF"
              }
            ]
          }
        ]
        ```

    !!! warning "Listing files cannot list other objects"
        If you pass `--files` to the list command, all other 
        arguments are ignored. As a result, 

        ```bash
        deeporigin data list --files --databases
        ```
        will only list files.

=== "Python"

    First, we start off by importing the necessary modules:

    We can list all files on Deep Origin using:

    ```python
    from deeporigin.data_hub import api
    api.list_files()
    ```

    To find only unassigned files, we can use:

    ```python
    api.list_files(is_unassigned=True)
    ```

    To find files that are assigned to a specific row:

    ```python
    api.list_files(assigned_row_ids=["row-1"])
    ```
