# View data

!!! danger "Deprecated"
    The Deep Origin DataHub is deprecated. 


This page describes how to describe and show the details of data objects in Deep Origin. To simply list objects, see [List data](./list-data.md).

## Describe data

### Describe rows

Describing rows provides metadata about the row, such as its ID, parent, and status.

!!! info "Describe vs. Show"
    This does not show you information in that row. To see data contained in that row, use the [`show`](#show-data) command.

To describe a row in a database in Deep Origin, run:

=== "CLI"

    ```bash
    deeporigin data describe <row-id>
    ```

    This will display a table similar to:

    ```
    ╭──────────────────┬────────────────────────────────────╮
    │ Column           │ Value                              │
    ├──────────────────┼────────────────────────────────────┤
    │ id               │ _row:WORR9xeGvG6mSg0yyDRlk         │
    │ parentId         │ _database:kyEws11L4KagGAaqRwONv    │
    │ type             │ row                                │
    │ dateCreated      │ 2024-05-08 17:59:32.512306         │
    │ dateUpdated      │ 2024-05-08 18:08:13.103            │
    │ createdByUserDrn │ scientist@deeporigin.com           │
    │ editedByUserDrn  │ scientist@deeporigin.com           │
    │ hid              │ data-2                             │
    │ validationStatus │ valid                              │
    ╰──────────────────┴────────────────────────────────────╯
    ```

    ??? tip "JSON output with `--json`"
        [JSON](https://www.json.org/) output can be requested by adding `--json`, and allows
        you to pipe out to a JSON processor like [jq](https://jqlang.github.io/jq/):

        ```bash
        deeporigin data describe <row-id> --json | jq
        ```

        ```json
        {
          "id": "_row:WORR9xeGvG6mSg0yyDRlk",
          "parentId": "_database:kyEws11L4KagGAaqRwONv",
          "type": "row",
          "dateCreated": "2024-05-08 17:59:32.512306",
          "dateUpdated": "2024-05-08 18:08:13.103",
          "createdByUserDrn": "scientist@deeporigin.com",
          "editedByUserDrn": "scientist@deeporigin.com",
          "hid": "data-2",
          "validationStatus": "valid"
        }
        ```

=== "Python"

    ```python
    from deeporigin.data_hub import api
    api.describe_row("_row:WORR9xeGvG6mSg0yyDRlk")
    ```

### Describe files

To describe a file in a database in Deep Origin, run:

=== "CLI"

    ```bash
    deeporigin data describe <file-id>
    ```

    This will display a table similar to:

    ```
    ╭──────────────────┬───────────────────────────────────╮
    │ Property         │ Value                             │
    ├──────────────────┼───────────────────────────────────┤
    │ id               │ _file:gBAK9tzFC5Cegx4NmSETc       │
    │ uri              │ s3://_file:gBAK9tzFC5Cegx4NmSETc  │
    │ name             │ db-dna.csv                        │
    │ status           │ ready                             │
    │ contentLength    │ 234                               │
    │ contentType      │ text/csv                          │
    │ dateCreated      │ 2024-05-08 01:01:48.925           │
    │ dateUpdated      │ 2024-05-08 01:01:48.925           │
    │ createdByUserDrn │ scientist@deeporigin.com          │
    ╰──────────────────┴───────────────────────────────────╯
    ```

    ??? tip "JSON output with `--json`"
        [JSON](https://www.json.org/) output can be requested by adding `--json`, and allows
        you to pipe out to a JSON processor like [jq](https://jqlang.github.io/jq/):

        ```bash
        deeporigin data describe <file-id> --json | jq
        ```

        ```json
        {
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
        ```

=== "Python"

    ```python
    from deeporigin.data_hub import api
    api.describe_file("file-id")
    ```

## Show data

### Show rows

To show the data in a row in a database in Deep Origin, run:

=== "CLI"

    ```bash
    deeporigin data show <row-id>
    ```

    This will display a table similar to:

    ```
    ╭───────────┬─────────────────────────────╮
    │ Column    │ Value                       │
    ├───────────┼─────────────────────────────┤
    │ File      │ _file:hnU7F62xeW8j0l1kR7YP1 │
    │ Float Num │ 112                         │
    │ selctcol  │ sdsd                        │
    ╰───────────┴─────────────────────────────╯
    ```

    ??? tip "JSON output with `--json`"
        [JSON](https://www.json.org/) output can be requested by adding `--json`, and allows
        you to pipe out to a JSON processor like [jq](https://jqlang.github.io/jq/):

        ```bash
        deeporigin data describe <file-id> --json | jq
        ```

        ```json
        {
          "File": "_file:hnU7F62xeW8j0l1kR7YP1",
          "Float Num": 112,
          "selctcol": "sdsd"
        }
        ```

=== "Python"

    ```python
    from deeporigin.data_hub import  api
    api.get_row_data("row-id")
    ```

    The data will be returned as a dictionary, where the keys are the column names and the values are values of the cells.

### Show databases

=== "CLI"

    ```bash
    deeporigin data show <database-id>
    ```

    This will display a table similar to:

    ```
    ╭───────────┬───────────┬──────────────────┬────────────┬───────────────────╮
    │  Status   │ stag-id   │ Customer Name    │ Status     │  Output Files     │
    ├───────────┼───────────┼──────────────────┼────────────┼───────────────────┤
    │ valid     │ stag-1    │ Blue Sun Corp    │ Processing │                   │
    │ valid     │ stag-2    │ Veridian Dynamics│ Complete   │   report.tar.gz   │
    ╰───────────┴───────────┴──────────────────┴────────────┴───────────────────╯
    ```

    ??? tip "JSON output with `--json`"
        [JSON](https://www.json.org/) output can be requested by adding `--json`, and allows
        you to pipe out to a JSON processor like [jq](https://jqlang.github.io/jq/):

        ```bash
        deeporigin data describe <database-id> --json | jq
        ```

        ```json
        {
          "Status": [
            "valid",
            "valid"
          ],
          "stag-id": [
            "stag-1",
            "stag-2"
          ],
          "Customer Name": [
            "Blue Sun Corp ",
            "Veridian Dynamics"
          ],
          "Status": [
            "Processing",
            "Complete"
          ],
          "Output Files": [
            null,
            "report.tar.gz"
          ]
        }
        ```

=== "Python"

    ```python
    from deeporigin.data_hub import api
    api.get_dataframe("database-id")
    ```

    The data will be returned as a [Pandas DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html).
