This page describes how to write data to the Deep Origin data hub.

## Write data to a cell in a database

### Write text data

Run the following to write a text value to a cell in a database:

=== "CLI"

    ```bash
    deeporigin data write "some-string" \
        --database $database_id \
        --column $column_id \
        --row $row_id

    ```

=== "Python"

    ```py
    from deeporigin.managed_data import api
    api.set_cell_data(
        "some string",
        database_id="database ID or name",
        column_id="column ID or programmatic key",
        row_id="row ID",
    )
    ```

### Write numeric data

Run the following to write a numeric value to a cell in a database:

=== "CLI"

    ```bash
    deeporigin data write 123 \
        --database $database_id \
        --column $column_id \
        --row $row_id

    ```

=== "Python"

    ```python
    from deeporigin.managed_data import api
    api.set_cell_data(
        1,
        database_id="database ID or name",
        column_id="column ID or programmatic key",
        row_id="row ID",
    )
    ```

Numeric data will be coerced to the data type of the underlying cell.

### Write select data

This section describes how to write data to a cell in a database that is of type `select`. This includes cells where a single option can be selected from a list, and cells where more than one value can be selected from a list of options.

=== "CLI"

    ```bash
    deeporigin data write "option A" \
        --database $database_id \
        --column $column_id \
        --row $row_id

    ```

=== "Python"

    ```py
    from deeporigin.managed_data import api
    api.set_cell_data(
        "option A",
        database_id="database ID or name",
        column_id="column ID or programmatic key",
        row_id="row ID",
    )
    ```

The value must be one of the options in the list. If it is not, an error will be raised showing the list of legal options.

### Write Boolean data

Run the following to write a Boolean value to a cell in a database:

=== "CLI"

    ```bash
    deeporigin data write "option A" \
        --database $database_id \
        --column $column_id \
        --row $row_id

    ```

=== "Python"

    ```python
    from deeporigin.managed_data import api
    api.set_cell_data(
        True,
        database_id="database ID or name",
        column_id="column ID or programmatic key",
        row_id="row ID",
    )
    ```

    ??? tip "Unsetting Boolean data"
        To unset a cell, so that it contains no data, write `None` to the cell as follows: 

        ```py
        from deeporigin.managed_data import api
        api.set_cell_data(
            None,
            database_id="database ID or name",
            column_id="column ID or programmatic key",
            row_id="row ID",
        )
        ```
