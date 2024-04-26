


## Fetch databases 


=== "CLI"

    Databases can be viewed in the command line using:

    ```bash
    deeporigin data show-db <id-of-database>
    ```

    The CLI can also return JSON output that can be piped
    to a JSON parser:

    ```bash
    deeporigin data show-db <id-of-database> --json | jq
    ```

    This is useful to quickly view the data in the database. 
    To save this database to a CSV, use:

    ```bash
    deeporigin data cp deeporigin://<id-of-database> <local-folder>
    ```

    The `deeporigin://` prefix identifies a resource as
    existing on Deep Origin, and the `cp` command downloads 
    the database as a CSV and all included files are downloaded 
    to the `<local-folder>`.




=== "Python Client"
    
    Databases can be retrieved as [`pandas.DataFrame`](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html)s using:

    ```python
    from deeporigin.managed_data.api import get_dataframe
    df = get_dataframe("id-of-database")
    ```

    The dataframe contains all rows and columns of the 
    database as seen on Deep Origin. This dataframe can be saved
    to disk, or used for further downstream analysis.

    By default, files in databases are referred by file names. 



## Downloading files