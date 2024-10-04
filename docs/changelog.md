# Changelog


## `deeporigin v2.2.1`

- added the ability to interact with the platform, including
    - fetching user secrets from the platform
    - fetching user info from the platform
    - fetching workstation info from the platform
- removed support for graphQL endpoints for secrets and variables
- support for new (breaking) syntax for deleting databases, folders, rows, and columns
- better support for caching tokens in filesystem
- numerous bug fixes and small improvements
- better syntax for interacting with the `api.list_files` function
- Improvements to Deep Origin DataFrames, including:
    - ability to create new columns in dataframes and push them to a database
    - DataFrames now show names of user who created and last updated a database, not just a ID
    - DataFrame pretty printing methods now never interrupt core logic, and fall back to printing methods of superclass on failure
    - DataFrames now support `df.tail()` and `df.head()`

## `deeporigin v2.1.1` 


- fixed a bug where dataframe header links were sometimes broken
- improved display of dataframe headers
- improved typing in dataframes
- added ability to check for newest versions on PyPI
- version numbers are now [PEP 440](https://peps.python.org/pep-0440/) compliant

## `deeporigin v2.1.0` 


- dropped support for column keys
- added the [deeporigin.DataFrame][src.data_hub.dataframe.DataFrame] class, which:
    - is a drop-in replacement for pandas.DataFrame
    - supports all pandas methods
    - supports automatic syncing with a Deep Origin database