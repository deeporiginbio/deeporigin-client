# Changelog

## `deeporigin v3.14.0`

## `deeporigin v3.1.0`

!!! warning "Deleted functions"
    The `data_hub.api.download_file`  function has been deleted. Use `data_hub.api.download_files` instead.

### Bugfixes

- Fixed a bug where documentation led pixi installs to use conda packages, not pypi packages
- Fixed a bug where the package could be installed on old and unsupported versions of python
- Fixed a bug where the default config file was always used
- Fixed a bug where expired token wasn't being refreshed correctly
- Fixed a bug where `download_files` wasn't working correctly

### New features

- Ability to run Deep Origin Tools
- Ability to fetch only some rows of a database as a dataframe
- Added modules for interacting with the platform API

### Miscellaneous

- Dropped support for Python 3.9, added support for Python 3.13
- Switched auth server to login.deeporgin.io


## `deeporigin v3.0.0`

!!! warning "Breaking changes"
    This release is a breaking change. Responses to all functions are no longer objects, but are now dictionaries wrapped in python [Boxes](https://pypi.org/project/python-box/). 

### Bugfixes

- Fixed a bug where writing more than 1000 rows failed
- Fixed a bug where constructing a dataframe from an empty DB failed


### New features

- Ability to add new data to existing dataframes
- Support for Deep Origin platform SDK
- Deep Origin dataframes now disallow adding rows
- Better support for writing data back to databases

### Miscellaneous 

- Initial support for integration tests with Deep Origin Data Hub
- Large refactor of testing

## `deeporigin v2.4.0`

### Bugfixes

- fixed a bug in how tokens were refreshed
- fixed a bug in CLI config use


### New features

- added the ability to download files in parallel
- added ability to upload files with file hashes 
- added the ability to decode user tokens
- support for a async client for the Deep Origin Data Hub API
- reference IDs now tracked in dataframes 
- dataframes now support expressions and users 
- dataframe slicing now restricts addition of new rows to slice
- dataframe updates more atomic: a single cell change only triggers updating that cell
- ability to delete terminated workstations 
- package wide ability to check for latest versions on PyPI

### Miscellaneous

- testing now includes python 3.13
- Data Hub API route now set to `datahub-api/api/`
- improved test coverage

## `deeporigin v2.2.1`

- added the ability to interact with the platform, including
    - fetching user secrets from the platform
    - fetching user info from the platform
    - fetching workstation info from the platform
- removed support for graphQL endpoints for secrets and variables
- support for new (breaking) syntax for deleting databases, folders, rows, and columns
- better support for caching tokens in filesystem
- numerous bug fixes and small improvements
- dropped support for Python 3.9
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