# Changelog

### `deeporigin v2.1.1` 


- fixed a bug where dataframe header links were sometimes broken
- improved display of dataframe headers
- improved typing in dataframes
- added ability to check for newest versions on PyPI
- version numbers are now [PEP 440](https://peps.python.org/pep-0440/) compliant

### `deeporigin v2.1.0` 


- dropped support for column keys
- added the [deeporigin.DataFrame][src.data_hub.dataframe.DataFrame] class, which:
    - is a drop-in replacement for pandas.DataFrame
    - supports all pandas methods
    - supports automatic syncing with a Deep Origin database