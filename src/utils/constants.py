"""this module contains constants used in the rest of this library"""

from dataclasses import dataclass

from beartype.typing import Literal

ObjectType = Literal["row", "database", "workspace"]
"""Type of a row. In Deep Origin, a row can be a database row, a database or a workspace"""

FileStatus = Literal["ready", "archived"]
"""Status of a file. Ready files are ready to be used, downloaded, and operated on."""

DataType = Literal[
    "integer",
    "select",
    "date",
    "text",
    "file",
    "reference",
    "editor",
    "float",
    "boolean",
]
"""Type of a column in a Deep Origin database. See [this page in the documentation](https://docs.deeporigin.io/docs/os/data-hub/databases/columns) for more information."""

DATAFRAME_ATTRIBUTE_KEYS = {
    "metadata",
    "id",
    "reference_ids",
    "last_updated_row",
}


Cardinality = Literal["one", "many"]
"""The cardinality defines whether a cell in a database can contain one or multiple objects"""

IDFormat = Literal["human-id", "system-id"]
"""Format of an ID"""

DatabaseReturnType = Literal["dataframe", "dict"]
"""Return type for [api.get_dataframe][src.data_hub.api.get_dataframe]"""


@dataclass
class PREFIXES:
    """Prefixes for CLI and Python client"""

    DO = "do://"
    FILE = "_file"
    DB = "_database"
    ROW = "_row"
    FOLDER = "_workspace"


number = int | float

ENVS = Literal["edge", "prod", "staging"]

API_ENDPOINT = {
    "prod": "https://api.deeporigin.io",
    "staging": "https://api.staging.deeporigin.io",
    "edge": "https://api.edge.deeporigin.io",
}
