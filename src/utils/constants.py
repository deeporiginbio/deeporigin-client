"""this module contains constants used in the rest of this library"""

from beartype.typing import Literal

FileStatus = Literal["ready", "archived"]
"""Status of a file. Ready files are ready to be used, downloaded, and operated on."""


DATAFRAME_ATTRIBUTE_KEYS = {
    "metadata",
    "id",
    "reference_ids",
    "last_updated_row",
}


number = int | float

ENVS = Literal["edge", "prod", "staging"]

API_ENDPOINT = {
    "prod": "https://api.deeporigin.io",
    "staging": "https://api.staging.deeporigin.io",
    "edge": "https://api.edge.deeporigin.io",
}


ENV_VARIABLES = {
    "access_token": "DEEPORIGIN_TOKEN",
    "refresh_token": "DEEPORIGIN_REFRESH_TOKEN",
    "org_key": "DEEPORIGIN_ORG_KEY",
    "env": "DEEPORIGIN_ENV",
}

POCKETS_BASE_DIR = "~/.deeporigin/pockets"
"""Base directory for storing pocket files."""

TERMINAL_STATES = ["Failed", "Succeeded", "Cancelled"]
"""Terminal states for a job."""


# Supported atom symbols for small-molecule ligands in this toolkit.
# These intentionally exclude certain elements (e.g., metals, boron) that are
# not supported by downstream tools in typical docking workflows.
SUPPORTED_ATOM_SYMBOLS = {
    "H",
    "C",
    "N",
    "O",
    "F",
    "P",
    "S",
    "Cl",
    "Br",
    "I",
}
"""Set of supported atom symbols for ligands.

This set is used to validate ligands before docking and related operations.
Ligands containing atoms outside this set should be rejected by preparation
utilities.
"""
