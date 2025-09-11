"""This module contains utility functions that are used internally by the Python Client and the CLI, reliant on the config module"""

from beartype.typing import TypeVar

from deeporigin.config import get_value

T = TypeVar("T")


def _get_domain_name() -> str:
    """utility function to get domain name based on env"""

    env = get_value()["env"]
    if env == "prod":
        return "https://os.deeporigin.io"
    else:
        return f"https://os.{env}.deeporigin.io"
