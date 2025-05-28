"""This module contains utility functions that are used internally by the Python Client and the CLI, reliant on the config module"""

from urllib.parse import urljoin

from beartype import beartype
from beartype.typing import TypeVar

from deeporigin.config import get_value
from deeporigin.utils.constants import ObjectType

T = TypeVar("T")


@beartype
def construct_resource_url(
    *,
    name: str,
    row_type: ObjectType,
) -> str:
    """Constructs the URL for a resource

    Args:
        name (str): name of the resource
        row_type (ObjectType): type of the resource

    Returns:
        str: URL for the resource
    """

    org = get_value()["organization_id"]

    return urljoin(_get_domain_name(), f"org/{org}/data/{row_type}/{name}")


@beartype
def _nucleus_url() -> str:
    """Returns URL for nucleus API endpoint"""
    url = urljoin(
        get_value()["api_endpoint"],
        get_value()["nucleus_api_route"],
    )
    if not url.endswith("/"):
        url += "/"

    return url


def _get_domain_name() -> str:
    """utility function to get domain name based on env"""

    env = get_value()["env"]
    if env == "prod":
        return "https://os.deeporigin.io"
    else:
        return f"https://os.{env}.deeporigin.io"
