"""module to interact with the platform API"""

from urllib.parse import urljoin

import diskcache as dc
import requests
from beartype import beartype
from deeporigin import auth
from deeporigin.config import get_value
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils.core import _ensure_do_folder
from deeporigin.utils.network import _get_domain_name


def _make_get_request(endpoint: str) -> dict:
    """make a request to the platform API"""

    if not endpoint.startswith("/"):
        endpoint = "/" + endpoint

    tokens = auth.get_tokens(refresh=False)
    access_token = tokens["access"]

    headers = {
        "accept": "application/json",
        "authorization": f"Bearer {access_token}",
        "cache-control": "no-cache",
    }

    url = urljoin(_get_domain_name(), f"api/{endpoint}")

    response = requests.get(
        url,
        headers=headers,
    )

    if response.status_code != 200:
        raise DeepOriginException(message=response.text, fix=url)

    data = response.json()
    if "data" in data.keys():
        return data["data"]

    else:
        return data


@beartype
def _get_org_id() -> str:
    value = get_value()
    return value["organization_id"]


def resolve_user(user_id: str):
    """get details about a user in the platform"""

    endpoint = f"/organizations/{_get_org_id()}/users/{user_id}"

    return _make_get_request(endpoint)


def whoami() -> dict:
    """get details about currently signed in user"""

    return _make_get_request("/users/me")


def get_workstations() -> list[dict]:
    """get information about all workstations in the organization"""
    return _make_get_request(f"/computebenches/{_get_org_id()}")


@beartype
def get_user_name(user_id: str) -> str:
    """get user name from user ID"""

    CACHE_PATH = _ensure_do_folder() / "user_ids"

    cache = dc.Cache(CACHE_PATH)

    if cache.get(user_id) is not None:
        return cache.get(user_id)

    response = resolve_user(user_id)

    name = response["data"]["attributes"]["name"]
    cache.set(user_id, name)
    return name


def get_last_edited_user_name(row):
    """given a row object, return a human readable user name
    of the user who last edited the row"""

    if row.edited_by_user_drn is None:
        # fall back to created by user
        user_id = row.created_by_user_drn.split(":")[-1]
    else:
        user_id = row.edited_by_user_drn.split(":")[-1]
    return get_user_name(user_id)


@beartype
def get_variables_and_secrets() -> list[dict]:
    """get variables and secrets for user and org"""

    # we don't need to provide a bench ID, so pass a placeholder
    response = _make_get_request(f"/computebenches/{_get_org_id()}/placeholder/secrets")

    return response["data"]
