"""this module handles authentication actions and interactions
with tokens"""

import functools
import json
import os
import time
from urllib.parse import urljoin

import requests
from beartype import beartype
from deeporigin.utils.core import _get_api_tokens_filepath, read_cached_tokens

from .config import get_value as get_config
from .exceptions import DeepOriginException

__all__ = [
    "get_tokens",
    "cache_tokens",
    "remove_cached_tokens",
    "authenticate",
    "refresh_tokens",
]


@beartype
def tokens_exist() -> bool:
    """Check if the cached API tokens exist"""

    return os.path.isfile(_get_api_tokens_filepath())


@functools.cache
def get_tokens(
    *,
    refresh: bool = True,
) -> dict:
    """Get access token for accessing the Deep Origin API

    Gets tokens to access Deep Origin API.


    If an access token exists in the ENV, then it is used before
    anything else. If not, then tokens file is
    checked for access tokens, and used if they exist.
    On first use (within a Python session), tokens are
    refreshed if refresh_tokens is `True`. On subsequent uses,
    tokens are not refreshed even if refresh_tokens is `True`.

    Args:
        refresh: whether to refresh the token (default: `True`)

    Returns:
        :obj:`tuple`: API access and refresh tokens
    """

    if "DEEP_ORIGIN_ACCESS_TOKEN" in os.environ:
        return dict(
            access=os.environ["DEEP_ORIGIN_ACCESS_TOKEN"],
        )

    if tokens_exist():
        # tokens exist on disk
        tokens = read_cached_tokens()

        if refresh:
            tokens["access"] = refresh_tokens(tokens["refresh"])
            cache_tokens(tokens)

    else:
        # no tokens on disk. have to sign into the platform to get tokens
        tokens = authenticate()

    return tokens


@beartype
def cache_tokens(tokens: dict) -> None:
    """Save access and refresh tokens to a local file, for example, to
    enable variables/secrets to be regularly pulled without the user
    needing to regularly re-login.

    Args:
        token: dictionary with access and refresh tokens
    """

    with open(_get_api_tokens_filepath(), "w") as file:
        json.dump(tokens, file)


def remove_cached_tokens():
    """Remove cached API tokens"""

    if os.path.isfile(_get_api_tokens_filepath()):
        os.remove(_get_api_tokens_filepath())


@beartype
def authenticate() -> dict:
    """Get an access token for use with the Deep Origin API.
    The tokens are also cached to file

    Returns:
        :obj:`tuple`: API access token, API refresh token
    """

    config = get_config()

    # Get a link for the user to sign into the Deep Origin platform
    endpoint = urljoin(config.auth_domain, config.auth_device_code_endpoint)
    body = {
        "client_id": config.auth_client_id,
        "scope": "offline_access",
        "audience": config.auth_audience,
    }

    response = requests.post(endpoint, json=body)
    response.raise_for_status()
    response_json = response.json()
    device_code = response_json["device_code"]
    user_code = response_json["user_code"]
    verification_url = response_json["verification_uri_complete"]
    sign_in_poll_interval_sec = response_json["interval"]

    # Prompt the user to sign into the Deep Origin platform
    print(
        (
            "To connect to the Deep Origin OS, "
            f"navigate your browser to \n\n{verification_url}\n\n"
            f'and verify the confirmation code is "{user_code}", '
            'and click the "Confirm" button.'
        )
    )

    # Wait for the user to sign into the Deep Origin platform
    endpoint = urljoin(config.auth_domain, config.auth_token_endpoint)
    body = {
        "grant_type": config.auth_grant_type,
        "device_code": device_code,
        "client_id": config.auth_client_id,
    }
    while True:
        response = requests.post(endpoint, json=body)
        if response.status_code == 200:
            break
        if (
            response.status_code != 403
            or response.json().get("error", None) != "authorization_pending"
        ):
            raise DeepOriginException(
                message="Sign in to Deep Origin failed. Please try again."
            )
        time.sleep(sign_in_poll_interval_sec)

    response_json = response.json()
    api_access_token = response_json["access_token"]
    api_refresh_token = response_json["refresh_token"]

    tokens = dict(access=api_access_token, refresh=api_refresh_token)

    cache_tokens(tokens)

    return tokens


def refresh_tokens(api_refresh_token: str) -> str:
    """Refresh the access token for the Deep Origin OS

    Args:
        api_refresh_token (:obj:`str`): API refresh token

    Returns:
        :obj:`str`: new API access token
    """
    config = get_config()

    endpoint = urljoin(config.auth_domain, config.auth_token_endpoint)
    body = {
        "grant_type": "refresh_token",
        "client_id": config.auth_client_id,
        "client_secret": config.auth_client_secret,
        "refresh_token": api_refresh_token,
    }
    response = requests.post(endpoint, json=body)
    response.raise_for_status()
    response_json = response.json()
    api_access_token = response_json["access_token"]

    return api_access_token
