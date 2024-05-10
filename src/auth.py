import functools
import json
import os
import time
from urllib.parse import urljoin

import requests
from beartype import beartype

from .config import get_value as get_config
from .exceptions import DeepOriginException

__all__ = [
    "get_tokens",
    "cache_tokens",
    "read_cached_tokens",
    "remove_cached_tokens",
    "authenticate",
    "refresh_tokens",
]


@beartype
def tokens_exist() -> bool:
    """Check if the cached API tokens exist"""
    config = get_config()
    return os.path.isfile(config.api_tokens_filename)


@functools.cache
def get_tokens(
    *,
    refresh: bool = True,
) -> dict:
    """Get access token for accessing the Deep Origin API

    Gets tokens to access Deep Origin API. If the tokens exist
    on disk, those tokens are used first. On first use (within
    a python session), tokens are refreshed if refresh_tokens
    is `True`. On subsequent uses, tokens are not refreshed
    even if refresh_tokens is `True`.

    Args:
        refresh: whether to refresh the token (default: `True`)

    Returns:
        :obj:`tuple`: API access and refresh tokens
    """

    if tokens_exist():
        # tokens exist on disk
        tokens = read_cached_tokens()

        if refresh:
            tokens["access"] = refresh_tokens(tokens["refresh"])

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
    config = get_config()

    os.makedirs(os.path.dirname(config.api_tokens_filename), exist_ok=True)
    with open(config.api_tokens_filename, "w") as file:
        json.dump(tokens, file)


def read_cached_tokens() -> dict[str, str]:
    """Read cached API tokens"""
    config = get_config()
    with open(config.api_tokens_filename, "r") as file:
        tokens = json.load(file)
    return tokens


def remove_cached_tokens():
    """Remove cached API tokens"""
    config = get_config()
    if os.path.isfile(config.api_tokens_filename):
        os.remove(config.api_tokens_filename)


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
            "To connect to the Deep Origin platform, "
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
                "Sign in to the Deep Origin platform failed. Please try again."
            )
        time.sleep(sign_in_poll_interval_sec)

    response_json = response.json()
    api_access_token = response_json["access_token"]
    api_refresh_token = response_json["refresh_token"]

    tokens = dict(access=api_access_token, refresh=api_refresh_token)

    cache_tokens(tokens)

    return tokens


def refresh_tokens(api_refresh_token: str) -> str:
    """Refresh the access token for the DO platform

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
