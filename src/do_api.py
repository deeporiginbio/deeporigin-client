import functools
import json
import os
import time

import requests

from .config import get_value as get_config
from .exceptions import DeepOriginException

__all__ = [
    "get_do_api_tokens",
    "cache_do_api_tokens",
    "read_cached_do_api_tokens",
    "remove_cached_do_api_tokens",
    "sign_into_do_platform",
    "refresh_access_to_do_platform",
]


@functools.cache
def get_do_api_tokens() -> tuple[str, str]:
    """Get a token for accessing the Deep Origin API

    If the user already has a token, refresh it. If not, sign into the Deep Origin platform.
    Then cache the new or refreshed token.

    Returns:
        :obj:`tuple`: API access and refresh tokens
    """
    config = get_config()

    if os.path.isfile(config.api_tokens_filename):
        tokens = read_cached_do_api_tokens()
        refresh_token = tokens["refresh"]

        access_token = refresh_access_to_do_platform(refresh_token)

    else:
        print("No cached file, signing in...")
        access_token, refresh_token = sign_into_do_platform()

    return access_token, refresh_token


def cache_do_api_tokens(access_token: str, refresh_token: str) -> None:
    """Save access and refresh tokens to a local file, for example, to
    enable variables/secrets to be regularly pulled without the user
    needing to regularly re-login.

    Args:
        access_token (:obj:`str`): access token
        refresh_token (:obj:`str`): refresh token
    """
    config = get_config()

    os.makedirs(os.path.dirname(config.api_tokens_filename), exist_ok=True)
    with open(config.api_tokens_filename, "w") as file:
        tokens = {
            "access": access_token,
            "refresh": refresh_token,
        }
        json.dump(tokens, file)


def read_cached_do_api_tokens() -> dict[str, str]:
    """Read cached API tokens"""
    config = get_config()
    with open(config.api_tokens_filename, "r") as file:
        tokens = json.load(file)
    return tokens


def remove_cached_do_api_tokens():
    """Remove cached API tokens"""
    config = get_config()
    if os.path.isfile(config.api_tokens_filename):
        os.remove(config.api_tokens_filename)


def sign_into_do_platform() -> tuple[str, str]:
    """Get an access token for use with the Deep Origin API.
    The tokens are also cached to file

    Returns:
        :obj:`tuple`: API access token, API refresh token
    """

    config = get_config()

    # Get a link for the user to sign into the Deep Origin platform
    endpoint = f"{config.auth_domain}{config.auth_device_code_endpoint}"
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
    endpoint = f"{config.auth_domain}{config.auth_token_endpoint}"
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

    return api_access_token, api_refresh_token


def refresh_access_to_do_platform(api_refresh_token: str) -> str:
    """Refresh the access token for the DO platform

    Args:
        api_refresh_token (:obj:`str`): API refresh token

    Returns:
        :obj:`str`: new API access token
    """
    config = get_config()

    endpoint = f"{config.auth_domain}{config.auth_token_endpoint}"
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
