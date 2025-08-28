"""this module handles authentication actions and interactions
with tokens"""

import json
import os
import time
from typing import Optional
from urllib.parse import urljoin

from beartype import beartype
import jwt
from jwt.algorithms import RSAAlgorithm
import requests

from deeporigin.config import get_value as get_config
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils.constants import ENVS
from deeporigin.utils.core import _get_api_tokens_filepath, read_cached_tokens

__all__ = [
    "get_tokens",
    "cache_tokens",
    "remove_cached_tokens",
    "authenticate",
    "refresh_tokens",
]

AUTH_DOMAIN = {
    "prod": "https://formicbio.us.auth0.com",
    "staging": "https://formicbio.us.auth0.com",
    "edge": "https://edge.login.deeporigin.io",
}

AUTH_DEVICE_CODE_ENDPOINT = "/oauth/device/code"
AUTH_TOKEN_ENDPOINT = "/oauth/token"
AUTH_AUDIENCE = "https://os.deeporigin.io/api"


AUTH_GRANT_TYPE = "urn:ietf:params:oauth:grant-type:device_code"

AUTH_CLIENT_ID = {
    "prod": "m3iyUcrANcIap2ogzWKpnYxCNujOrW3s",
    "staging": "2AMGd2bJnKjMtd7QBvJYlGPqb9vtntsY",
    "edge": "jbYq4bkeX2UJb1bePH1ci172KbUMVyak",
}


@beartype
def tokens_exist() -> bool:
    """Check if the cached API tokens exist"""

    return os.path.isfile(_get_api_tokens_filepath())


def get_tokens(
    *,
    refresh: bool = False,
) -> dict:
    """Get access token for accessing the Deep Origin API

    Gets tokens to access Deep Origin API.


    If an access token exists in the ENV, then it is used before
    anything else. If not, then tokens file is
    checked for access tokens, and used if they exist.

    Args:
        refresh: whether to refresh the token (default: `False`)

    Returns:
        API access and refresh tokens
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

    # check if the access token is expired
    if is_token_expired(decode_access_token(tokens["access"])):
        tokens["access"] = refresh_tokens(tokens["refresh"])
        cache_tokens(tokens)

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


def remove_cached_tokens() -> None:
    """Remove cached API tokens"""

    if os.path.isfile(_get_api_tokens_filepath()):
        os.remove(_get_api_tokens_filepath())


@beartype
def authenticate(
    *,
    env: Optional[ENVS] = None,
    save_tokens: bool = True,
) -> dict:
    """Get an access token for use with the Deep Origin API.
    The tokens are also cached to file

    Returns:
        :obj:`tuple`: API access token, API refresh token
    """

    if env is None:
        env = get_config()["env"]

    # Get a link for the user to sign into the Deep Origin platform

    body = {
        "client_id": AUTH_CLIENT_ID[env],
        "scope": "offline_access",
        "audience": AUTH_AUDIENCE,
    }

    response = requests.post(AUTH_DOMAIN[env] + AUTH_DEVICE_CODE_ENDPOINT, json=body)
    response.raise_for_status()
    response_json = response.json()
    device_code = response_json["device_code"]
    user_code = response_json["user_code"]
    verification_url = response_json["verification_uri_complete"]
    sign_in_poll_interval_sec = response_json["interval"]

    # Prompt the user to sign into the Deep Origin platform
    print(
        (
            "To connect to the Deep Origin Platform API, "
            f"navigate your browser to \n\n{verification_url}\n\n"
            f'and verify the confirmation code is "{user_code}", '
            'and click the "Confirm" button.'
        )
    )

    body = {
        "grant_type": AUTH_GRANT_TYPE,
        "device_code": device_code,
        "client_id": AUTH_CLIENT_ID[env],
    }
    # Wait for the user to sign into the Deep Origin platform
    while True:
        response = requests.post(AUTH_DOMAIN[env] + AUTH_TOKEN_ENDPOINT, json=body)
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

    tokens = {
        "access": api_access_token,
        "refresh": api_refresh_token,
    }

    if save_tokens:
        cache_tokens(tokens)

    return tokens


@beartype
def refresh_tokens(api_refresh_token: str, *, env: Optional[ENVS] = None) -> str:
    """Refresh the access token for the Deep Origin OS

    Args:
        api_refresh_token: API refresh token

    Returns:
        new API access token
    """
    config = get_config()

    if env is None:
        env = get_config()["env"]

    body = {
        "grant_type": "refresh_token",
        "client_id": config.auth_client_id,
        "client_secret": config.auth_client_secret,
        "refresh_token": api_refresh_token,
    }
    response = requests.post(AUTH_DOMAIN[env] + AUTH_TOKEN_ENDPOINT, json=body)
    response.raise_for_status()
    response_json = response.json()
    api_access_token = response_json["access_token"]

    return api_access_token


@beartype
def is_token_expired(token: dict) -> bool:
    """
    Check if the JWT token is expired. The token is expected to have an 'exp' field as a Unix timestamp. This dict can be obtained from the `decode_access_token` function.

    Args:
        token (dict): The JWT token with an 'exp' field as a Unix timestamp.

    Returns:
        bool: True if the token is expired, False otherwise.
    """
    # Get the expiration time from the token, defaulting to 0 if not found.
    exp_time = token.get("exp", 0)
    current_time = time.time()  # Get current time in seconds since the epoch.

    # If current time is greater than the expiration time, it's expired.
    return current_time > exp_time


@beartype
def get_public_keys(env: Optional[ENVS] = None) -> list[dict]:
    """get public keys from public endpoint"""

    if env is None:
        env = get_config()["env"]

    jwks_url = urljoin(AUTH_DOMAIN[env], ".well-known/jwks.json")
    data = requests.get(jwks_url).json()
    return data["keys"]


@beartype
def decode_access_token(token: Optional[str] = None) -> dict:
    """decode access token into human readable data"""

    if token is None:
        tokens = get_tokens()
        token = tokens["access"]

    # Get the JWT header to extract the Key ID
    unverified_header = jwt.get_unverified_header(token)
    kid = unverified_header["kid"]

    # Get the public key using the Key ID
    public_keys = get_public_keys()
    for key in public_keys:
        if key["kid"] == kid:
            public_key = RSAAlgorithm.from_jwk(key)
            break
        raise DeepOriginException(f"Key ID {kid} not found in JWKS.")

    # Decode the JWT using the public key
    return jwt.decode(
        token,
        public_key,
        algorithms=["RS256"],
        options={
            "verify_aud": False,  # matches what platform does
            "verify_exp": False,  # we want to decode this no matter what, because we'll check the expiration in the caller
        },
    )
