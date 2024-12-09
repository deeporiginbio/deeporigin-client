"""experimental module to interact with the platform API"""

import concurrent.futures
import functools
from typing import Optional
from urllib.parse import urljoin

import diskcache as dc
import jwt
import requests
from beartype import beartype
from beartype.typing import Literal
from deeporigin import auth
from deeporigin.config import get_value
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils.core import _ensure_do_folder
from deeporigin.utils.network import _get_domain_name, check_for_updates
from jwt.algorithms import RSAAlgorithm
from tqdm import tqdm

check_for_updates()


HTTP_VERB = Literal["GET", "POST", "PUT", "PATCH", "DELETE"]


@beartype
def _make_request(
    endpoint: str,
    *,
    verb: HTTP_VERB = "GET",
) -> dict | list[dict]:
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

    response = requests.request(
        verb,
        url,
        headers=headers,
    )

    if response.status_code >= 400:
        raise DeepOriginException(
            title=f"Error: {response.status_code}",
            message=response.text,
            fix=url,
        )

    data = response.json()
    if "data" in data.keys():
        return data["data"]

    else:
        return data


@beartype
def delete_workstation(workstation: str):
    """delete a workstation"""

    return _make_request(
        f"computebenches/{_get_org_id()}/{workstation}:delete", verb="POST"
    )


@beartype
def delete_terminated_workstations():
    """delete all terminated workstation"""

    workstations = get_workstations()

    names = [
        ws["attributes"]["name"]
        for ws in workstations
        if ws["attributes"]["status"] == "TERMINATED"
    ]

    if len(names) == 0:
        print("No workstations to delete")
        return

    org_id = _get_org_id()

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for name in names:
            futures.append(
                executor.submit(
                    _make_request,
                    f"computebenches/{org_id}/{name}:delete",
                    "POST",
                )
            )

        for _ in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures),
            desc="Deleting terminated workstations...",
        ):
            pass


@beartype
def _get_org_id() -> str:
    value = get_value()
    return value["organization_id"]


@beartype
def resolve_user(user_id: str):
    """get details about a user in the platform"""

    endpoint = f"/organizations/{_get_org_id()}/users/{user_id}"

    return _make_request(endpoint)


@beartype
def whoami() -> dict:
    """get details about currently signed in user"""

    return _make_request("/users/me", verb="PATCH")


@beartype
def get_workstations() -> list[dict]:
    """get information about all workstations in the organization"""
    return _make_request(f"/computebenches/{_get_org_id()}")


@beartype
def get_user_name(user_id: str) -> str:
    """get user name from user ID"""

    CACHE_PATH = _ensure_do_folder() / "user_ids"

    cache = dc.Cache(CACHE_PATH)

    if cache.get(user_id) is not None:
        return cache.get(user_id)

    response = resolve_user(user_id)

    name = response["attributes"]["name"]
    cache.set(user_id, name)
    return name


def get_last_edited_user_name(row: dict):
    """given a row object, return a human readable user name
    of the user who last edited the row"""

    if getattr(row, "editedByUserDrn", None) is None:
        # fall back to created by user
        user_id = row.createdByUserDrn.split(":")[-1]
    else:
        user_id = row.editedByUserDrn.split(":")[-1]
    return get_user_name(user_id)


@beartype
def get_variables_and_secrets() -> list[dict]:
    """get variables and secrets for user and org"""

    # we don't need to provide a bench ID, so pass a placeholder
    response = _make_request(f"/computebenches/{_get_org_id()}/placeholder/secrets")

    return response["data"]


@functools.cache
@beartype
def get_public_keys() -> list[dict]:
    """get public keys from public endpoint"""

    jwks_url = urljoin(get_value()["auth_domain"], ".well-known/jwks.json")
    data = requests.get(jwks_url).json()
    return data["keys"]


@beartype
def decode_access_token(token: Optional[str] = None) -> dict:
    """decode access token into human readable data"""

    if token is None:
        tokens = auth.get_tokens()
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
        raise Exception(f"Key ID {kid} not found in JWKS.")

    # Decode the JWT using the public key
    return jwt.decode(
        token,
        public_key,
        algorithms=["RS256"],
        options={"verify_aud": False},  # matches what platform does
    )
