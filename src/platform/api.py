"""module to interact with the platform API"""

import requests
from deeporigin import auth
from deeporigin.config import get_value


def resolve_user(user_id: str):
    """get details about a user in the platform"""

    tokens = auth.get_tokens(refresh=False)
    access_token = tokens["access"]

    value = get_value()
    org_id = value["organization_id"]

    headers = {
        "accept": "application/json",
        "authorization": f"Bearer {access_token}",
        "cache-control": "no-cache",
    }

    env = get_value()["env"]
    if env == "prod":
        url = f"https://os.deeporigin.io/api/organizations/{org_id}/users/{user_id}"
    else:
        f"https://{env}.deeporigin.io/api/organizations/{org_id}/users/{user_id}"

    response = requests.get(
        url,
        headers=headers,
    )

    return response.json()


def get_user_name(user_id: str):
    """get user name from user DRN"""

    response = resolve_user(user_id)

    return response["data"]["attributes"]["name"]


def get_last_edited_user_name(row):
    """given a row object, return a human readable user name
    of the user who last edited the row"""

    if row.edited_by_user_drn is None:
        # fall back to created by user
        user_id = row.created_by_user_drn.split(":")[-1]
    else:
        user_id = row.edited_by_user_drn.split(":")[-1]
    return get_user_name(user_id)
