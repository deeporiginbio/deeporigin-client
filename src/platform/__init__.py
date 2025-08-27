from beartype import beartype

from deeporigin.utils.constants import API_ENDPOINT, ENVS


class Client:
    """Deep Origin client, to interact with the Deep Origin platform"""

    @beartype
    def __init__(
        self,
        *,
        token: str,
        env: ENVS = "prod",
        org_key: str,
    ):
        self.token = token
        self.org_key = org_key
        self.env = env
        self.api_endpoint = API_ENDPOINT[env]

    def __repr__(self):
        return f"DeepOrigin Platform Client(token={self.token[:5]}..., org_key={self.org_key}, env={self.env})"
