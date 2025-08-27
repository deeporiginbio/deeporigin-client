from beartype import beartype

from deeporigin.utils.constants import ENVS


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
        if env == "prod":
            self.api_endpoint = "https://api.deeporigin.io"
        elif env == "staging":
            self.api_endpoint = "https://api.staging.deeporigin.io"
        elif env == "edge":
            self.api_endpoint = "https://api.edge.deeporigin.io"
        else:
            raise ValueError(f"Invalid environment: {env}")

    def __repr__(self):
        return f"DeepOrigin Platform Client(token={self.token[:5]}..., org_key={self.org_key}, env={self.env})"
