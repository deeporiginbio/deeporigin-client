"""Platform client module.

Provides the `Client` used to interact with the Deep Origin platform.

This module supports configuration via keyword arguments or the following
environment variables when keywords are omitted:

- `DEEPORIGIN_TOKEN`
- `DEEPORIGIN_ENV` (defaults to "prod" if not provided)
- `DEEPORIGIN_ORG_KEY`
"""

from beartype import beartype

from deeporigin.auth import get_tokens
from deeporigin.config import get_value
from deeporigin.utils.constants import API_ENDPOINT, ENVS


class Client:
    """Deep Origin client, to interact with the Deep Origin platform"""

    @beartype
    def __init__(
        self,
        *,
        token: str | None = None,
        env: ENVS | None = None,
        org_key: str | None = None,
    ):
        """Initialize a platform client.

        Args:
            token: API token. If not provided, falls back to `DEEPORIGIN_TOKEN`.
            env: Target environment (e.g., "prod", "staging"). If not provided,
                falls back to `DEEPORIGIN_ENV`, defaulting to "prod" when unset.
            org_key: Organization key. If not provided, falls back to
                `DEEPORIGIN_ORG_KEY`.

        Raises:
            ValueError: If a required value is missing after applying fallbacks.
            KeyError: If the resolved environment is not a known API endpoint.
        """

        tokens = get_tokens()
        resolved_token = token or tokens["access"]
        resolved_env = env or get_value()["env"]
        resolved_org_key = org_key or get_value()["org_key"]

        if not resolved_org_key:
            raise ValueError(
                "Missing org_key. Provide `org_key` or set DEEPORIGIN_ORG_KEY in the environment."
            )

        # Assign resolved values
        self.token = resolved_token
        self.org_key = resolved_org_key
        self.env = resolved_env
        self.api_endpoint = API_ENDPOINT[resolved_env]

    def __repr__(self):
        return f"DeepOrigin Platform Client(token={self.token[:5]}..., org_key={self.org_key}, env={self.env})"
