"""Tests for environment variable fallbacks in platform Client."""

import os
from typing import Generator

import pytest

from deeporigin.platform import Client


@pytest.fixture(autouse=True)
def clear_env() -> Generator[None, None, None]:
    """Clear relevant env vars for each test to avoid cross-test contamination."""
    keys = ["DEEPORIGIN_TOKEN", "DEEPORIGIN_ENV", "DEEPORIGIN_ORG_KEY"]
    old = {k: os.environ.get(k) for k in keys}
    for k in keys:
        os.environ.pop(k, None)
    try:
        yield
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def test_env_fallbacks_with_defaults() -> None:
    os.environ["DEEPORIGIN_TOKEN"] = "tok_abc"
    os.environ["DEEPORIGIN_ORG_KEY"] = "org_123"
    # DEEPORIGIN_ENV intentionally not set -> defaults to prod

    client = Client()

    assert client.token == "tok_abc"
    assert client.org_key == "org_123"
    assert client.env == "prod"
    assert client.api_endpoint.endswith("deeporigin.io")


def test_kwarg_overrides_env() -> None:
    os.environ["DEEPORIGIN_TOKEN"] = "tok_env"
    os.environ["DEEPORIGIN_ORG_KEY"] = "org_env"
    os.environ["DEEPORIGIN_ENV"] = "staging"

    client = Client(token="tok_kw", org_key="org_kw", env="edge")

    assert client.token == "tok_kw"
    assert client.org_key == "org_kw"
    assert client.env == "edge"


def test_missing_token_raises() -> None:
    os.environ["DEEPORIGIN_ORG_KEY"] = "org_123"
    with pytest.raises(ValueError):
        Client()


def test_missing_org_key_raises() -> None:
    os.environ["DEEPORIGIN_TOKEN"] = "tok_123"
    with pytest.raises(ValueError):
        Client()
