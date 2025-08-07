"""helper module to set up tests"""

from contextlib import redirect_stderr, redirect_stdout
import io

from beartype import beartype
from mock_client import MockClient
import pytest

from deeporigin import cli


@beartype
def _run_cli_command(
    argv: list[str],
    *,
    client,
) -> str:
    """helper function to run a CLI command, parse output and return"""
    stdout = io.StringIO()
    stderr = io.StringIO()

    with redirect_stdout(stdout), redirect_stderr(stderr):
        with cli.App(argv=argv) as app:
            app.client = client
            app.run()

    return stdout.getvalue().strip()


@pytest.fixture(scope="session", autouse=True)
def config(pytestconfig):
    """this fixture performs some setup tasks
    before all tests are run, and runs only once"""

    data = dict(
        databases=["kitchen-sink"],
        stash=pytestconfig.getoption("responses") == "stash",
    )

    # set up client
    if pytestconfig.getoption("client") == "mock":
        client = MockClient()
        data["mock"] = True

    else:
        data["mock"] = False
        client = None

    data["client"] = client

    # tests run on yield
    yield data


@pytest.fixture(scope="session", autouse=True)
def minimal_config(pytestconfig):
    """this fixture is a minimal config for working
    with DBs on a live instance where the config happens
    in the test"""

    data = dict(
        stash=pytestconfig.getoption("responses") == "stash",
    )

    # set up client
    if pytestconfig.getoption("client") == "mock":
        data["client"] = MockClient()
        data["mock"] = True

    else:
        data["mock"] = False
        data["client"] = None

    # tests run on yield
    yield data
