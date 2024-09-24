import io
from contextlib import redirect_stderr, redirect_stdout

import pytest
from beartype import beartype
from deeporigin import cli
from deeporigin.data_hub import api
from mock_client import MockClient


@beartype
def _run_cli_command(argv: list[str], client) -> str:
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

    data = dict()

    # set up client
    if pytestconfig.getoption("client") == "mock":
        data["client"] = MockClient()
        data["mock"] = True

        # unpack mock data from client
        data["folders"] = data["client"].folders
        data["databases"] = data["client"].databases
        data["rows"] = data["client"].rows
        data["file"] = data["client"].file
    else:
        data["mock"] = False
        client = api._get_default_client()

        data["client"] = client

        # if we're going to be making requests to a live
        # instance, we need to make sensible requests
        rows = api.list_rows()

        data["databases"] = [row.hid for row in rows if row.type == "database"]
        data["folders"] = [row.hid for row in rows if row.type == "workspace"]
        data["rows"] = [row.hid for row in rows if row.type == "row"]

        # get a list of all files
        files = api.list_files()
        if len(files) > 0:
            data["file"] = files[0].file

    # tests run on yield
    yield data

    # teardown tasks, if any
