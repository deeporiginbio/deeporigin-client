import io
import json
from contextlib import redirect_stderr, redirect_stdout

import pytest
from beartype import beartype
from deeporigin import cli
from deeporigin.managed_data import _api
from deeporigin.managed_data.client import (
    Client,
    DeepOriginClient,
    MockClient,
    file_description,
)

describe_file_inputs = [
    dict(
        file_id="_file:placeholder",
        error="The requested resource was not found",
    ),
]


@pytest.fixture(scope="session", autouse=True)
def config(pytestconfig):
    """this fixture performs some setup tasks
    before all tests are run, and runs only once"""

    data = dict()

    # set up client
    if pytestconfig.getoption("client") == "mock":
        data["client"] = MockClient()

        data["databases"] = ["db-sample"]
        data["rows"] = ["row-1"]
        data["file"] = file_description()
    else:
        client = DeepOriginClient()
        client.authenticate()
        data["client"] = client

        # if we're going to be making requests to a live
        # instance, we need to make sensible requests
        databases = _api.list_rows(row_type="database")
        data["databases"] = [db["hid"] for db in databases]
        rows = _api.list_rows(row_type="row")
        data["rows"] = [row["id"] for row in rows]

        # get a list of all files
        files = _api.list_files()
        if len(files) > 0:
            data["file"] = files[0]["file"]

    # tests run on yield
    yield data

    # teardown tasks, if any


def test_data(config):
    stdout = _run_cli_command(
        ["data"],
        config["client"],
    )

    assert "List data in managed data on Deep Origin" in stdout, "Unexpected output"


def test_describe_file(config):
    file_id = config["file"]["id"]

    stdout = _run_cli_command(
        ["data", "describe-file", file_id],
        config["client"],
    )
    assert file_id in stdout, "Expected to see file_id in output"


def test_describe_row(config):
    row_id = config["rows"][0]

    stdout = _run_cli_command(
        ["data", "describe-row", row_id],
        config["client"],
    )
    assert row_id in stdout, "Expected to see row_id in output"

    # try JSON output
    stdout = _run_cli_command(
        ["data", "describe-row", row_id, "--json"],
        config["client"],
    )

    # check that we can parse into JSON
    data = json.loads(stdout)

    assert isinstance(data, dict), "Expected data to be a dictionary"

    assert data["id"] == row_id, "Expected row ID to match"


def test_list_rows(config):
    db_id = config["databases"][0]

    stdout = _run_cli_command(
        ["data", "list-rows", db_id],
        config["client"],
    )
    assert "Parent ID" in stdout, "Unexpected output"

    # JSON

    stdout = _run_cli_command(
        ["data", "list-rows", db_id, "--json"],
        config["client"],
    )
    # check that we can parse into JSON
    data = json.loads(stdout)
    assert isinstance(data, list), "Expected data to be a list"


def test_ls(config):
    _run_cli_command(
        ["data", "ls"],
        config["client"],
    )


def test_show_db(config):
    _run_cli_command(
        ["data", "show-db", config["databases"][0]],
        config["client"],
    )


def test_row(config):
    _run_cli_command(
        ["data", "row", config["rows"][0]],
        config["client"],
    )


@beartype
def _run_cli_command(argv: list[str], client: Client) -> str:
    """helper function to run a CLI command, parse output and return"""
    stdout = io.StringIO()
    stderr = io.StringIO()

    with redirect_stdout(stdout), redirect_stderr(stderr):
        with cli.App(argv=argv) as app:
            app.client = client
            app.run()

    stdout = stdout.getvalue().strip()
    return stdout
