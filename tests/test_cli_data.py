import io
import json
from contextlib import redirect_stderr, redirect_stdout
from typing import Union

import pytest
from beartype import beartype
from deeporigin import cli
from deeporigin.managed_data import _api
from deeporigin.managed_data.client import (
    Client,
    DeepOriginClient,
)
from deeporigin.managed_data.schema import DescribeFileResponse, ListRowsResponse
from mock_client import MockClient

# this allows us to try every CLI command with both
# multiple options using pytest.mark.parametrize
JSON_OPTIONS = (
    ["--json"],
    [],
)
LIST_OPTIONS = (
    [],
    ["--rows"],
    ["--databases"],
    ["--files"],
    ["--workspaces"],
    ["--workspaces", "--databases"],
    ["--workspaces", "--rows"],
    ["--databases", "--rows"],
)


@pytest.fixture(scope="session", autouse=True)
def config(pytestconfig):
    """this fixture performs some setup tasks
    before all tests are run, and runs only once"""

    data = dict()

    # set up client
    if pytestconfig.getoption("client") == "mock":
        data["client"] = MockClient()

        # unpack mock data from client
        data["workspaces"] = data["client"].workspaces
        data["databases"] = data["client"].databases
        data["rows"] = data["client"].rows
        data["file"] = data["client"].file
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


@pytest.mark.parametrize("json_option", JSON_OPTIONS)
def test_describe_file(config, json_option):
    file_id = config["file"]["id"]

    stdout = _run_cli_command(
        ["data", "describe", file_id] + json_option,
        config["client"],
    )

    if json_option == ["--json"]:
        data = _check_json(stdout)
        assert data["id"] == file_id, "Expected file ID to match"


@pytest.mark.parametrize("option", JSON_OPTIONS)
def test_describe_row(config, option):
    row_id = config["rows"][0]

    stdout = _run_cli_command(
        ["data", "describe", row_id] + option,
        config["client"],
    )
    assert row_id in stdout, "Expected to see row_id in output"

    # check that we can parse into JSON
    if option == ["--json"]:
        data = _check_json(stdout)

        assert data["id"] == row_id, "Expected row ID to match"


@pytest.mark.parametrize("json_option", JSON_OPTIONS)
@pytest.mark.parametrize("list_option", LIST_OPTIONS)
def test_list(config, list_option, json_option):
    stdout = _run_cli_command(
        ["data", "list"] + list_option + json_option,
        config["client"],
    )

    # check that we can parse into JSON
    if json_option == ["--json"]:
        data = _check_json(stdout)

        for item in data:
            if list_option == ["--files"]:
                DescribeFileResponse(**(item["file"]))
            else:
                ListRowsResponse(**item)


@pytest.mark.parametrize("json_option", JSON_OPTIONS)
def test_show_db(config, json_option):
    stdout = _run_cli_command(
        ["data", "show", config["databases"][0]] + json_option,
        config["client"],
    )

    if json_option == ["--json"]:
        _check_json(stdout)


@beartype
def _run_cli_command(argv: list[str], client: Client) -> str:
    """helper function to run a CLI command, parse output and return"""
    stdout = io.StringIO()
    stderr = io.StringIO()

    with redirect_stdout(stdout), redirect_stderr(stderr):
        with cli.App(argv=argv) as app:
            app.client = client
            app.run()

    return stdout.getvalue().strip()


@beartype
def _check_json(txt: str) -> Union[dict, list[dict]]:
    """helper function to check that a string is valid JSON"""
    data = json.loads(txt)

    return data
