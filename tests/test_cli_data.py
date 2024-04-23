import io
from contextlib import redirect_stderr, redirect_stdout

import pytest
from deeporigin import cli
from deeporigin.managed_data import _api
from deeporigin.managed_data.client import (
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
        data["rows"] = ["sample-1"]
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
        data["rows"] = [row["hid"] for row in rows]

        # get a list of all files
        files = _api.list_files()
        if len(files) > 0:
            data["file"] = files[0]["file"]

    # tests run on yield
    yield data

    # teardown tasks, if any


def test_data():
    stdout = io.StringIO()
    stderr = io.StringIO()

    with redirect_stdout(stdout), redirect_stderr(stderr):
        with cli.App(argv=["data"]) as app:
            app.run()

    stdout = stdout.getvalue().strip()

    assert "List data in managed data on Deep Origin" in stdout, "Unexpected output"


def test_describe_file(config):
    stdout = io.StringIO()
    stderr = io.StringIO()

    file_id = config["file"]["id"]

    with redirect_stdout(stdout), redirect_stderr(stderr):
        with cli.App(argv=["data", "describe-file", file_id]) as app:
            app.client = config["client"]
            app.run()
