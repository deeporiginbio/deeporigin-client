"""helper module to set up tests"""

from contextlib import redirect_stderr, redirect_stdout
import io

from beartype import beartype
from mock_client import MockClient
import pytest

from deeporigin import cli
from deeporigin.data_hub import api

TEST_PREFIX = "tc-4Qzkrn57rM-"
TEST_DB_NAME = TEST_PREFIX + "db"
TEST_WS_NAME = TEST_PREFIX + "ws"


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


def clean_up_test_objects(test_prefix: str = None):
    """utility function to clean up objects that have been created by tests"""

    print("Cleaning up...")
    rows = api.list_rows()

    if test_prefix is None:
        test_prefix = TEST_PREFIX

    for row in rows:
        if test_prefix in row.hid:
            try:
                if row.type == "database":
                    api.delete_database(database_id=row.hid)
                elif row.type == "workspace":
                    api.delete_workspace(workspace_id=row.hid)
            except Exception:
                # it's possible it doesn't exist
                pass


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
        client = api._get_default_client()

    data["client"] = client

    # generate lists of objects used in tests
    rows = api.list_rows(
        row_type="workspace",
        client=client,
    )
    data["folders"] = [row.hid for row in rows]

    rows = api.list_rows(
        row_type="row",
        client=client,
    )
    data["rows"] = [row.hid for row in rows]

    # get a list of all files
    files = api.list_files(client=client)
    if len(files) > 0:
        data["file"] = files[0].file

    # tests run on yield
    yield data

    # clean up all the object we created
    if pytestconfig.getoption("client") != "mock":
        clean_up_test_objects()


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
        data["client"] = api._get_default_client()

    # tests run on yield
    yield data
