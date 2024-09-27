import pytest
from deeporigin.config import get_value
from deeporigin.data_hub import api
from mock_client import MockClient

from tests.utils import _run_cli_command


@pytest.fixture(scope="session", autouse=True)
def config(pytestconfig):
    """this fixture performs some setup tasks
    before all tests are run, and runs only once"""

    data = dict()

    # set up client
    if pytestconfig.getoption("client") == "mock":
        data["client"] = MockClient()

        # unpack mock data from client
        data["folders"] = data["client"].folders
        data["databases"] = data["client"].databases
        data["rows"] = data["client"].rows
        data["file"] = data["client"].file
    else:
        client = api._get_default_client()
        data["client"] = client

        # if we're going to be making requests to a live
        # instance, we need to make sensible requests
        databases = api.list_rows(row_type="database")
        data["databases"] = [db.hid for db in databases]
        rows = api.list_rows(row_type="row")
        data["rows"] = [row.id for row in rows]

        # get a list of all files
        files = api.list_files()
        if len(files) > 0:
            data["file"] = files[0].file

    # tests run on yield
    yield data

    # teardown tasks, if any


def test_set_config(config):
    org_id = "foo_1234"
    try:
        org_id = get_value()["organization_id"]
    except Exception:
        # this may not work on github actions
        pass

    stdout = _run_cli_command(
        ["config", "set", "organization_id", "foo_3423"],
        config["client"],
    )

    assert "foo_3423" in stdout, "Failed to set organization_id"

    # reset to what it was before
    stdout = _run_cli_command(
        ["config", "set", "organization_id", org_id],
        config["client"],
    )

    assert org_id in stdout, "Failed to set organization_id"
