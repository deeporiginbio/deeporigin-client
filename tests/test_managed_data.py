"""this tests low level functions in the data API"""

import os
import uuid

import pandas as pd
import pytest
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data import _api, api
from deeporigin.managed_data.client import (
    DeepOriginClient,
)
from deeporigin.managed_data.schema import (
    DATAFRAME_ATTRIBUTE_KEYS,
    CreateDatabaseResponse,
    CreateWorkspaceResponse,
    DescribeFileResponse,
    DescribeRowResponseDatabase,
    DescribeRowResponseRow,
    ListRowsResponse,
)
from deeporigin.utils import PREFIXES
from mock_client import MockClient


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
        data["workspaces"] = data["client"].workspaces
        data["databases"] = data["client"].databases
        data["rows"] = data["client"].rows
        data["file"] = data["client"].file
    else:
        data["mock"] = False
        client = DeepOriginClient()
        client.authenticate()
        data["client"] = client

        # if we're going to be making requests to a live
        # instance, we need to make sensible requests
        rows = _api.list_rows()

        data["databases"] = [row["hid"] for row in rows if row["type"] == "database"]
        data["workspaces"] = [row["hid"] for row in rows if row["type"] == "workspace"]
        data["rows"] = [row["hid"] for row in rows if row["type"] == "row"]

        # get a list of all files
        files = _api.list_files()
        if len(files) > 0:
            data["file"] = files[0]["file"]

    # tests run on yield
    yield data

    # teardown tasks, if any


def test_create_workspace(config):
    data = _api.create_workspace(
        name="test-" + str(uuid.uuid4())[:6],
        client=config["client"],
    )

    CreateWorkspaceResponse(**data)


def test_create_database(config):
    unique_id = str(uuid.uuid4())[:6]
    data = _api.create_database(
        name="test-" + unique_id,
        hid_prefix="test" + unique_id,
        client=config["client"],
        parent_id=config["workspaces"][0],
    )

    CreateDatabaseResponse(**data)


def test_delete_rows(config):
    """delete workspaces and databases"""

    for row_type in ["workspace", "database"]:
        rows = _api.list_rows(row_type=row_type, client=config["client"])
        row_ids = [row["id"] for row in rows if "test" in row["hid"]]

        if len(row_ids) > 0:
            _api.delete_rows(row_ids, client=config["client"])


def test_list_rows(config):
    rows = _api.list_rows(
        parent_id=config["databases"][0],
        client=config["client"],
    )

    for row in rows:
        assert row["type"] == "row"

        # check type
        ListRowsResponse(**row)


def test_list_rows_root_parent(config):
    root = _api.list_rows(
        parent_is_root=True,
        client=config["client"],
    )

    root = root[0]

    assert root["parentId"] is None, "Expected root to have no parent"

    # check type
    ListRowsResponse(**root)


def test_list_rows_by_type(config):
    rows = _api.list_rows(
        row_type="workspace",
        client=config["client"],
    )

    assert len(rows) > 0, "Expected at least one workspace"

    for row in rows:
        assert (
            row["type"] == "workspace"
        ), f"Expected to get a list of workspaces, but {row} is not a workspace"

        # check type
        ListRowsResponse(**row)


def test_list_files_unassigned(config):
    files = _api.list_files(
        is_unassigned=True,
        client=config["client"],
    )

    assert len(files) > 0, "Expected to find at least 1 unassigned files"

    for file in files:
        assert (
            "assignments" not in file.keys()
        ), f"Expected not to see an assignments key for this file, but instead found {file}"


def test_list_files_assigned(config):
    files = _api.list_files(
        is_unassigned=False,
        client=config["client"],
    )

    assert len(files) > 0, "Expected to find at least 1 assigned files"

    for file in files:
        assert (
            "assignments" in file.keys()
        ), f"Expected to see an assignments key for this file, but instead found {file}"

        assignments = file["assignments"]
        assert (
            len(assignments) > 0
        ), "Expected assignments to be a list with at least 1 element"

        for assignment in assignments:
            assert (
                "rowId" in assignment.keys()
            ), f"Expected to find a rowId in assignments, but instead found {assignment}"


def test_describe_database_stats(config):
    stats = _api.describe_database_stats(
        config["databases"][0],
        client=config["client"],
    )

    assert (
        "rowCount" in stats.keys()
    ), f"Expected stats to be a dictionary with a key called `rowCount`, instead got {stats}"

    assert stats["rowCount"] >= 0, "Expected a positive number of rows"


def test_describe_row(config):
    """test describe_row, using mocked response"""

    row = _api.describe_row(
        config["rows"][0],
        client=config["client"],
        fields=True,
    )

    assert "fields" in row.keys(), "Expected to find fields in response"

    # check if we can coerce into response type
    DescribeRowResponseRow(**row)

    row = _api.describe_row(
        config["rows"][0],
        client=config["client"],
        fields=False,
    )

    assert "fields" not in row.keys(), "Expected to NOT find fields in response"

    # check if we can coerce into response type
    DescribeRowResponseRow(**row)

    # database
    row = _api.describe_row(
        config["databases"][0],
        client=config["client"],
    )

    # check if we can coerce into response type
    DescribeRowResponseDatabase(**row)


def test_convert_id_format(config):
    conversions = _api.convert_id_format(
        hids=config["rows"],
        client=config["client"],
    )

    system_ids = []
    for conversion in conversions:
        system_ids.append(conversion["id"])
        assert {"id", "hid"} == set(conversion.keys())

    conversions = _api.convert_id_format(
        ids=system_ids,
        client=config["client"],
    )

    for conversion in conversions:
        assert {"id", "hid"} == set(conversion.keys())

    with pytest.raises(DeepOriginException, match="non-None and a list of strings"):
        _api.convert_id_format()


def test_list_database_rows(config):
    rows = _api.list_database_rows(
        config["databases"][0],
        client=config["client"],
    )

    # coerce into type
    for row in rows:
        DescribeRowResponseRow(**row)


def test_get_dataframe(config):
    df = api.get_dataframe(
        config["databases"][0],
        client=config["client"],
    )

    assert isinstance(df, pd.DataFrame), "Expected return type to be a pandas Dataframe"

    assert (
        set(df.attrs.keys()) == DATAFRAME_ATTRIBUTE_KEYS
    ), f"Expected to find a dictionary in `df.attrs` with these keys: {DATAFRAME_ATTRIBUTE_KEYS}, instead found a dictionary with these keys: {df.attrs.keys()}"

    assert (
        "Validation Status" in df.columns
    ), f"Expected to find a column called `Validation Status` in the dataframe. Instead, the columns in this dataframe are: {df.columns}"

    assert (
        df.attrs["primary_key"] in df.columns
    ), "Expected to find the primary key as a column"

    data = api.get_dataframe(
        config["databases"][0],
        client=config["client"],
        return_type="dict",
    )

    assert isinstance(data, dict), "Expected return type to be a dict"


def test_list_mentions(config):
    data = _api.list_mentions(
        config["rows"][0],
        client=config["client"],
    )

    assert (
        "mentions" in data.keys()
    ), "Expected to find a dictionary with the key `mentions`"

    assert isinstance(data["mentions"], list), "Expected `mentions` to be a list"


def test_get_tree(config):
    tree = api.get_tree(client=config["client"])

    tree = tree[0]

    assert tree["parentId"] is None, "Expected the root of the tree to have no parent"

    tree.pop("children")
    ListRowsResponse(**tree)

    tree = api.get_tree(client=config["client"], include_rows=False)
    tree = tree[0]
    assert tree["parentId"] is None, "Expected the root of the tree to have no parent"

    tree.pop("children")
    ListRowsResponse(**tree)


def test_create_file_download_url(config):
    file_id = config["file"]["id"]
    data = _api.create_file_download_url(
        file_id,
        client=config["client"],
    )

    assert (
        "downloadUrl" in data.keys()
    ), "Expected to find `downloadUrl` in data response"


def test_download_file(config):
    if "file" not in config.keys():
        return

    file_id = config["file"]["id"]

    if config["mock"]:
        _api.download_file(file_id, client=config["client"])

        with pytest.raises(DeepOriginException, match="should be a path to a folder"):
            _api.download_file(
                file_id,
                client=config["client"],
                destination="non-existent-path",
            )

    else:
        _api.download_file(file_id, client=config["client"])
        data = _api.describe_file(file_id)
        os.remove(data["name"])


def test_describe_file(config):
    file_id = config["file"]["id"]

    data = _api.describe_file(file_id, client=config["client"])

    assert isinstance(data, dict), "Expected response to be a dict"

    DescribeFileResponse(**data)


def test_get_row_data(config):
    row_id = config["rows"][0]

    data = api.get_row_data(row_id, client=config["client"])
    assert isinstance(data, dict), "Expected response to be a dict"

    data = api.get_row_data(
        row_id,
        client=config["client"],
        use_column_keys=True,
    )
    assert isinstance(data, dict), "Expected response to be a dict"


def test_get_cell_data(config):
    row_id = config["rows"][0]
    data = api.get_row_data(row_id, client=config["client"])
    column_name = list(data.keys())[0]
    data = api.get_cell_data(
        row_id=row_id,
        column_name=column_name,
        client=config["client"],
    )


def test_download_database(config):
    db_id = config["databases"][0]
    api.download_database(
        db_id,
        include_files=False,
        client=config["client"],
    )

    file_name = db_id + ".csv"
    assert os.path.exists(file_name), f"Expected to find a CSV file called {file_name}"
    os.remove(file_name)


def test_download(config):
    db_id = config["databases"][0]
    api.download(
        db_id,
        destination=os.getcwd(),
        client=config["client"],
    )

    file_name = db_id + ".csv"
    assert os.path.exists(file_name), f"Expected to find a CSV file called {file_name}"
    os.remove(file_name)

    api.download(
        f"{PREFIXES.DO}{db_id}",
        destination=os.getcwd(),
        client=config["client"],
    )

    file_name = db_id + ".csv"
    assert os.path.exists(file_name), f"Expected to find a CSV file called {file_name}"
    os.remove(file_name)
