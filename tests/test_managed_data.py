"""this tests low level functions in the data API"""

from dataclasses import asdict

import pandas as pd
import pytest
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data import _api, api
from deeporigin.managed_data.client import (
    DeepOriginClient,
    MockClient,
    file_description,
)
from deeporigin.managed_data.schema import (
    DescribeFileResponse,
    ListRowsResponse,
    RowDescription,
)

# constants
dataframe_attr_keys = {
    "file_ids",
    "id",
    "primary_key",
    "reference_ids",
}


describe_file_keys = {
    "id",
    "uri",
    "name",
    "status",
    "contentLength",
    "contentType",
}

list_row_keys = {"id", "hid", "type", "name", "parentId"}

list_database_rows_keys = {
    "createdByUserDrn",
    "dateCreated",
    "dateUpdated",
    "editedByUserDrn",
    "fields",
    "hid",
    "hidNum",
    "id",
    "parentId",
    "submissionStatus",
    "type",
    "validationStatus",
}


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


def test_list_rows(config):
    rows = _api.list_rows(
        parent_id=config["databases"][0],
        client=config["client"],
    )

    for row in rows:
        assert set(row.keys()) == list_row_keys
        assert row["type"] == "row"


def test_list_rows_root_parent(config):
    root = _api.list_rows(
        parent_is_root=True,
        client=config["client"],
    )

    assert len(root) == 1, "Expected there to be exactly one root"

    assert root[0]["parentId"] is None, "Expected root to have no parent"


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

    keys_with_fields = set(asdict(RowDescription()).keys())

    assert set(row.keys()) == keys_with_fields

    row = _api.describe_row(
        config["rows"][0],
        client=config["client"],
        fields=False,
    )

    assert set(row.keys()) == keys_with_fields.difference({"fields"})


def test_convert_id_format(config):
    conversions = _api.convert_id_format(
        hids=[config["rows"][0]],
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

    for row in rows:
        assert set(row.keys()).difference({"name"}) == list_database_rows_keys


def test_get_dataframe(config):
    df = api.get_dataframe(
        config["databases"][0],
        client=config["client"],
    )

    assert isinstance(df, pd.DataFrame), "Expected return type to be a pandas Dataframe"

    assert (
        set(df.attrs.keys()) == dataframe_attr_keys
    ), f"Expected to find a dictionary in `df.attrs` with these keys: {dataframe_attr_keys}, instead found a dictionary with these keys: {df.attrs.keys()}"

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
    file_id = config["file"]["id"]

    if config["file"]["name"] != "placeholder":
        _api.download_file(file_id)


def test_describe_file(config):
    file_id = config["file"]["id"]

    data = _api.describe_file(file_id, client=config["client"])

    assert isinstance(data, dict), "Expected response to be a dict"

    DescribeFileResponse(**data)
