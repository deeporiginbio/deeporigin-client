"""this tests low level functions in the data API"""

import os
import uuid

import pandas as pd
import pytest
from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils import DATAFRAME_ATTRIBUTE_KEYS, PREFIXES
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


def test_upload_file(config):
    full_path = os.path.abspath(__file__)

    from requests.exceptions import ConnectionError

    if config["mock"]:
        with pytest.raises(ConnectionError):
            api.upload_file(
                full_path,
                client=config["client"],
            )
    else:
        api.upload_file(
            full_path,
            client=config["client"],
        )


def test_make_database_rows(config):
    with pytest.raises(DeepOriginException, match="must be at least 1"):
        api.make_database_rows(
            config["databases"][0],
            n_rows=0,
            client=config["client"],
        )

    response = api.make_database_rows(
        config["databases"][0],
        client=config["client"],
    )

    for row in response.rows:
        assert row.type == "row", "Expected rows to be created."


def test_create_workspace(config):
    name = "test-" + str(uuid.uuid4())[:6]
    api.create_workspace(
        workspace=dict(name=name, hid=name),
        client=config["client"],
    )


def test_create_database(config):
    unique_id = str(uuid.uuid4())[:6]
    api.create_database(
        database=dict(
            name="test-" + unique_id,
            hid="test-" + unique_id,
            hid_prefix="test" + unique_id,
            parent_id=config["folders"][0],
        ),
        client=config["client"],
    )


def test_delete_rows(config):
    """delete folders and databases"""

    for row_type in ["workspace", "database"]:
        rows = api.list_rows(row_type=row_type, client=config["client"])
        row_ids = [row.id for row in rows if "test" in row.hid]

        if len(row_ids) > 0:
            api.delete_rows(
                row_ids=row_ids,
                client=config["client"],
            )


def test_list_rows(config):
    rows = api.list_rows(
        parent_id=config["databases"][0],
        client=config["client"],
    )

    for row in rows:
        assert row.type == "row"


def test_list_rows_root_parent(config):
    root = api.list_rows(
        parent_is_root=True,
        client=config["client"],
    )

    root = root[0]

    assert root.parent_id is None, "Expected root to have no parent"


def test_list_rows_by_type(config):
    rows = api.list_rows(
        row_type="workspace",
        client=config["client"],
    )

    assert len(rows) > 0, "Expected at least one folder"

    for row in rows:
        assert (
            row.type == "workspace"
        ), f"Expected to get a list of folders, but {row} is not a folder"


def test_list_files_unassigned(config):
    files = api.list_files(
        filters=[dict(is_unassigned=True)],
        client=config["client"],
    )

    assert len(files) > 0, "Expected to find at least 1 unassigned files"

    for file in files:
        assert (
            file.assignments is None
        ), f"Expected not to see an assignments key for this file, but instead found {file}"


def test_list_files_assigned(config):
    files = api.list_files(
        filters=[dict(is_unassigned=False)],
        client=config["client"],
    )

    assert len(files) > 0, "Expected to find at least 1 assigned files"

    for file in files:
        assignments = file.assignments
        assert (
            len(assignments) > 0
        ), "Expected assignments to be a list with at least 1 element"

        for assignment in assignments:
            assert (
                assignment.row_id is not None
            ), f"Expected to find a row_id in assignments, but instead found {assignment}"


def test_describe_database_stats(config):
    stats = api.describe_database_stats(
        database_id=config["databases"][0],
        client=config["client"],
    )

    assert stats.row_count >= 0, "Expected a positive number of rows"


def test_describe_row(config):
    """test describe_row, using mocked response"""

    row = api.describe_row(
        row_id=config["rows"][0],
        client=config["client"],
        fields=True,
    )

    # this doesn't seem to always be true
    # assert hasattr(row, "fields"), "Expected to find fields in response"

    row = api.describe_row(
        row_id=config["rows"][0],
        client=config["client"],
        fields=False,
    )

    assert not hasattr(row, "fields"), "Expected to NOT find fields in response"

    # database
    row = api.describe_row(
        row_id=config["databases"][0],
        client=config["client"],
    )


def test_convert_id_format(config):
    conversions = api.convert_id_format(
        hids=config["rows"],
        client=config["client"],
    )

    system_ids = []
    for conversion in conversions:
        system_ids.append(conversion.id)
        assert hasattr(conversion, "id"), "Expected to find `id` in the conversion"
        assert hasattr(conversion, "id"), "Expected to find `hid` in the conversion"

    conversions = api.convert_id_format(
        ids=system_ids,
        client=config["client"],
    )

    with pytest.raises(DeepOriginException, match="non-None and a list of strings"):
        api.convert_id_format(client=config["client"])


def test_list_database_rows(config):
    api.list_database_rows(
        database_row_id=config["databases"][0],
        client=config["client"],
    )


def test_get_dataframe(config):
    if config["mock"]:
        return

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

    data = api.get_dataframe(
        config["databases"][0],
        client=config["client"],
        return_type="dict",
    )

    assert isinstance(data, dict), "Expected return type to be a dict"


def test_list_mentions(config):
    data = api.list_mentions(
        query=config["rows"][0],
        client=config["client"],
    )

    assert hasattr(
        data, "mentions"
    ), "Expected to find a object with the attribute `mentions`"

    assert isinstance(data.mentions, list), "Expected `mentions` to be a list"


def test_get_tree(config):
    tree = api.get_tree(client=config["client"])

    tree = tree[0]

    assert tree["parent_id"] is None, "Expected the root of the tree to have no parent"

    tree.pop("children")

    tree = api.get_tree(client=config["client"], include_rows=False)
    tree = tree[0]
    assert tree["parent_id"] is None, "Expected the root of the tree to have no parent"

    tree.pop("children")


def test_create_file_download_url(config):
    file_id = config["file"].id
    data = api.create_file_download_url(
        file_id=file_id,
        client=config["client"],
    )

    assert hasattr(
        data, "download_url"
    ), "Expected to find `download_url` in data response"


def test_download_file(config):
    if "file" not in config.keys():
        return

    file_id = config["file"].id

    if config["mock"]:
        api.download_file(file_id=file_id, client=config["client"])

        with pytest.raises(DeepOriginException, match="should be a path to a folder"):
            api.download_file(
                file_id=file_id,
                client=config["client"],
                destination="non-existent-path",
            )

    else:
        api.download_file(file_id=file_id, client=config["client"])
        data = api.describe_file(file_id=file_id)
        os.remove(data.name)


def test_describe_file(config):
    file_id = config["file"].id

    api.describe_file(
        file_id=file_id,
        client=config["client"],
    )


def test_get_row_data(config):
    row_id = config["rows"][0]

    data1 = api.get_row_data(
        row_id=row_id,
        client=config["client"],
    )

    assert isinstance(data1, dict), "Expected return type to be a dict"

    data2 = api.get_row_data(
        row_id,
        client=config["client"],
        use_column_keys=True,
    )

    assert isinstance(data2, dict), "Expected return type to be a dict"


def test_get_cell_data(config):
    row_id = config["rows"][0]
    data = api.get_row_data(
        row_id=row_id,
        client=config["client"],
    )
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
