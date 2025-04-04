"""this tests functions used to interact with the data hub"""

import os
import shutil
import uuid

import pandas as pd
import pytest

from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from tests.utils import TEST_DB_NAME, TEST_PREFIX, TEST_WS_NAME, config  # noqa: F401

row_types = ["row", "workspace", "database"]


def test_upload_file(config):  # noqa: F811
    full_path = os.path.abspath(__file__)

    if config["mock"]:
        pass
    else:
        api.upload_file(
            full_path,
            client=config["client"],
            _stash=config["stash"],
        )


def test_make_database_rows(config):  # noqa: F811
    with pytest.raises(DeepOriginException, match="must be at least 1"):
        api.make_database_rows(
            config["databases"][0],
            n_rows=0,
            client=config["client"],
            _stash=config["stash"],
        )

    response = api.make_database_rows(
        config["databases"][0],
        client=config["client"],
        _stash=config["stash"],
    )

    for row in response.rows:
        assert row.type == "row", "Expected rows to be created."


def test_create_and_delete_database(config):  # noqa: F811
    api.create_database(
        name=TEST_DB_NAME,
        hid=TEST_DB_NAME,
        hid_prefix=TEST_DB_NAME,
        parent_id=config["folders"][0],
        client=config["client"],
        _stash=config["stash"],
    )

    api.delete_database(
        database_id=TEST_DB_NAME,
        _stash=config["stash"],
        client=config["client"],
    )


def test_create_and_delete_workspace(config):  # noqa: F811
    api.create_workspace(
        name=TEST_WS_NAME,
        client=config["client"],
        _stash=config["stash"],
    )

    api.delete_workspace(
        workspace_id=TEST_WS_NAME,
        client=config["client"],
        _stash=config["stash"],
    )


@pytest.mark.parametrize("row_type", row_types)
def test_list_rows(config, row_type):  # noqa: F811
    print(config["stash"])
    rows = api.list_rows(
        client=config["client"],
        _stash=config["stash"],
        row_type=row_type,
    )

    for row in rows:
        assert row.type == row_type


def test_list_rows_parent(config):  # noqa: F811
    rows = api.list_rows(
        parent_id=config["databases"][0],
        client=config["client"],
        _stash=config["stash"],
    )

    for row in rows:
        assert row.type == "row"


def test_list_rows_root_parent(config):  # noqa: F811
    root = api.list_rows(
        parent_is_root=True,
        client=config["client"],
        _stash=config["stash"],
    )

    root = root[0]

    assert "parentId" not in root.keys() or root.parentId is None, (
        "Expected root to have no parent"
    )


def test_list_rows_by_type(config):  # noqa: F811
    rows = api.list_rows(
        row_type="workspace",
        client=config["client"],
        _stash=config["stash"],
    )

    assert len(rows) > 0, "Expected at least one folder"

    for row in rows:
        assert row.type == "workspace", (
            f"Expected to get a list of folders, but {row} is not a folder"
        )


def test_list_files(config):  # noqa: F811
    """simple listing of all files with no args"""
    api.list_files(
        client=config["client"],
        _stash=config["stash"],
    )


def test_list_files_unassigned(config):  # noqa: F811
    files = api.list_files(
        is_unassigned=True,
        client=config["client"],
        _stash=config["stash"],
    )

    assert len(files) > 0, "Expected to find at least 1 unassigned files"

    for file in files:
        assert "assignments" not in file.keys() or file.assignments is None, (
            f"Expected not to see an assignments key for this file, but instead found {file}"
        )


def test_list_files_assigned(config):  # noqa: F811
    files = api.list_files(
        is_unassigned=False,
        client=config["client"],
        _stash=config["stash"],
    )

    assert len(files) > 0, "Expected to find at least 1 assigned files"

    for file in files:
        assignments = file.assignments
        assert len(assignments) > 0, (
            "Expected assignments to be a list with at least 1 element"
        )

        for assignment in assignments:
            assert assignment.rowId is not None, (
                f"Expected to find a rowId in assignments, but instead found {assignment}"
            )


def test_describe_database_stats(config):  # noqa: F811
    stats = api.describe_database_stats(
        database_id=config["databases"][0],
        client=config["client"],
        _stash=config["stash"],
    )

    assert stats.rowCount >= 0, "Expected a positive number of rows"


def test_describe_row(config):  # noqa: F811
    """test describe_row, using mocked response"""

    row = api.describe_row(
        row_id=config["rows"][0],
        client=config["client"],
        fields=True,
        _stash=config["stash"],
    )

    # this doesn't seem to always be true
    # assert hasattr(row, "fields"), "Expected to find fields in response"

    row = api.describe_row(
        row_id=config["rows"][0],
        client=config["client"],
        fields=False,
        _stash=config["stash"],
    )

    assert not hasattr(row, "fields"), "Expected to NOT find fields in response"

    # database
    row = api.describe_row(
        row_id=config["databases"][0],
        client=config["client"],
        _stash=config["stash"],
    )


def test_convert_id_format(config):  # noqa: F811
    conversions = api.convert_id_format(
        hids=config["rows"][:3],
        client=config["client"],
        _stash=config["stash"],
    )

    system_ids = []
    for conversion in conversions:
        system_ids.append(conversion.id)
        assert hasattr(conversion, "id"), "Expected to find `id` in the conversion"
        assert hasattr(conversion, "id"), "Expected to find `hid` in the conversion"

    conversions = api.convert_id_format(
        ids=system_ids,
        client=config["client"],
        _stash=config["stash"],
    )


def test_list_database_rows(config):  # noqa: F811
    api.list_database_rows(
        database_row_id=config["databases"][0],
        client=config["client"],
        _stash=config["stash"],
    )


def test_get_dataframe(config):  # noqa: F811
    if config["mock"]:
        return

    print(config["databases"][0])

    df = api.get_dataframe(
        config["databases"][0],
        reference_format="system-id",
        client=config["client"],
        _stash=config["stash"],
    )

    assert isinstance(df, pd.DataFrame), "Expected return type to be a pandas Dataframe"

    assert "Validation Status" in df.columns, (
        f"Expected to find a column called `Validation Status` in the dataframe. Instead, the columns in this dataframe are: {df.columns}"
    )

    data = api.get_dataframe(
        config["databases"][0],
        reference_format="system-id",
        client=config["client"],
        return_type="dict",
        _stash=config["stash"],
    )

    assert isinstance(data, dict), "Expected return type to be a dict"


def test_list_mentions(config):  # noqa: F811
    data = api.list_mentions(
        query=config["rows"][0],
        client=config["client"],
        _stash=config["stash"],
    )

    assert hasattr(data, "mentions"), (
        "Expected to find a object with the attribute `mentions`"
    )

    assert isinstance(data.mentions, list), "Expected `mentions` to be a list"


def test_get_tree(config):  # noqa: F811
    tree = api.get_tree(
        client=config["client"],
        _stash=config["stash"],
    )

    tree = tree[0]

    assert "parentId" not in tree.keys() or tree.parentId is None, (
        "Expected the root of the tree to have no parent"
    )

    tree.pop("children")

    tree = api.get_tree(client=config["client"], include_rows=False)
    tree = tree[0]
    assert "parentId" not in tree.keys() or tree.parentId is None, (
        "Expected the root of the tree to have no parent"
    )

    tree.pop("children")


def test_create_file_download_url(config):  # noqa: F811
    file_id = config["file"].id
    data = api.create_file_download_url(
        file_id=file_id,
        client=config["client"],
        _stash=config["stash"],
    )

    assert "downloadUrl" in data.keys(), (
        "Expected to find `downloadUrl` in data response"
    )


def test_download_files(config):  # noqa: F811
    if "file" not in config.keys():
        return

    file_id = config["file"].id

    if config["mock"]:
        with pytest.raises(DeepOriginException, match="should be a path for a folder"):
            api.download_files(
                file_ids=[file_id],
                client=config["client"],
                save_to_dir="non-existent-path",
            )

    else:
        destination = str(uuid.uuid4())[:5]
        destination = os.path.join(os.getcwd(), destination)
        os.makedirs(destination)

        api.download_files(
            file_ids=[file_id],
            client=config["client"],
            _stash=config["stash"],
            save_to_dir=destination,
        )

        # clean up
        shutil.rmtree(destination)


def test_describe_file(config):  # noqa: F811
    file_id = config["file"].id

    api.describe_file(
        file_id=file_id,
        client=config["client"],
        _stash=config["stash"],
    )


def test_get_row_data(config):  # noqa: F811
    rowId = config["rows"][0]

    data1 = api.get_row_data(
        row_id=rowId,
        client=config["client"],
        _stash=config["stash"],
    )

    assert isinstance(data1, dict), "Expected return type to be a dict"

    data2 = api.get_row_data(
        rowId,
        client=config["client"],
        use_column_keys=True,
        _stash=config["stash"],
    )

    assert isinstance(data2, dict), "Expected return type to be a dict"


def test_get_cell_data(config):  # noqa: F811
    row_id = config["rows"][0]
    data = api.get_row_data(
        row_id=row_id,
        client=config["client"],
        _stash=config["stash"],
    )
    column_name = list(data.keys())[0]
    data = api.get_cell_data(
        row_id=row_id,
        column_name=column_name,
        client=config["client"],
        _stash=config["stash"],
    )


def test_download_database(config):  # noqa: F811
    db_id = config["databases"][0]
    api.download_database(
        db_id,
        include_files=False,
        client=config["client"],
        _stash=config["stash"],
    )

    file_name = db_id + ".csv"
    assert os.path.exists(file_name), f"Expected to find a CSV file called {file_name}"
    os.remove(file_name)
