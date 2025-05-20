"""this module contains tests for the DataFrame class"""

import uuid

from box import BoxList
import numpy as np
import pandas as pd
import pytest

from deeporigin.data_hub import api
from deeporigin.data_hub.dataframe import DataFrame
from deeporigin.exceptions import DeepOriginException
from tests.utils import clean_up_test_objects

TEST_PREFIX = "tc-"

NUMERIC_COLUMNS = ["float", "integer"]

SKIP_MSG = "Skipping this test because we can't test this well using a mocked client"


@pytest.fixture(scope="function", autouse=True)
def config(pytestconfig):
    """this config block sets up a database that we can operate on in the tests"""

    salt = str(uuid.uuid4())[:8]

    TEST_DB_NAME = TEST_PREFIX + salt + "-db"

    data = dict()
    if pytestconfig.getoption("client") == "mock":
        data["mock"] = True

        # Create the DataFrame
        data["df"] = DataFrame(
            pd.DataFrame(
                {
                    "integer": np.random.randint(1, 100, size=10),
                    "float": np.random.random(10),
                },
                index=[f"x-{i}" for i in range(1, 11)],
            )
        )

        data["df"].attrs["id"] = "placeholder"
        data["df"].attrs["metadata"] = {
            "id": "_database:jeWMK8pVA6XHD9Ht6tGrM",
            "type": "database",
            "hid": "xy",
            "name": "XY",
            "dateCreated": "2024-10-23 20:17:32.602123",
            "dateUpdated": "2024-10-23 20:17:32.602123",
            "createdByUserDrn": "drn:identity::user:google-apps|user@deeporigin.com",
            "parentId": "_workspace:gRrmQ9z14diV39ZwgvRx6",
            "hidPrefix": "xy",
            "cols": BoxList(
                [
                    {
                        "id": "_column:TbQYk64PIxSjk5QcaO95V",
                        "parentId": "_database:jeWMK8pVA6XHD9Ht6tGrM",
                        "name": "float",
                        "type": "float",
                        "cardinality": "one",
                        "configNumeric": {},
                    },
                    {
                        "id": "_column:Vesbs3ViQGJC47N8JgGRh",
                        "parentId": "_database:jeWMK8pVA6XHD9Ht6tGrM",
                        "name": "integer",
                        "type": "integer",
                        "cardinality": "one",
                        "configNumeric": {},
                    },
                ]
            ),
        }

    else:
        data["mock"] = False

        rows = api.list_rows(row_type="database")
        if TEST_DB_NAME in [row.hid for row in rows]:
            api.delete_database(database_id=TEST_DB_NAME)

        api.create_database(name=TEST_DB_NAME)

        # make columns for values and squared values
        for column in NUMERIC_COLUMNS:
            api.add_database_column(
                database_id=TEST_DB_NAME,
                name=column,
                type=column,
            )

            api.add_database_column(
                database_id=TEST_DB_NAME,
                name="sq_" + column,
                type=column,
            )

        # make rows. we'll only fill some rows
        response = api.make_database_rows(
            database_id=TEST_DB_NAME,
            n_rows=10,
        )
        row_ids = [row.hid for row in response.rows]

        api.set_data_in_cells(
            values=np.random.random(5),
            row_ids=row_ids[0:5],
            column_id="float",
            database_id=TEST_DB_NAME,
        )

        api.set_data_in_cells(
            values=np.random.randint(100, size=5),
            row_ids=row_ids[0:5],
            column_id="integer",
            database_id=TEST_DB_NAME,
        )

        data["db-name"] = TEST_DB_NAME

    # tests run on yield
    yield data

    # teardown code
    if pytestconfig.getoption("client") != "mock":
        clean_up_test_objects(TEST_DB_NAME)


@pytest.mark.parametrize("add_column", [True, False])
def test_empty_db(config, add_column):  # noqa: F811
    """check that we can create a dataframe from an empty database"""

    if config["mock"]:
        pytest.skip(SKIP_MSG)

    name = "tc-" + str(uuid.uuid4())[:8]

    try:
        api.create_database(name=name)

        if add_column:
            api.add_database_column(
                database_id=name,
                name="float",
                type="float",
            )

        df = DataFrame.from_deeporigin(name)
        html = df._repr_html_()

        assert "was last edited" in html, "Malformed view of empty dataframe"

        assert len(df) == 0, "Empty dataframe should have no rows"

    finally:
        # clean up
        api.delete_database(database_id=name)


def test_df_loc_indexer_2(config):  # noqa: F811
    """check that we can modify 2 rows using the loc indexer"""

    if config["mock"]:
        df = config["df"]
    else:
        df = DataFrame.from_deeporigin(config["db-name"])

    # should be able to modify an existing row
    first_row = df.index[0]
    last_row = df.index[-1]
    df.loc[[first_row, last_row]] = list(df.loc[last_row])

    assert df._modified_columns != {}, (
        "Failed to successfully modify a row using the loc indexer"
    )

    if not config["mock"]:
        df.to_deeporigin()


def test_df_loc_indexer_3(config):  # noqa: F811
    """check that we cannot add rows to a dataframe"""

    if config["mock"]:
        df = config["df"]
    else:
        df = DataFrame.from_deeporigin(config["db-name"])

    # should be able to modify an existing row
    first_row = df.index[0]

    with pytest.raises(
        DeepOriginException,
        match="Adding rows to Deep Origin DataFrames is not allowed",
    ):
        df.loc["new-row"] = list(df.loc[first_row])


@pytest.mark.parametrize("column", NUMERIC_COLUMNS)
def test_dataframe_read_modify(config, column):  # noqa: F811
    """this function tests our ability to fetch data, modify it, and write it back"""

    if config["mock"]:
        pytest.skip(SKIP_MSG)

    df = DataFrame.from_deeporigin(config["db-name"])
    df["sq_" + column] = df[column] ** 2
    df.to_deeporigin()


@pytest.mark.parametrize("column", NUMERIC_COLUMNS)
def test_dataframe_write_new_columns(config, column):  # noqa: F811
    """this function tests our ability to write new columns to a database"""

    if config["mock"]:
        pytest.skip(SKIP_MSG)

    df = DataFrame.from_deeporigin(config["db-name"])
    df["cube_" + column] = df[column] ** 3
    df.to_deeporigin()


@pytest.mark.parametrize("row", [0, -1])
def test_modify_cell(config, row):
    """test that we can modify a cell, write to it, and check that the changes stick on the remote DB

    we are parameterizing over first and last row because
    the first row has non-missing data, and the last row has missing data"""

    if config["mock"]:
        pytest.skip(SKIP_MSG)

    df = DataFrame.from_deeporigin(config["db-name"])
    row_id = df.index[row]

    # modify
    salt = np.random.random()
    df.at[row_id, "float"] = salt
    df.to_deeporigin()

    df = DataFrame.from_deeporigin(config["db-name"])
    assert df.at[row_id, "float"] == salt, (
        "Failed to successfully modify a cell using the at syntax and write to DB"
    )
