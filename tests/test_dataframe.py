"""this module contains tests for the DataFrame class"""

import uuid

import numpy as np
import pandas as pd
import pytest
from deeporigin.data_hub import api
from deeporigin.data_hub.dataframe import DataFrame

from tests.utils import clean_up_test_objects

TEST_PREFIX = "tc-"

COLUMNS = ["float", "integer"]

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

    else:
        data["mock"] = False

        rows = api.list_rows(row_type="database")
        if TEST_DB_NAME in [row.hid for row in rows]:
            api.delete_database(database_id=TEST_DB_NAME)

        api.create_database(name=TEST_DB_NAME)

        # make columns for values and squared values
        for column in COLUMNS:
            api.add_database_column(
                database_id=TEST_DB_NAME,
                name=column,
                type=column,
                key=column,
            )

            api.add_database_column(
                database_id=TEST_DB_NAME,
                name="sq_" + column,
                type=column,
                key="sq_" + column,
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
        clean_up_test_objects(TEST_PREFIX + salt)


@pytest.mark.parametrize("column", COLUMNS)
def test_dataframe_read_modify(config, column):  # noqa: F811
    """this function tests our ability to fetch data, modify it, and write it back"""

    if config["mock"]:
        pytest.skip(SKIP_MSG)

    df = DataFrame.from_deeporigin(config["db-name"])
    df["sq_" + column] = df[column] ** 2
    df.to_deeporigin()


@pytest.mark.parametrize("column", COLUMNS)
def test_dataframe_write_new_columns(config, column):  # noqa: F811
    """this function tests our ability to write new columns to a database"""

    if config["mock"]:
        pytest.skip(SKIP_MSG)

    df = DataFrame.from_deeporigin(config["db-name"])
    df["cube_" + column] = df[column] ** 3
    df.to_deeporigin()


def test_slicing_restrictions(config):
    """check that we control what happens when a dataframe is sliced by rows"""

    if config["mock"]:
        df = DataFrame(config["df"])
    else:
        df = DataFrame.from_deeporigin(config["db-name"])

    assert df._allow_adding_rows is True, "Expected _allow_adding_rows to be True"

    # slice to 2 rows
    df = df.loc[df.index[:2]]

    assert df._allow_adding_rows is False, "Expected _allow_adding_rows to be False"


def test_slice_and_modify(config):
    """check that we control what happens when a dataframe is sliced by rows"""

    if config["mock"]:
        df = DataFrame(config["df"])
    else:
        df = DataFrame.from_deeporigin(config["db-name"])

    assert df._allow_adding_rows is True, "Expected _allow_adding_rows to be True"

    # slice to 2 rows
    df = df.loc[df.index[:2]]

    assert df._allow_adding_rows is False, "Expected _allow_adding_rows to be False"

    # we should be allowed to modify a slice
    df.iloc[0]["integer"] = 100

    if config["mock"]:
        return

    # should be allowed to write back to DB
    df.to_deeporigin()


def test_slice_and_extend_loc(config):
    """test that we can add a row to a whole df, but not to a slice of it"""

    if config["mock"]:
        df = DataFrame(config["df"])
        row_prefix = "x"
    else:
        df = DataFrame.from_deeporigin(config["db-name"])
        row_prefix = df.attrs["metadata"]["hid_prefix"]

    # should be possible to add a new row
    df.loc[row_prefix + "-" + str(len(df) + 1)] = list(df.loc[df.index[0]])

    if not config["mock"]:
        df.to_deeporigin()

    # slice to 2 rows
    df = df.loc[df.index[:2]]

    # should not be possible to add a new row
    with pytest.raises(ValueError, match="Adding rows is not allowed"):
        df.loc[row_prefix + "-" + str(len(df) + 1)] = list(df.loc[df.index[0]])
