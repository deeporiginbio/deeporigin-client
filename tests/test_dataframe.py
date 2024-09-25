"""this module contains tests for the DataFrame class"""

import numpy as np
import pytest
from deeporigin.data_hub import api
from deeporigin.data_hub.dataframe import DataFrame

from tests.utils import clean_up_test_objects

TEST_PREFIX = "tc-x1zrncgf6F-"
TEST_DB_NAME = TEST_PREFIX + "db"
COLUMNS = ["float", "integer"]


@pytest.fixture(scope="session", autouse=True)
def config(pytestconfig):
    """this config block sets up a database that we can operate on in the tests"""
    data = dict()
    if pytestconfig.getoption("client") == "mock":
        data["mock"] = True
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

    # tests run on yield
    yield data

    # teardown code
    if pytestconfig.getoption("client") != "mock":
        clean_up_test_objects(TEST_PREFIX)


def test_dataframe_operations(config):
    """this function tests our ability to fetch data, modify it, and write it back"""

    if config["mock"]:
        pytest.skip(
            "Skipping this test because we can't test this well using a mocked client"
        )

    # we are explicitly using a loop in here rather than
    # pytest.parameterize because these tests run in parallel.
    # this can cause race conditions where the teardown/setup code
    # can be called simultaneously from two workers
    for column in COLUMNS:
        df = DataFrame.from_deeporigin(TEST_DB_NAME)
        df["sq_" + column] = df[column] ** 2
        df.to_deeporigin()
