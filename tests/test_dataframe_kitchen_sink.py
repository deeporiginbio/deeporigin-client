"""this test tests against a live instance which contains a kitchen sink database"""

import random

import numpy as np
import pytest

from deeporigin.data_hub import api
from deeporigin.data_hub.dataframe import DataFrame
from tests.utils import minimal_config as config  # noqa: F401


def skip_if_unavailable(cfg):
    if cfg["mock"]:
        pytest.skip(
            "Skipping this test because we can't test this well using a mocked client"
        )

    dbs = api.list_rows(row_type="database")
    db = [db for db in dbs if db.hid == "kitchen-sink"]
    if len(db) != 1:
        pytest.skip(
            "Skipping this test because there is no DB called `kitchen-sink` on Deep Origin"
        )


def test_kitchen_sink_db_int(config):  # noqa: F811
    """this test assumes that there is a DB called `kitchen-sink` on Deep Origin, and that we're running against a live instance"""

    skip_if_unavailable(config)

    # get the DB
    df = DataFrame.from_deeporigin("kitchen-sink")

    # test writing a single value
    df.at["ks-34", "Int"] = np.random.randint(len(df))
    df.to_deeporigin()

    # test writing entire columns
    df["Int"] = np.random.randint(50, 999, len(df))
    df.to_deeporigin()


def test_kitchen_sink_db_float(config):  # noqa: F811
    """this test assumes that there is a DB called `kitchen-sink` on Deep Origin, and that we're running against a live instance"""

    skip_if_unavailable(config)

    # get the DB
    df = DataFrame.from_deeporigin("kitchen-sink")

    # test writing a single value
    df.at["ks-34", "Float"] = np.random.random()
    df.to_deeporigin()

    # test writing entire columns
    df["Float"] = np.random.random(len(df)) * 100
    df.to_deeporigin()


def test_kitchen_sink_db_bool(config):  # noqa: F811
    """this test assumes that there is a DB called `kitchen-sink` on Deep Origin, and that we're running against a live instance"""

    skip_if_unavailable(config)

    # get the DB
    df = DataFrame.from_deeporigin("kitchen-sink")

    # test writing a single value
    df.at["ks-34", "Bool"] = True
    df.to_deeporigin()

    # test writing entire columns
    df["Bool"] = [random.choice([True, False]) for _ in range(len(df))]

    df.to_deeporigin()
