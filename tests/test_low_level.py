"""this tests low level functions in the data API"""

import sys

import pytest
import requests_mock
from deeporigin import read_cached_do_api_tokens
from deeporigin.config import get_value
from deeporigin.managed_data import _api as api

# constants
row_description_keys = {
    "id",
    "parentId",
    "type",
    "name",
    "dateCreated",
    "dateUpdated",
    "editedByUserDrn",
    "hid",
    "validationStatus",
    "parent",
}


API_URL = get_value()["nucleus_api_endpoint"]
DB_NAME = "db-dna"
ROW_NAME = "dna-1"


@pytest.fixture
def mock_describe_rows():
    """mock endpoint to mimic what ListRows would return"""

    data = dict()
    for key in row_description_keys:
        data[key] = "placeholder"

    tokens = read_cached_do_api_tokens()

    with requests_mock.Mocker() as m:
        m.post(
            f"{API_URL}DescribeRow",
            json={"data": data},
        )
        m.post(
            "https://formicbio-dev.us.auth0.com/oauth/token",
            json=dict(access_token=tokens["access"]),
        )
        yield m


def test_describe_row(mock_describe_rows):
    """test describe_row, using mocked response"""

    print(sys.executable)

    row = api.describe_row(ROW_NAME)
    assert set(row.keys()) == row_description_keys
