"""this tests low level functions in the data API"""

import pytest
import requests_mock
from deep_origin import sign_into_do_platform
from deep_origin.config import get_value
from deep_origin.data_api import low_level as ll

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

    with requests_mock.Mocker() as m:
        m.post(
            f"{API_URL}DescribeRow",
            json={"data": data},
        )
        yield m


def test_describe_row(mock_describe_rows):
    """test describe_row, using mocked response"""

    sign_into_do_platform(mock=True)

    row = ll.describe_row(ROW_NAME)
    assert set(row.keys()) == row_description_keys
