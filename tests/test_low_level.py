"""this tests low level functions in the data API"""

import sys

import pytest
import requests_mock
from deeporigin import read_cached_do_api_tokens
from deeporigin.config import get_value
from deeporigin.managed_data import _api

# constants
row_description_keys = {
    "id",
    "parentId",
    "type",
    "name",
    "dateCreated",
    "dateUpdated",
    "createdByUserDrn",
    "editedByUserDrn",
    "submissionStatus",
    "hid",
    "hidNum",
    "validationStatus",
    "cols",
    "parent",
    "rowJsonSchema",
}


API_URL = get_value()["nucleus_api_endpoint"]
DB_NAME = "db-dna"
ROW_NAME = "dna-1"
AUTH_DOMAIN = get_value()["auth_domain"]


@pytest.fixture
def mocker(pytestconfig):
    """mock endpoints of nucleus service"""

    no_mock = pytestconfig.getoption("nomock")

    data = dict()
    for key in row_description_keys:
        data[key] = "placeholder"

    tokens = read_cached_do_api_tokens()

    if no_mock:
        print("NO MOCKING! RuNNING LIVe")
        with requests_mock.Mocker(real_http=True) as m:
            yield m
    else:
        print("Mocking..")

        with requests_mock.Mocker() as m:
            m.post(
                f"{API_URL}DescribeRow",
                json={"data": data},
            )
            m.post(
                f"{AUTH_DOMAIN}/oauth/token",
                json=dict(access_token=tokens["access"]),
            )

            # file download
            m.post(
                f"{API_URL}CreateFileDownloadUrl",
                json=dict(
                    data=dict(downloadUrl="https://deeporigin-amazonaws-GetObject")
                ),
            )

            yield m


def test_describe_row(mocker):
    """test describe_row, using mocked response"""

    row = _api.describe_row(ROW_NAME)
    assert set(row.keys()) == row_description_keys


def test_create_file_download_url(mocker):
    """test create_file_download_url"""

    data = _api.create_file_download_url("_file:4Jii3X9t7PF3IcG6SDmo4")

    assert "downloadUrl" in data.keys(), "Expected to find `downloadUrl` in response"

    expected_strings = ["deeporigin", "amazonaws", "GetObject"]

    for string in expected_strings:
        assert (
            string in data["downloadUrl"]
        ), f"Expected to find {string} in returned downloadUrl"
