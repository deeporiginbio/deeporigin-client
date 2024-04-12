"""this tests low level functions in the data API"""

import os
from urllib.parse import urljoin

import pytest
import requests_mock
from deeporigin import read_cached_do_api_tokens
from deeporigin.config import get_value
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data import _api, api
from deeporigin.utils import _nucleus_url

API_URL = _nucleus_url()
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


describe_file_keys = {
    "id",
    "uri",
    "name",
    "status",
    "contentLength",
    "contentType",
}

list_row_keys = {"id", "hid", "type", "name"}

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

MOCK_URL = "https://deeporigin.mock/nucleus-api/api/"
DB_NAME = "db-dna"
ROW_NAME = "dna-1"
AUTH_DOMAIN = get_value()["auth_domain"]

# if we're running against a real instance, determine
# if the database has any files in it
if API_URL != MOCK_URL:
    df = api.get_dataframe(DB_NAME)
    file_ids = df.attrs["file_ids"]
    if len(file_ids) == 0:
        raise RuntimeError("No files in database, cannot run tests")
    FILE_ID = file_ids[0]
    being_mocked = False
    print("Testing against live instance")
else:
    FILE_ID = "_file:placeholder"
    being_mocked = True
    print("Using mocked responses")


@pytest.fixture
def mocker():
    """mock endpoints of nucleus service. This mocker
    mocks every request going to `MOCK_URL` but passes
    through real requests going to a live instance. Thus,
    the same mocker can handle fully mocked and fully
    live tests"""

    data = dict()
    for key in row_description_keys:
        data[key] = "placeholder"

    tokens = read_cached_do_api_tokens()

    with requests_mock.Mocker(real_http=True) as m:
        m.post(
            urljoin(MOCK_URL, "DescribeRow"),
            json={"data": data},
        )
        m.post(
            urljoin(AUTH_DOMAIN, "/oauth/token/"),
            json=dict(access_token=tokens["access"]),
        )

        # list database rows
        row = dict()
        for key in list_database_rows_keys:
            row[key] = "placeholder"
        data = [row for _ in range(10)]
        m.post(urljoin(MOCK_URL, "ListDatabaseRows"), json=dict(data=data))

        # convert ID formats
        data = [{"id": "_row:3352QeWzQab5nxmvLQwvo", "hid": "db-dna"}]
        m.post(urljoin(MOCK_URL, "ConvertIdFormat"), json=dict(data=data))

        # list rows
        data = [
            {
                "id": "_row:ZP9sCgQHESaMf9jFUP3VC",
                "hid": "dna-1",
                "type": "row",
                "name": "sdsdf",
            },
            {
                "id": "_row:Ovhr6Kg9Ksczb10TwmBJi",
                "hid": "dna-3",
                "type": "row",
                "name": None,
            },
            {
                "id": "_row:92wgjsSY7jM74TmmmBTko",
                "hid": "dna-2",
                "type": "row",
                "name": "fddsf",
            },
        ]
        m.post(urljoin(MOCK_URL, "ListRows"), json=dict(data=data))

        # describe file
        data = {
            "id": "_file:4Jii3X9t7PF3IcG6SDmo4",
            "uri": "s3://deeporigin-nucleus-local-uploads/files/_file:4Jii3X9t7PF3IcG6SDmo4",
            "name": "test.png",
            "status": "ready",
            "contentLength": 554588,
            "contentType": "image/png",
        }
        m.post(urljoin(MOCK_URL, "DescribeFile"), json=dict(data=data))

        # file download
        m.post(
            urljoin(MOCK_URL, "CreateFileDownloadUrl"),
            json=dict(data=dict(downloadUrl="https://deeporigin-amazonaws-GetObject")),
        )

        yield m


def test_describe_file(mocker):
    data = _api.describe_file(FILE_ID)

    assert set(data.keys()) == describe_file_keys


def test_describe_row(mocker):
    """test describe_row, using mocked response"""

    row = _api.describe_row(ROW_NAME)
    assert set(row.keys()) == row_description_keys


def test_create_file_download_url(mocker):
    """test create_file_download_url"""

    data = _api.create_file_download_url(FILE_ID)

    assert "downloadUrl" in data.keys(), "Expected to find `downloadUrl` in response"

    expected_strings = ["deeporigin", "amazonaws", "GetObject"]

    for string in expected_strings:
        assert (
            string in data["downloadUrl"]
        ), f"Expected to find {string} in returned downloadUrl"


def test_list_rows(mocker):
    rows = _api.list_rows(DB_NAME)

    for row in rows:
        assert set(row.keys()) == list_row_keys
        assert row["type"] == "row"


def test_list_database_rows(mocker):
    rows = _api.list_database_rows(DB_NAME)

    for row in rows:
        assert set(row.keys()).difference({"name"}) == list_database_rows_keys


def test_convert_id_format(mocker):
    conversions = _api.convert_id_format(hids=[DB_NAME])
    for conversion in conversions:
        assert {"id", "hid"} == set(conversion.keys())

    with pytest.raises(DeepOriginException, match="non-None and a list of strings"):
        _api.convert_id_format()


def test_download_file(mocker):
    if being_mocked:
        return
    _api.download_file(FILE_ID, os.getcwd())

    file_name = _api.describe_file(FILE_ID)["name"]

    file_path = os.path.join(os.getcwd(), file_name)
    assert os.path.exists(file_path)

    # clean up
    os.remove(file_path)


def test_test_get_cell_data(mocker):
    if being_mocked:
        return
    data = api.get_cell_data(row_id="dna-1", column_name="Base Sequence")

    assert data == "AAA"
