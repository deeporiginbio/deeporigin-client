"""this tests low level functions in the data API"""

import os

import pytest
from deeporigin.config import get_value
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data import _api, api
from deeporigin.managed_data._api import DeepOriginClient

FILE_ID = "_file:V08GBdErNGqynC3O7bill"

# constants
row_description_keys = {
    "id",
    "parentId",
    "type",
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

list_row_keys = {"id", "hid", "type", "name", "parentId"}

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


class MockClient(DeepOriginClient):
    """mock client to respond with static data"""

    def invoke(self, endpoint, data):
        """overload this function so that we can return
        static data"""

        if endpoint == "ListRows":
            if data == dict(filters=[dict(parent=dict(id="db-sample"))]):
                return [
                    {
                        "id": "_row:0sJjiHf18ZtdzRyt1uKY5",
                        "parentId": "_row:J5FiZ1Z202GuiF78dxhMr",
                        "hid": "sample-2",
                        "type": "row",
                        "name": None,
                    },
                    {
                        "id": "_row:W6DjtaCrZ201EGLpmZtGO",
                        "parentId": "_row:J5FiZ1Z202GuiF78dxhMr",
                        "hid": "sample-1",
                        "type": "row",
                        "name": None,
                    },
                ]

        elif endpoint == "DescribeRow":
            return {
                "id": "_row:W6DjtaCrZ201EGLpmZtGO",
                "parentId": "_row:J5FiZ1Z202GuiF78dxhMr",
                "type": "row",
                "dateCreated": "2024-04-05 19:04:04.094428",
                "dateUpdated": "2024-04-08 14:58:18.376",
                "createdByUserDrn": "drn:identity::user:auth0|65ca7d6f5a130b87df994e5c",
                "editedByUserDrn": "drn:identity::user:auth0|65ca7d6f5a130b87df994e5c",
                "submissionStatus": "draft",
                "hid": "sample-1",
                "hidNum": 1,
                "validationStatus": "valid",
                "cols": [],
                "parent": {"id": "_row:J5FiZ1Z202GuiF78dxhMr"},
                "rowJsonSchema": {"type": "object", "required": [], "properties": {}},
            }

        elif endpoint == "ConvertIdFormat":
            return [{"id": "_row:W6DjtaCrZ201EGLpmZtGO", "hid": "sample-1"}]

        elif endpoint == "ListDatabaseRows":
            return [
                {
                    "id": "_row:W6DjtaCrZ201EGLpmZtGO",
                    "parentId": "_row:J5FiZ1Z202GuiF78dxhMr",
                    "type": "row",
                    "dateCreated": "2024-04-05 19:04:04.094428",
                    "dateUpdated": "2024-04-08 14:58:18.376",
                    "createdByUserDrn": "drn:identity::user:auth0|65ca7d6f5a130b87df994e5c",
                    "editedByUserDrn": "drn:identity::user:auth0|65ca7d6f5a130b87df994e5c",
                    "submissionStatus": "draft",
                    "hid": "sample-1",
                    "hidNum": 1,
                    "validationStatus": "valid",
                    "fields": [],
                },
                {
                    "id": "_row:0sJjiHf18ZtdzRyt1uKY5",
                    "parentId": "_row:J5FiZ1Z202GuiF78dxhMr",
                    "type": "row",
                    "dateCreated": "2024-04-08 15:23:26.530019",
                    "dateUpdated": "2024-04-08 15:23:40.748",
                    "createdByUserDrn": "drn:identity::user:auth0|65ca7d6f5a130b87df994e5c",
                    "editedByUserDrn": "drn:identity::user:auth0|65ca7d6f5a130b87df994e5c",
                    "submissionStatus": "draft",
                    "hid": "sample-2",
                    "hidNum": 2,
                    "validationStatus": "valid",
                    "fields": [],
                },
            ]
        elif endpoint == "DescribeFile":
            return {
                "id": "_file:V08GBdErNGqynC3O7bill",
                "uri": "s3://deeporigin-nucleus-local-uploads/files/_file:V08GBdErNGqynC3O7bill",
                "name": "pbr322_egfr (1).gb",
                "status": "ready",
                "contentLength": 9757,
                "contentType": "",
            }


mock_client = MockClient()


def test_list_rows():
    rows = _api.list_rows(parent_id="db-sample", client=mock_client)

    for row in rows:
        assert set(row.keys()) == list_row_keys
        assert row["type"] == "row"


def test_describe_row():
    """test describe_row, using mocked response"""

    row = _api.describe_row("sample-1", client=mock_client)
    assert set(row.keys()) == row_description_keys


def test_convert_id_format():
    conversions = _api.convert_id_format(hids=["sample-1"], client=mock_client)
    for conversion in conversions:
        assert {"id", "hid"} == set(conversion.keys())

    with pytest.raises(DeepOriginException, match="non-None and a list of strings"):
        _api.convert_id_format()


def test_list_database_rows():
    rows = _api.list_database_rows("db-sample")

    for row in rows:
        assert set(row.keys()).difference({"name"}) == list_database_rows_keys


# def test_download_file():
#     if being_mocked:
#         return
#     _api.download_file(FILE_ID, os.getcwd())

#     file_name = _api.describe_file(FILE_ID)["name"]

#     file_path = os.path.join(os.getcwd(), file_name)
#     assert os.path.exists(file_path)

#     # clean up
#     os.remove(file_path)


# def test_get_cell_data():
#     data = api.get_cell_data(row_id="dna-1", column_name="Base Sequence")

#     assert data == "AAA"


def test_describe_file():
    data = _api.describe_file(FILE_ID, client=mock_client)

    assert set(data.keys()) == describe_file_keys


# def test_create_file_download_url():
#     """test create_file_download_url"""

#     data = _api.create_file_download_url(FILE_ID)

#     assert "downloadUrl" in data.keys(), "Expected to find `downloadUrl` in response"

#     expected_strings = ["deeporigin", "amazonaws", "GetObject"]

#     for string in expected_strings:
#         assert (
#             string in data["downloadUrl"]
#         ), f"Expected to find {string} in returned downloadUrl"
