"""this tests low level functions in the data API"""


import pytest
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data import _api, api
from deeporigin.managed_data._api import DeepOriginClient

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


@pytest.fixture(scope="session", autouse=True)
def config(pytestconfig):
    """this fixture performs some setup tasks
    before all tests are run, and runs only once"""

    config = dict()

    # set up client
    if pytestconfig.getoption("client") == "mock":
        config["client"] = MockClient()

        config["databases"] = ["db-sample"]
        config["rows"] = ["sample-1"]
    else:
        config["client"] = DeepOriginClient()

        # if we're going to be making requests to a live
        # instance, we need to make sensible requests
        databases = _api.list_rows(row_type="database")
        config["databases"] = [db["hid"] for db in databases]
        rows = _api.list_rows(row_type="row")
        config["rows"] = [row["hid"] for row in rows]

    # tests run on yield
    yield config

    # teardown tasks, if any


def test_list_rows(config):
    rows = _api.list_rows(
        parent_id=config["databases"][0],
        client=config["client"],
    )

    for row in rows:
        assert set(row.keys()) == list_row_keys
        assert row["type"] == "row"


def test_describe_row(config):
    """test describe_row, using mocked response"""

    row = _api.describe_row(
        config["rows"][0],
        client=config["client"],
    )
    assert set(row.keys()) == row_description_keys


def test_convert_id_format(config):
    conversions = _api.convert_id_format(
        hids=[config["rows"][0]],
        client=config["client"],
    )
    for conversion in conversions:
        assert {"id", "hid"} == set(conversion.keys())

    with pytest.raises(DeepOriginException, match="non-None and a list of strings"):
        _api.convert_id_format()


def test_list_database_rows(config):
    rows = _api.list_database_rows(
        config["databases"][0],
        client=config["client"],
    )

    for row in rows:
        assert set(row.keys()).difference({"name"}) == list_database_rows_keys


def test_get_dataframe(config):
    df = api.get_dataframe(config["databases"][0])

    assert (
        set(df.attrs.keys()) == {"file_ids", "reference_ids"}
    ), "Expected to find a dictionary in `df.attrs` with keys called `file_ids` and `reference_ids`"

    assert (
        "validation_status" in df.columns
    ), f"Expected to find a column called `validaton_status` in the dataframe. Instead, the columns in this dataframe are: {df.columns}"

    assert df.index.name == "row", "Expected the index to be called `row`"
