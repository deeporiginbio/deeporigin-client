"""this tests low level functions in the data API"""

from typing import Optional

import pytest
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data import _api, api
from deeporigin.managed_data._api import DeepOriginClient

# constants
dataframe_attr_keys = {
    "file_ids",
    "id",
    "primary_key",
    "reference_ids",
}


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


def _row_object(
    *,
    parentId: Optional[str] = "_row:placeholder",
    id: str = "row:placeholder",
    type: str = "row",
    name: Optional[str] = None,
    hid: str = "placeholder-1",
) -> dict:
    """helper function to create dummy responses"""
    return dict(
        id=id,
        parentId=parentId,
        type=type,
        name=name,
        hid=hid,
    )


def _database_row_object(
    submission_status: str = "draft",
    hid: str = "sample-1",
    hidNum: int = 1,
    validation_status: str = "valid",
):
    return {
        "id": "_row:W6DjtaCrZ201EGLpmZtGO",
        "parentId": "_row:J5FiZ1Z202GuiF78dxhMr",
        "type": "row",
        "dateCreated": "2024-04-05 19:04:04.094428",
        "dateUpdated": "2024-04-08 14:58:18.376",
        "createdByUserDrn": "drn:identity::user:auth0|65ca7d6fb87df994e5c",
        "editedByUserDrn": "drn:identity::user:auth0|65ca787df994e5c",
        "submissionStatus": submission_status,
        "hid": hid,
        "hidNum": hidNum,
        "validationStatus": validation_status,
        "fields": [],
    }


class MockClient(DeepOriginClient):
    """mock client to respond with static data"""

    def invoke(self, endpoint, data):
        """overload this function so that we can return
        static data"""

        if endpoint == "ListRows":
            if data == dict(filters=[dict(parent=dict(id="db-sample"))]):
                return [
                    _row_object(hid="sample-1"),
                    _row_object(hid="sample-2"),
                ]
            elif data == dict(filters=[dict(parent=dict(isRoot=True))]) or data == dict(
                filters=[dict(rowType="workspace")]
            ):
                return [
                    _row_object(
                        type="workspace",
                        parentId=None,
                    )
                ]

        elif endpoint == "DescribeRow":
            if data["fields"]:
                return {
                    "id": "_row:J5FiZ1Z202GuiF78dxhMr",
                    "parentId": "_row:sWOOkUQ3GEUPH3NzKt2f1",
                    "type": "database",
                    "name": "Sample",
                    "dateCreated": "2024-04-04 17:03:33.033115",
                    "dateUpdated": "2024-04-04 17:03:33.033115",
                    "createdByUserDrn": "drn:identity::user:auth0|65ca7d6f5a130b87df994e5c",
                    "hid": "db-sample",
                    "hidPrefix": "sample",
                    "cols": [
                        {
                            "id": "_col:2uMlXWSjvBaeYCiR2znkP",
                            "name": "Order ID",
                            "parentId": "_row:J5FiZ1Z202GuiF78dxhMr",
                            "type": "reference",
                            "dateCreated": "2024-04-04T17:03:33.033115",
                            "cardinality": "one",
                            "referenceDatabaseRowId": "_row:Zl8k9fwC47gf3q4846CLe",
                        },
                        {
                            "id": "_col:rFfMhRDSgiQuaY55roMIy",
                            "name": "Status",
                            "parentId": "_row:J5FiZ1Z202GuiF78dxhMr",
                            "type": "select",
                            "dateCreated": "2024-04-04T17:03:33.033115",
                            "cardinality": "one",
                            "configSelect": {
                                "options": [
                                    "Report sent to client",
                                    "Clinical interpretation completed",
                                    "Secondary analysis completed",
                                    "Primary analysis completed",
                                    "Sample processed by CRO",
                                    "Ordered",
                                ],
                                "canCreate": False,
                            },
                        },
                        {
                            "id": "_col:MYIsZA1z5K5pD7s2o4ouY",
                            "name": "To client tracking",
                            "parentId": "_row:J5FiZ1Z202GuiF78dxhMr",
                            "type": "text",
                            "dateCreated": "2024-04-04T17:03:33.033115",
                            "cardinality": "one",
                        },
                        {
                            "id": "_col:3Ge9GS0eej8cuHZ26uISb",
                            "name": "From client tracking",
                            "parentId": "_row:J5FiZ1Z202GuiF78dxhMr",
                            "type": "text",
                            "dateCreated": "2024-04-04T17:03:33.033115",
                            "cardinality": "one",
                        },
                    ],
                    "parent": {"id": "_row:sWOOkUQ3GEUPH3NzKt2f1"},
                    "fields": [],
                    "rowJsonSchema": {},
                }

            else:
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
                    "rowJsonSchema": {
                        "type": "object",
                        "required": [],
                        "properties": {},
                    },
                }

        elif endpoint == "ConvertIdFormat":
            return [{"id": "_row:W6DjtaCrZ201EGLpmZtGO", "hid": "sample-1"}]

        elif endpoint == "ListDatabaseRows":
            return [
                _database_row_object(),
                _database_row_object(),
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

        elif endpoint == "DescribeDatabaseStats":
            return {"rowCount": 5}
        elif endpoint == "ListFiles":
            if data == dict(filters=[dict(isUnassigned=True)]):
                return [
                    {
                        "file": {
                            "id": "_file:2n5jHmnbLC4tJShNrk6Df",
                            "uri": "s3://deeporigin-nucleus-local-uploads/files/_file:2n5jHmnbLC4tJShNrk6Df",
                            "name": "QC report (1).pdf",
                            "status": "archived",
                            "contentLength": 237478,
                            "contentType": "application/pdf",
                        }
                    },
                ]
            elif dict(filters=[dict(isUnassigned=False)]):
                return [
                    {
                        "file": {
                            "id": "_file:Fi7dHZJHgA3nqT1y1Ro5u",
                            "uri": "s3://deeporigin-nucleus-local-uploads/files/_file:Fi7dHZJHgA3nqT1y1Ro5u",
                            "name": "sequence_preprocessing.pdf",
                            "status": "ready",
                            "contentLength": 698255,
                            "contentType": "application/pdf",
                        },
                        "assignments": [
                            {"rowId": "_row:ZEaEUIgsbHmGLVlgnxfvU"},
                            {"rowId": "_row:aCWxUxumDFDnu8ZhmhQ0X"},
                            {"rowId": "_row:WZVb1jsebafhfLgrHtz2l"},
                            {"rowId": "_row:3A3okCbvuaZvEkOZLqLwY"},
                        ],
                    },
                ]


@pytest.fixture(scope="session", autouse=True)
def config(pytestconfig):
    """this fixture performs some setup tasks
    before all tests are run, and runs only once"""

    data = dict()

    # set up client
    if pytestconfig.getoption("client") == "mock":
        data["client"] = MockClient()

        data["databases"] = ["db-sample"]
        data["rows"] = ["sample-1"]
    else:
        client = DeepOriginClient()
        client.authenticate()
        data["client"] = client

        # if we're going to be making requests to a live
        # instance, we need to make sensible requests
        databases = _api.list_rows(row_type="database")
        data["databases"] = [db["hid"] for db in databases]
        rows = _api.list_rows(row_type="row")
        data["rows"] = [row["hid"] for row in rows]

    # tests run on yield
    yield data

    # teardown tasks, if any


def test_list_rows(config):
    rows = _api.list_rows(
        parent_id=config["databases"][0],
        client=config["client"],
    )

    for row in rows:
        assert set(row.keys()) == list_row_keys
        assert row["type"] == "row"


def test_list_rows_root_parent(config):
    root = _api.list_rows(
        parent_is_root=True,
        client=config["client"],
    )

    assert len(root) == 1, "Expected there to be exactly one root"

    assert root[0]["parentId"] is None, "Expected root to have no parent"


def test_list_rows_by_type(config):
    rows = _api.list_rows(
        row_type="workspace",
        client=config["client"],
    )

    assert len(rows) > 0, "Expected at least one workspace"

    for row in rows:
        assert (
            row["type"] == "workspace"
        ), f"Expected to get a list of workspaces, but {row} is not a workspace"


def test_list_files_unassigned(config):
    files = _api.list_files(
        is_unassigned=True,
        client=config["client"],
    )

    assert len(files) > 0, "Expected to find at least 1 unassigned files"

    for file in files:
        assert (
            "assignments" not in file.keys()
        ), f"Expected not to see an assignments key for this file, but instead found {file}"


def test_list_files_assigned(config):
    files = _api.list_files(
        is_unassigned=False,
        client=config["client"],
    )

    assert len(files) > 0, "Expected to find at least 1 assigned files"

    for file in files:
        assert (
            "assignments" in file.keys()
        ), f"Expected to see an assignments key for this file, but instead found {file}"

        assignments = file["assignments"]
        assert (
            len(assignments) > 0
        ), "Expected assignments to be a list with at least 1 element"

        for assignment in assignments:
            assert (
                "rowId" in assignment.keys()
            ), f"Expected to find a rowId in assignments, but instead found {assignment}"


def test_describe_database_stats(config):
    stats = _api.describe_database_stats(
        config["databases"][0],
        client=config["client"],
    )

    assert (
        "rowCount" in stats.keys()
    ), f"Expected stats to be a dictionary with a key called `rowCount`, instead got {stats}"

    assert stats["rowCount"] >= 0, "Expected a positive number of rows"


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
    df = api.get_dataframe(
        config["databases"][0],
        client=config["client"],
    )

    assert (
        set(df.attrs.keys()) == dataframe_attr_keys
    ), f"Expected to find a dictionary in `df.attrs` with these keys: {dataframe_attr_keys}, instead found a dictionary with these keys: {df.attrs.keys()}"

    assert (
        "Validation Status" in df.columns
    ), f"Expected to find a column called `Validation Status` in the dataframe. Instead, the columns in this dataframe are: {df.columns}"

    assert (
        df.attrs["primary_key"] in df.columns
    ), "Expected to find the primary key as a column"
