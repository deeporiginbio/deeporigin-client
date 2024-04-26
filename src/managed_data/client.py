from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Union

import requests
from beartype import beartype
from deeporigin import cache_do_api_tokens, get_do_api_tokens
from deeporigin.config import get_value
from deeporigin.do_api import read_cached_do_api_tokens
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data.schema import (
    ColumnItem,
    DescribeRowResponseDatabase,
    DescribeRowResponseRow,
    FieldItem,
    ListRowsResponse,
)
from deeporigin.utils import _nucleus_url


@dataclass
class Client(ABC):
    @abstractmethod
    def authenticate(self, refresh_tokens: bool = True):
        pass  # pragma: no cover

    @abstractmethod
    def invoke(self, endpoint: str, data: dict):
        pass  # pragma: no cover


@dataclass
class DeepOriginClient(Client):
    api_url = _nucleus_url()
    org_id = get_value()["organization_id"]

    headers = dict()

    def authenticate(self, refresh_tokens: bool = True):
        if refresh_tokens:
            api_access_token, api_refresh_token = get_do_api_tokens()
            cache_do_api_tokens(api_access_token, api_refresh_token)
        else:
            tokens = read_cached_do_api_tokens()
            api_access_token = tokens["access"]

        self.headers = {
            "accept": "application/json",
            "authorization": f"Bearer {api_access_token}",
            "content-type": "application/json",
            "x-org-id": self.org_id,
        }

    @beartype
    def invoke(
        self,
        endpoint: str,
        data: dict,
    ) -> Union[dict, list]:
        """core call to API endpoint"""

        response = requests.post(
            f"{self.api_url}{endpoint}",
            headers=self.headers,
            json=data,
        )

        return _check_response(response)


class MockClient(Client):
    """mock client to respond with static data for testing
    purposes"""

    workspaces = ["ws-placeholder"]
    databases = ["db-placeholder"]
    rows = [f"row-placeholder-{idx}" for idx in range(10)]
    columns = [f"column-{idx}" for idx in range(5)]

    def authenticate(self, refresh_tokens: bool = True):
        """no need to do anything here"""
        pass

    def invoke_list_rows(self, data: dict):
        """callback when we invoke ListRows"""
        if data == dict(filters=[]):
            # list_rows called with no filters. return
            # a workspace, a database, and some rows
            rows = []
            rows += [
                ListRowsResponse(
                    hid=self.workspaces[0],
                    id=self.workspaces[0],
                    type="workspace",
                    parentId=None,
                ).dict()
            ]
            rows += [
                ListRowsResponse(
                    hid=self.databases[0],
                    id=self.databases[0],
                    type="database",
                    parentId=self.workspaces[0],
                ).dict()
            ]
            rows += [
                ListRowsResponse(
                    hid=f"row-placeholder-{idx}",
                    id=f"row-placeholder-{idx}",
                    type="row",
                    parentId=self.databases[0],
                ).dict()
                for idx in range(10)
            ]
            return rows

        elif (
            "parent" in data["filters"][0].keys()
            and "id" in data["filters"][0]["parent"].keys()
            and data["filters"][0]["parent"]["id"].startswith("db-")
        ):
            # we're asking for rows that belong to a database

            return [
                ListRowsResponse(
                    hid=f"row-placeholder-{idx}",
                    id=f"row-placeholder-{idx}",
                    type="row",
                    parentId=data["filters"][0]["parent"]["id"],
                ).dict()
                for idx in range(10)
            ]

        elif data == dict(filters=[dict(parent=dict(isRoot=True))]) or data == dict(
            filters=[dict(rowType="workspace")]
        ):
            # return the root workspace
            return [
                ListRowsResponse(
                    hid=self.workspaces[0],
                    id=self.workspaces[0],
                    type="workspace",
                    parentId=None,
                ).dict()
            ]
        elif data == dict(filters=[dict(rowType="database")]):
            # return the example database
            return [
                ListRowsResponse(
                    hid=self.databases[0],
                    id=self.databases[0],
                    type="database",
                    parentId=self.workspaces[0],
                ).dict()
            ]

    def invoke_describe_row(self, data):
        if data["rowId"].startswith("db-"):
            # we are likely asking for a database
            name = self.databases[0]
            return DescribeRowResponseDatabase(
                id=name,
                hid=name,
                name=name,
                type="database",
                parentId="ws-placeholder",
                cols=[
                    dict(
                        ColumnItem(
                            id=f"column-{idx}",
                            key=f"column-{idx}",
                        )
                    )
                    for idx in range(5)
                ],
                hidPrefix="placeholder",
            ).dict()

        else:
            # we are asking for a row in a database
            name = data["rowId"]
            fields = [
                FieldItem(
                    columnId=f"column-{idx}",
                    cellId=f"column-{idx}",
                )
                for idx in range(5)
            ]
            row = DescribeRowResponseRow(
                id=name,
                hid=name,
                parentId=self.databases[0],
                fields=fields,
            ).dict()

            if not data["fields"]:
                row.pop("fields", None)

        return row

    def invoke_list_database_rows(self, data):
        fields = [
            FieldItem(
                columnId=f"column-{idx}",
                cellId=f"column-{idx}",
            )
            for idx in range(5)
        ]
        return [
            DescribeRowResponseRow(
                id=f"row-placeholder-{idx}",
                hid=f"row-placeholder-{idx}",
                parentId=data["databaseRowId"],
                fields=fields,
            ).dict()
            for idx in range(5)
        ]

    def invoke(self, endpoint: str, data: dict):
        """simple returns data without making any network requests"""

        if endpoint == "ListRows":
            return self.invoke_list_rows(data)

        elif endpoint == "DescribeRow":
            return self.invoke_describe_row(data)

        elif endpoint == "ConvertIdFormat":
            return [{"id": "_row:W6DjtaCrZ201EGLpmZtGO", "hid": "sample-1"}]

        elif endpoint == "ListDatabaseRows":
            return self.invoke_list_database_rows(data)

        elif endpoint == "DescribeFile":
            return file_description()

        elif endpoint == "DescribeDatabaseStats":
            return {"rowCount": 5}

        elif endpoint == "CreateFileDownloadUrl":
            return {
                "downloadUrl": "https://github.com/formiclabs/deeporigin-client/archive/refs/tags/0.0.3.zip"
            }

        elif endpoint == "ListMentions":
            return {
                "mentions": [
                    {
                        "type": "row",
                        "id": "_row:W6DjtaCrZ201EGLpmZtGO",
                        "hid": "sample-1",
                    }
                ]
            }
        elif endpoint == "ListFiles":
            if data == dict(filters=[dict(isUnassigned=True)]):
                return [{"file": file_description()}]
            elif dict(filters=[dict(isUnassigned=False)]):
                return [
                    {
                        "file": file_description(),
                        "assignments": [
                            {"rowId": "_row:ZEaEUIgsbHmGLVlgnxfvU"},
                            {"rowId": "_row:aCWxUxumDFDnu8ZhmhQ0X"},
                            {"rowId": "_row:WZVb1jsebafhfLgrHtz2l"},
                            {"rowId": "_row:3A3okCbvuaZvEkOZLqLwY"},
                        ],
                    },
                ]


@beartype
def _check_response(response: requests.models.Response) -> Union[dict, list]:
    """utility function to check responses"""

    if response.status_code == 404:
        raise DeepOriginException(
            f"[Error 404] The requested resource was not found. The response was: {response.json()}"
        )

    response.raise_for_status()
    response = response.json()

    if "error" in response:
        raise DeepOriginException(response["error"])

    if "data" in response:
        return response["data"]
    else:
        raise KeyError("`data` not in response")


def file_description():
    return dict(
        id="_file:placeholder-file",
        uri="s3://placeholder/uri",
        name="placeholder",
        status="ready",
        contentLength=123,
        contentType="application/foo",
    )
