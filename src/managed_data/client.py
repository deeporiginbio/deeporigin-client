import json
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Union
from urllib.parse import urljoin

import requests
from beartype import beartype
from deeporigin import auth
from deeporigin.config import get_value
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data.schema import (
    ColumnItem,
    CreateDatabaseResponse,
    CreateWorkspaceResponse,
    DescribeFileResponse,
    DescribeRowResponseDatabase,
    DescribeRowResponseRow,
    FieldItem,
    ListRowsResponse,
)
from deeporigin.utils import _nucleus_url


def file_description():
    """mock file description"""

    return DescribeFileResponse(
        id="_file:placeholder-file",
        uri="s3://placeholder/uri",
        name="placeholder",
        status="ready",
        contentLength=123,
        contentType="application/foo",
        dateCreated="2024-05-08 19:16:09.26078",
        dateUpdated="2024-05-08 19:16:09.26078",
        createdByUserDrn="drn:identity::user:google-apps|foo@deeporigin.com",
    ).model_dump()


@dataclass
class Client(ABC):
    """client abstract base class. Actual clients inherit from this base class. Actual clients need to implement a `authenticate` and `invoke` method"""

    @abstractmethod
    def authenticate(self, refresh_tokens: bool = True):
        """authenticate method to be called before making any requests. This needs to be implemented in derived classes"""

        pass  # pragma: no cover

    @abstractmethod
    def invoke(self, endpoint: str, data: dict):
        """invoke method to be called to make requests. This needs to be implemented in derived classes"""

        pass  # pragma: no cover

    def __props__(self):
        """helper method to show all properties of the client. Do not use."""
        props = dict()
        for attribute in dir(self):
            if attribute.startswith("_"):
                continue

            try:
                if not callable(getattr(self, attribute)):
                    props[str(attribute)] = getattr(self, attribute)
            except Exception:
                pass

        return json.dumps(props, indent=4)

    def __repr__(self):
        """helper method to show all properties of the client. Do not use."""
        return self.__props__()

    def _repr_html_(self):
        """helper method to show all properties of the client. Do not use. This method is called when you invoke a variable in a Jupyter instance"""
        return self.__props__()

    def __str__(self):
        """helper method to show all properties of the client. Do not use."""
        return self.__props__()


@dataclass
class DeepOriginClient(Client):
    """client to interact with DeepOrigin API"""

    api_url = _nucleus_url()
    org_id = get_value()["organization_id"]

    headers = dict()

    def authenticate(self, refresh_tokens: bool = True):
        """authenticate to DeepOrigin API

        Authenticate to Deep Origin API. This needs to be called
        before making any requests. refresh_tokens=False allows
        for turning off refresh (saving one network request on
        startup, ), making the CLI faster

        Args:
            refresh_tokens (bool, optional): Whether to
            refresh tokens. Defaults to True.


        """

        tokens = auth.get_tokens(refresh=refresh_tokens)

        self.headers = {
            "accept": "application/json",
            "authorization": f"Bearer {tokens['access']}",
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
            urljoin(self.api_url, endpoint),
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
    file = file_description()

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
                ).model_dump()
            ]
            rows += [
                ListRowsResponse(
                    hid=self.databases[0],
                    id=self.databases[0],
                    type="database",
                    parentId=self.workspaces[0],
                ).model_dump()
            ]
            rows += [
                ListRowsResponse(
                    hid=f"row-placeholder-{idx}",
                    id=f"row-placeholder-{idx}",
                    type="row",
                    parentId=self.databases[0],
                ).model_dump()
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
                ).model_dump()
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
                ).model_dump()
            ]
        elif data == dict(filters=[dict(rowType="database")]):
            # return the example database
            return [
                ListRowsResponse(
                    hid=self.databases[0],
                    id=self.databases[0],
                    type="database",
                    parentId=self.workspaces[0],
                ).model_dump()
            ]
        elif data == dict(filters=[dict(rowType="row")]):
            # return the example row
            return [
                ListRowsResponse(
                    hid=self.rows[0],
                    id=self.rows[0],
                    type="row",
                    parentId=self.databases[0],
                ).model_dump()
            ]
        else:
            print(data)
            raise NotImplementedError(
                "This specific request type hasn't been mocked yet"
            )

    def invoke_describe_row(self, data):
        """callback when we invoke DescribeRow"""

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
            ).model_dump()

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
            ).model_dump()

            if not data["fields"]:
                row.pop("fields", None)

        return row

    def invoke_list_database_rows(self, data):
        """callback when we invoke ListDatabaseRows"""

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
            ).model_dump()
            for idx in range(5)
        ]

    def invoke(self, endpoint: str, data: dict):
        """simple returns data without making any network requests"""

        if endpoint == "ListRows":
            return self.invoke_list_rows(data)

        elif endpoint == "CreateWorkspace":
            name = data["workspace"]["name"]
            return CreateWorkspaceResponse(
                id="_workspace:" + name,
                hid="_workspace:" + name,
                name=name,
                dateCreated="2024-05-08 19:16:09.26078",
                dateUpdated="2024-05-08 19:16:09.26078",
                createdByUserDrn="drn:identity::user:google-apps|foo@deeporigin.com",
            ).model_dump()

        elif endpoint == "DeleteRows":
            return dict(data=dict())

        elif endpoint == "CreateDatabase":
            name = data["database"]["name"]
            return CreateDatabaseResponse(
                id="_database:" + name,
                hid="_database:" + name,
                hidPrefix="placeholder",
                name=name,
                dateCreated="2024-05-08 19:16:09.26078",
                dateUpdated="2024-05-08 19:16:09.26078",
                createdByUserDrn="drn:identity::user:google-apps|foo@deeporigin.com",
                parentId=data["database"]["parentId"],
                cols=[],
            ).model_dump()

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
    elif response.status_code == 400:
        raise DeepOriginException(f"[Error 400] The response was: {response.json()}")

    response.raise_for_status()
    response = response.json()

    if "error" in response:
        raise DeepOriginException(response["error"])

    if "data" in response:
        return response["data"]
    else:
        raise KeyError("`data` not in response")
