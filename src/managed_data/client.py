import random
import string
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass
from typing import Union

import requests
from beartype import beartype
from deeporigin import cache_do_api_tokens, get_do_api_tokens
from deeporigin.config import get_value
from deeporigin.do_api import read_cached_do_api_tokens
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data.schema import (
    DatabaseListing,
    DatabaseRowDescription,
    RowDescription,
    RowListing,
    WorkspaceListing,
)
from deeporigin.utils import _nucleus_url


@dataclass
class Client(ABC):
    @abstractmethod
    def authenticate(self):
        pass  # pragma: no cover

    @abstractmethod
    def invoke(self):
        pass  # pragma: no cover


@dataclass
class DeepOriginClient(Client):
    api_url = _nucleus_url()
    org_id = get_value()["organization_id"]

    headers = dict()

    def authenticate(self, refresh_tokens=True):
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
    """mock client to respond with static data"""

    def authenticate():
        """no need to do anything here"""
        pass

    def invoke(self, endpoint, data):
        """overload this function so that we can return
        static data"""

        if endpoint == "ListRows":
            if data == dict(filters=[dict(parent=dict(id="db-sample"))]):
                return [
                    asdict(RowListing(hid="sample-1", id="row-1")),
                    asdict(RowListing(hid="sample-2", id="row-2")),
                ]
            elif data == dict(filters=[dict(parent=dict(isRoot=True))]) or data == dict(
                filters=[dict(rowType="workspace")]
            ):
                # return the root workspace
                return [asdict(WorkspaceListing())]
            elif data == dict(filters=[]):
                # list_rows called with no filters. return
                # a workspace, a database, and some rows
                rows = []
                rows += [asdict(WorkspaceListing(id="_row:workspace"))]
                rows += [
                    asdict(
                        DatabaseListing(parentId="_row:workspace", id="_row:database")
                    )
                ]
                rows += [
                    asdict(RowListing(parentId="_row:database", id=f"_row:{row}"))
                    for row in range(10)
                ]
                return rows

        elif endpoint == "DescribeRow":
            if data["rowId"].startswith("db-"):
                # we are likely asking for a database
                row = asdict(DatabaseRowDescription())

            else:
                # we are asking for a row in a database
                row = asdict(RowDescription())

                if not data["fields"]:
                    row.pop("fields", None)

            return row

        elif endpoint == "ConvertIdFormat":
            return [{"id": "_row:W6DjtaCrZ201EGLpmZtGO", "hid": "sample-1"}]

        elif endpoint == "ListDatabaseRows":
            row = asdict(RowDescription())
            row.pop("cols", None)
            row.pop("parent", None)
            row.pop("rowJsonSchema", None)

            return [row for _ in range(5)]

        elif endpoint == "DescribeFile":
            return file_description()

        elif endpoint == "DescribeDatabaseStats":
            return {"rowCount": 5}

        elif endpoint == "CreateFileDownloadUrl":
            return {"downloadUrl": "https://local/data"}

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
        raise DeepOriginException("[Error 404] The requested resource was not found.")

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
