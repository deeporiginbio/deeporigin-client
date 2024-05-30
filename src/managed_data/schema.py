"""The `deeporigin.managed_data.schema` module contains Pydantic `BaseModel`s that describe
responses from Deep Origin's managed data API, and literals that
describe possible values for certain queries.

These models are used both to validate responses and to generate
mock data for testing.
"""

from typing import Literal, Optional, Union

from pydantic import BaseModel, ConfigDict

RowType = Literal["row", "database", "workspace"]
"""Type of a row"""

FileStatus = Literal["ready", "archived"]
"""Status of a file"""

DataType = Literal[
    "integer",
    "str",
    "select",
    "date",
    "text",
    "file",
    "reference",
    "editor",
    "float",
]
"""Type of a column"""

DATAFRAME_ATTRIBUTE_KEYS = {
    "file_ids",
    "id",
    "primary_key",
    "reference_ids",
}


Cardinality = Literal["one", "many"]

IDFormat = Literal["human-id", "system-id"]
"""Format of an ID"""

DatabaseReturnType = Literal["dataframe", "dict"]
"""Return type of a database"""


class CreateDatabaseResponse(BaseModel):
    """Schema for responses from the
    `CreateDatabase` endpoint."""

    id: str
    type: RowType = "database"
    name: str
    parentId: str = "placeholder"  # sometimes not returned
    dateCreated: str
    dateUpdated: str
    createdByUserDrn: str
    hid: str
    hidPrefix: str
    cols: list = []
    rowJsonSchema: dict = dict()

    model_config = ConfigDict(extra="forbid")


class CreateWorkspaceResponse(BaseModel):
    """Schema for responses from the
    `CreateWorkspace` endpoint."""

    id: str
    hid: str
    name: str
    dateCreated: str
    dateUpdated: str
    createdByUserDrn: str
    type: RowType = "workspace"
    rowJsonSchema: dict = dict()

    model_config = ConfigDict(extra="forbid")


class DescribeFileResponse(BaseModel):
    """Schema for responses from the
    `DescribeFile` endpoint."""

    id: str
    uri: str
    name: str = "placeholder"
    status: FileStatus
    contentLength: int
    contentType: str
    dateCreated: str
    dateUpdated: str
    createdByUserDrn: str

    model_config = ConfigDict(extra="forbid")


class ListRowsResponse(BaseModel):
    """Schema for responses from the
    `ListRows` endpoint."""

    id: str
    parentId: Optional[str]
    hid: str
    name: Optional[str] = "placeholder"
    type: RowType

    model_config = ConfigDict(extra="forbid")


class FieldItem(BaseModel):
    """Schema for items in the `fields` attribute of responses from
    the `DescribeRow` endpoint."""

    columnId: str
    cellId: str
    validationStatus: str = "valid"
    type: DataType = "text"
    value: Union[str, dict, int, float] = "placeholder-text"
    systemType: Optional[str] = None

    model_config = ConfigDict(extra="forbid")


class ColumnItem(BaseModel):
    """Schema for items in the `cols` attribute of responses from the `DescribeRow` endpoint for a database."""

    id: str
    name: str = "Placeholder Name"
    key: str
    parentId: str = "db-placeholder-1"
    type: DataType = "text"
    dateCreated: str = "2024-04-04T17:03:33.033115"
    cardinality: Cardinality = "one"
    canCreate: Optional[bool] = False
    configSelect: Optional[dict] = None

    model_config = ConfigDict(extra="forbid")


class AddDatabaseColumnResponse(BaseModel):
    """Schema for responses from the
    `AddDatabaseColumn` endpoint."""

    column: ColumnItem
    database: dict


class DescribeRowResponse(BaseModel):
    """Schema for responses from the `DescribeRow` endpoint. This schema is complex because
    the response schema depends on whether `DescribeRow` is called for a row or database."""

    id: str
    hid: str

    parentId: str
    type: RowType = "row"
    dateCreated: str = "2024-04-04 16:33:58.622469"
    dateUpdated: str = "2024-04-04 16:33:58.622469"
    createdByUserDrn: str = "placeholder"
    rowJsonSchema: dict = {"type": "object", "required": [], "properties": {}}


class DescribeRowResponseRow(DescribeRowResponse):
    """Schema for responses from the `DescribeRow` endpoint for a row.
    This is also the schema for responses from the `DescribeDatabaseRows` endpoint."""

    fields: Optional[list[FieldItem]] = None
    editedByUserDrn: str = "placeholder"
    hidNum: int = 1
    submissionStatus: str = "draft"
    validationStatus: str = "valid"

    model_config = ConfigDict(extra="forbid")


class DescribeRowResponseDatabase(DescribeRowResponse):
    """Schema for responses for the `DescribeRow` endpoint for a database."""

    cols: list
    hidPrefix: str
    name: str

    model_config = ConfigDict(extra="forbid")
