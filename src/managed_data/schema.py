"""this module contains classes to generate data that
matches the schema of real responses"""

from datetime import datetime
from typing import Literal, Optional, Union

from pydantic import BaseModel

RowType = Literal["row", "database", "workspace"]
FileStatus = Literal["ready", "archived"]
DataType = Literal["integer", "str", "select", "date", "text", "file"]

DATAFRAME_ATTRIBUTE_KEYS = {
    "file_ids",
    "id",
    "primary_key",
    "reference_ids",
}


class DescribeFileResponse(BaseModel):
    """schema for file description, as returned by
    DescribeFile"""

    id: str
    uri: str
    name: str = "placeholder"
    status: FileStatus
    contentLength: int
    contentType: str

    class Config:
        extra = "forbid"


class ListRowsResponse(BaseModel):
    """response schema for ListRows"""

    id: str
    parentId: Optional[str]
    hid: str
    name: str = "placeholder"
    type: RowType

    class Config:
        extra = "forbid"


class FieldItem(BaseModel):
    """item in fields in DescribeRow"""

    columnId: str
    cellId: str
    validationStatus: str = "valid"
    type: DataType = "text"
    value: Union[str, dict, int] = "placeholder-text"

    class Config:
        extra = "forbid"


class ColumnItem(BaseModel):
    """item in cols in DescribeRow called on a database"""

    id: str
    name: str = "Placeholder Name"
    key: str
    parentId: str = "db-placeholder-1"
    type: DataType = "text"
    dateCreated: str = "2024-04-04T17:03:33.033115"
    cardinality: str = "one"
    canCreate: Optional[bool] = False
    configSelect: Optional[dict] = None

    class Config:
        extra = "forbid"


class DescribeRowResponse(BaseModel):
    """response schema for DescribeRow. This is complex because
    the response schema differs based on whether you call
    DescribeRow on a normal row or a database."""

    id: str
    hid: str

    parentId: str
    type: RowType = "row"
    dateCreated: datetime = "2024-04-04 16:33:58.622469"
    dateUpdated: datetime = "2024-04-04 16:33:58.622469"
    createdByUserDrn: str = "placeholder"
    rowJsonSchema: dict = {"type": "object", "required": [], "properties": {}}


class DescribeRowResponseRow(DescribeRowResponse):
    """schema for responses for DescribeRow called on a row
    this schema also works for DescribeDatabaseRows"""

    fields: Optional[list[FieldItem]] = None
    editedByUserDrn: str = "placeholder"
    hidNum: int = 1
    submissionStatus: str = "draft"
    validationStatus: str = "valid"

    class Config:
        extra = "forbid"


class DescribeRowResponseDatabase(DescribeRowResponse):
    """schema for responses for DescribeRow called on a database"""

    cols: list
    hidPrefix: str
    name: str

    class Config:
        extra = "forbid"
