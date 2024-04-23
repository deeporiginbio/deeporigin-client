"""this module contains classes to generate data that
matches the schema of real responses"""

from abc import ABC
from dataclasses import dataclass, field
from typing import Literal, Optional

from pydantic import BaseModel

RowType = Literal["row", "database", "workspace"]
FileStatus = Literal["ready", "archived"]


class DescribeFileResponse(BaseModel):
    """schema for file description, as returned by
    DescribeFile"""

    id: str
    uri: str
    name: str
    status: FileStatus
    contentLength: int
    contentType: str

    class Config:
        extra = "forbid"


@dataclass
class GenericRowListing(ABC):
    """abstract class for a generic row listing"""

    name: str
    id: str = "_row:placeholder"
    parentId: Optional[str] = None


@dataclass
class WorkspaceListing(GenericRowListing):
    """response for a workspace listing"""

    hid: str = "workspace"
    type: RowType = "workspace"
    name: str = "Placeholder Workspace"


@dataclass
class DatabaseListing(GenericRowListing):
    """response for a database listing"""

    hid: str = "database"
    type: RowType = "database"
    type: RowType = "database"
    name: str = "Placeholder DB"


@dataclass
class RowListing(GenericRowListing):
    """response for a true row (leaf node) listing"""

    hid: str = "row"
    type: RowType = "row"
    name: str = "Placeholder row"


@dataclass
class FieldItem:
    """item in fields in DescribeRow"""

    columnId: str = "_col:placeholder"
    cellId: str = "_cell:placeholder"
    validationStatus: str = "valid"
    type: str = "text"
    value: str = "placeholder-text"


@dataclass
class ColumnItem:
    """item in cols in DescribeRow"""

    id: str = "_col:placeholder"
    name: str = "Placeholder Name"
    key: str = "placeholder_key"
    parentId: str = "_row:placeholder_row"
    type: str = "text"
    dateCreated: str = "2024-04-04T17:03:33.033115"
    cardinality: str = "one"


@dataclass
class GenericRowDescription(ABC):
    id: str = "row-1"
    parentId: str = "_row:placeholder-parent"
    type: RowType = "row"
    dateCreated: str = "2024-04-04 16:34:19.968737"
    dateUpdated: str = "2024-04-04 16:34:19.968737"
    createdByUserDrn: str = "drn:identity::user:auth0|65ca1d6f5a130b87df994e5c"
    editedByUserDrn: str = "drn:identity::user:auth0|65ca1d6f5a130b87df994e5c"
    submissionStatus: str = "valid"
    hid: str = "row-1"
    hidNum: int = 1
    validationStatus: str = "valid"

    parent: dict = field(default_factory=lambda: {"placeholder": 42})
    rowJsonSchema: dict = field(default_factory=dict)


@dataclass
class RowDescription(GenericRowDescription):
    fields: list = field(default_factory=lambda: [FieldItem() for _ in range(5)])


@dataclass
class DatabaseRowDescription(GenericRowDescription):
    type: RowType = "database"
    hidPrefix: str = "placeholder"
    cols: list = field(
        default_factory=lambda: [
            ColumnItem(id=f"_col:column-{col}") for col in range(5)
        ]
    )
