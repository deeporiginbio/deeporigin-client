"""this module contains classes to generate data that
matches the schema of real responses"""

from dataclasses import dataclass, field
from typing import Literal

RowType = Literal["row", "database"]


@dataclass
class FieldItem:
    """item in fields in DescribeRow"""

    columnId: str = "_col:Q2mLSmW4LNGDXAXCUUa1j"
    cellId: str = "_cell:a6rAF44Xfyf9nwaGfo6oG"
    validationStatus: str = "valid"
    type: str = "text"
    value: str = "12356789"


@dataclass
class ColumnItem:
    """item in cols in DescribeRow"""

    id: str = "_col:generated"
    name: str = "Order ID"
    key: str = "order_id"
    parentId: str = "_row:J5FiZ1Z202GuiF78dxhMr"
    type: str = "text"
    dateCreated: str = "2024-04-04T17:03:33.033115"
    cardinality: str = "one"


@dataclass
class GenericRowDescription:
    id: str = "_row: foo"
    parentId: str = "_row: bar"
    type: RowType = "row"
    dateCreated: str = "2024-04-04 16:34:19.968737"
    dateUpdated: str = "2024-04-04 16:34:19.968737"
    createdByUserDrn: str = "drn:identity::user:auth0|65ca1d6f5a130b87df994e5c"
    editedByUserDrn: str = "drn:identity::user:auth0|65ca1d6f5a130b87df994e5c"
    submissionStatus: str = "valid"
    hid: str = "row-1"
    hidNum: int = 1
    validationStatus: str = "valid"

    parent: dict = field(default_factory=lambda: {"wow": 412})
    rowJsonSchema: dict = field(default_factory=dict)


@dataclass
class RowDescription(GenericRowDescription):
    fields: list = field(default_factory=lambda: [FieldItem() for _ in range(5)])


@dataclass
class DatabaseRowDescription(GenericRowDescription):
    type: RowType = "database"
    hidPrefix: str = "prefix"
    cols: list = field(
        default_factory=lambda: [
            ColumnItem(id=f"_col:column-{col}") for col in range(5)
        ]
    )
