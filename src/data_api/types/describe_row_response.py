# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from typing import List, Union, Optional
from typing_extensions import Literal, Annotated

from pydantic import Field as FieldInfo

from .._utils import PropertyInfo
from .._models import BaseModel
from .database import Database
from .workspace import Workspace

__all__ = [
    "DescribeRowResponse",
    "Data",
    "DataDatabaseRow",
    "DataDatabaseRowField",
    "DataDatabaseRowFieldUnionMember0",
    "DataDatabaseRowFieldUnionMember0InvalidData",
    "DataDatabaseRowFieldUnionMember1",
    "DataDatabaseRowFieldUnionMember1InvalidData",
    "DataDatabaseRowFieldUnionMember2",
    "DataDatabaseRowFieldUnionMember2InvalidData",
    "DataDatabaseRowFieldUnionMember3",
    "DataDatabaseRowFieldUnionMember3InvalidData",
    "DataDatabaseRowFieldUnionMember4",
    "DataDatabaseRowFieldUnionMember4InvalidData",
    "DataDatabaseRowFieldUnionMember4Value",
    "DataDatabaseRowFieldUnionMember5",
    "DataDatabaseRowFieldUnionMember5InvalidData",
    "DataDatabaseRowFieldUnionMember5Value",
    "DataDatabaseRowFieldUnionMember6",
    "DataDatabaseRowFieldUnionMember6InvalidData",
    "DataDatabaseRowFieldUnionMember6Value",
    "DataDatabaseRowFieldUnionMember7",
    "DataDatabaseRowFieldUnionMember7InvalidData",
    "DataDatabaseRowFieldUnionMember7Value",
    "DataDatabaseRowFieldUnionMember8",
    "DataDatabaseRowFieldUnionMember8InvalidData",
    "DataDatabaseRowFieldUnionMember9",
    "DataDatabaseRowFieldUnionMember9InvalidData",
    "DataDatabaseRowFieldUnionMember9Value",
    "DataDatabaseRowFieldUnionMember9ValueURL",
    "DataDatabaseRowFieldUnionMember10",
    "DataDatabaseRowFieldUnionMember10InvalidData",
    "DataDatabaseRowFieldUnionMember10Value",
    "DataDatabaseRowFieldUnionMember11",
    "DataDatabaseRowFieldUnionMember11Value",
    "DataDatabaseRowFieldUnionMember11InvalidData",
    "DataDatabaseRowFieldUnionMember12",
    "DataDatabaseRowFieldUnionMember12InvalidData",
    "DataDatabaseRowFieldUnionMember12Value",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember0",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember0InvalidData",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember1",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember1InvalidData",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember2",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember2InvalidData",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember3",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember3InvalidData",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember4",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember4InvalidData",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember4Value",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember5",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember5InvalidData",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember5Value",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember6",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember6InvalidData",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember6Value",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember7",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember7InvalidData",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember7Value",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember8",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember8InvalidData",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember9",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember9InvalidData",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember9Value",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember9ValueURL",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember10",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember10InvalidData",
    "DataDatabaseRowFieldUnionMember12ValueUnionMember10Value",
]


class DataDatabaseRowFieldUnionMember0InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember0(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["text"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember0InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[str] = None


class DataDatabaseRowFieldUnionMember1InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember1(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["integer"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember1InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[float] = None


class DataDatabaseRowFieldUnionMember2InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember2(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["float"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember2InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[float] = None


class DataDatabaseRowFieldUnionMember3InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember3(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["boolean"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember3InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[bool] = None


class DataDatabaseRowFieldUnionMember4InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember4Value(BaseModel):
    row_ids: List[str] = FieldInfo(alias="rowIds")


class DataDatabaseRowFieldUnionMember4(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["reference"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember4InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember4Value] = None


class DataDatabaseRowFieldUnionMember5InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember5Value(BaseModel):
    top_level_blocks: List[object] = FieldInfo(alias="topLevelBlocks")


class DataDatabaseRowFieldUnionMember5(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["editor"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember5InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember5Value] = None


class DataDatabaseRowFieldUnionMember6InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember6Value(BaseModel):
    file_ids: List[str] = FieldInfo(alias="fileIds")


class DataDatabaseRowFieldUnionMember6(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["file"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember6InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember6Value] = None


class DataDatabaseRowFieldUnionMember7InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember7Value(BaseModel):
    selected_options: List[str] = FieldInfo(alias="selectedOptions")


class DataDatabaseRowFieldUnionMember7(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["select"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember7InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember7Value] = None


class DataDatabaseRowFieldUnionMember8InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember8(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["date"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember8InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[str] = None


class DataDatabaseRowFieldUnionMember9InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember9ValueURL(BaseModel):
    url: str

    title: Optional[str] = None


class DataDatabaseRowFieldUnionMember9Value(BaseModel):
    urls: List[DataDatabaseRowFieldUnionMember9ValueURL]


class DataDatabaseRowFieldUnionMember9(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["url"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember9InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember9Value] = None


class DataDatabaseRowFieldUnionMember10InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember10Value(BaseModel):
    user_drns: List[str] = FieldInfo(alias="userDrns")


class DataDatabaseRowFieldUnionMember10(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["user"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember10InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember10Value] = None


class DataDatabaseRowFieldUnionMember11Value(BaseModel):
    error_message: Optional[str] = FieldInfo(alias="errorMessage", default=None)
    """The error message from executing the expression, if one occurred."""

    invalid_result: Optional[object] = FieldInfo(alias="invalidResult", default=None)
    """Expression result that is not a valid return value type."""

    result: Union[str, float, None] = None
    """The return value from executing the expression."""


class DataDatabaseRowFieldUnionMember11InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember11(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["expression"]

    validation_status: str = FieldInfo(alias="validationStatus")

    value: DataDatabaseRowFieldUnionMember11Value

    invalid_data: Optional[DataDatabaseRowFieldUnionMember11InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)


class DataDatabaseRowFieldUnionMember12InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember0InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember0(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["text"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember0InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember1InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember1(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["integer"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember1InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[float] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember2InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember2(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["float"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember2InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[float] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember3InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember3(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["boolean"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember3InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[bool] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember4InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember4Value(BaseModel):
    row_ids: List[str] = FieldInfo(alias="rowIds")


class DataDatabaseRowFieldUnionMember12ValueUnionMember4(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["reference"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember4InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember4Value] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember5InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember5Value(BaseModel):
    top_level_blocks: List[object] = FieldInfo(alias="topLevelBlocks")


class DataDatabaseRowFieldUnionMember12ValueUnionMember5(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["editor"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember5InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember5Value] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember6InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember6Value(BaseModel):
    file_ids: List[str] = FieldInfo(alias="fileIds")


class DataDatabaseRowFieldUnionMember12ValueUnionMember6(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["file"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember6InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember6Value] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember7InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember7Value(BaseModel):
    selected_options: List[str] = FieldInfo(alias="selectedOptions")


class DataDatabaseRowFieldUnionMember12ValueUnionMember7(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["select"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember7InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember7Value] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember8InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember8(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["date"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember8InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember9InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember9ValueURL(BaseModel):
    url: str

    title: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember9Value(BaseModel):
    urls: List[DataDatabaseRowFieldUnionMember12ValueUnionMember9ValueURL]


class DataDatabaseRowFieldUnionMember12ValueUnionMember9(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["url"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember9InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember9Value] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember10InvalidData(BaseModel):
    message: Optional[str] = None


class DataDatabaseRowFieldUnionMember12ValueUnionMember10Value(BaseModel):
    user_drns: List[str] = FieldInfo(alias="userDrns")


class DataDatabaseRowFieldUnionMember12ValueUnionMember10(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["user"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember10InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember12ValueUnionMember10Value] = None


DataDatabaseRowFieldUnionMember12Value = Union[
    DataDatabaseRowFieldUnionMember12ValueUnionMember0,
    DataDatabaseRowFieldUnionMember12ValueUnionMember1,
    DataDatabaseRowFieldUnionMember12ValueUnionMember2,
    DataDatabaseRowFieldUnionMember12ValueUnionMember3,
    DataDatabaseRowFieldUnionMember12ValueUnionMember4,
    DataDatabaseRowFieldUnionMember12ValueUnionMember5,
    DataDatabaseRowFieldUnionMember12ValueUnionMember6,
    DataDatabaseRowFieldUnionMember12ValueUnionMember7,
    DataDatabaseRowFieldUnionMember12ValueUnionMember8,
    DataDatabaseRowFieldUnionMember12ValueUnionMember9,
    DataDatabaseRowFieldUnionMember12ValueUnionMember10,
]


class DataDatabaseRowFieldUnionMember12(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["lookup"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataDatabaseRowFieldUnionMember12InvalidData] = FieldInfo(alias="invalidData", default=None)

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(alias="systemType", default=None)

    value: Optional[DataDatabaseRowFieldUnionMember12Value] = None


DataDatabaseRowField = Union[
    DataDatabaseRowFieldUnionMember0,
    DataDatabaseRowFieldUnionMember1,
    DataDatabaseRowFieldUnionMember2,
    DataDatabaseRowFieldUnionMember3,
    DataDatabaseRowFieldUnionMember4,
    DataDatabaseRowFieldUnionMember5,
    DataDatabaseRowFieldUnionMember6,
    DataDatabaseRowFieldUnionMember7,
    DataDatabaseRowFieldUnionMember8,
    DataDatabaseRowFieldUnionMember9,
    DataDatabaseRowFieldUnionMember10,
    DataDatabaseRowFieldUnionMember11,
    DataDatabaseRowFieldUnionMember12,
]


class DataDatabaseRow(BaseModel):
    id: str

    date_created: str = FieldInfo(alias="dateCreated")

    hid: str

    type: Literal["row"]

    created_by_user_drn: Optional[str] = FieldInfo(alias="createdByUserDrn", default=None)

    creation_block_id: Optional[str] = FieldInfo(alias="creationBlockId", default=None)

    creation_parent_id: Optional[str] = FieldInfo(alias="creationParentId", default=None)

    date_updated: Optional[str] = FieldInfo(alias="dateUpdated", default=None)

    edited_by_user_drn: Optional[str] = FieldInfo(alias="editedByUserDrn", default=None)

    fields: Optional[List[DataDatabaseRowField]] = None

    is_template: Optional[bool] = FieldInfo(alias="isTemplate", default=None)

    name: Optional[str] = None

    parent_id: Optional[str] = FieldInfo(alias="parentId", default=None)

    submission_status: Optional[Literal["draft", "final"]] = FieldInfo(alias="submissionStatus", default=None)

    validation_status: Optional[str] = FieldInfo(alias="validationStatus", default=None)


Data = Annotated[Union[Database, DataDatabaseRow, Workspace], PropertyInfo(discriminator="type")]


class DescribeRowResponse(BaseModel):
    data: Data
