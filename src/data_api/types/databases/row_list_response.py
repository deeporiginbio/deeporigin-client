# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from typing import List, Union, Optional
from typing_extensions import Literal

from pydantic import Field as FieldInfo

from ..._models import BaseModel

__all__ = [
    "RowListResponse",
    "Data",
    "DataField",
    "DataFieldUnionMember0",
    "DataFieldUnionMember0InvalidData",
    "DataFieldUnionMember1",
    "DataFieldUnionMember1InvalidData",
    "DataFieldUnionMember2",
    "DataFieldUnionMember2InvalidData",
    "DataFieldUnionMember3",
    "DataFieldUnionMember3InvalidData",
    "DataFieldUnionMember4",
    "DataFieldUnionMember4InvalidData",
    "DataFieldUnionMember4Value",
    "DataFieldUnionMember5",
    "DataFieldUnionMember5InvalidData",
    "DataFieldUnionMember5Value",
    "DataFieldUnionMember6",
    "DataFieldUnionMember6InvalidData",
    "DataFieldUnionMember6Value",
    "DataFieldUnionMember7",
    "DataFieldUnionMember7InvalidData",
    "DataFieldUnionMember7Value",
    "DataFieldUnionMember8",
    "DataFieldUnionMember8InvalidData",
    "DataFieldUnionMember9",
    "DataFieldUnionMember9InvalidData",
    "DataFieldUnionMember9Value",
    "DataFieldUnionMember9ValueURL",
    "DataFieldUnionMember10",
    "DataFieldUnionMember10InvalidData",
    "DataFieldUnionMember10Value",
    "DataFieldUnionMember11",
    "DataFieldUnionMember11Value",
    "DataFieldUnionMember11InvalidData",
    "DataFieldUnionMember12",
    "DataFieldUnionMember12InvalidData",
    "DataFieldUnionMember12Value",
    "DataFieldUnionMember12ValueUnionMember0",
    "DataFieldUnionMember12ValueUnionMember0InvalidData",
    "DataFieldUnionMember12ValueUnionMember1",
    "DataFieldUnionMember12ValueUnionMember1InvalidData",
    "DataFieldUnionMember12ValueUnionMember2",
    "DataFieldUnionMember12ValueUnionMember2InvalidData",
    "DataFieldUnionMember12ValueUnionMember3",
    "DataFieldUnionMember12ValueUnionMember3InvalidData",
    "DataFieldUnionMember12ValueUnionMember4",
    "DataFieldUnionMember12ValueUnionMember4InvalidData",
    "DataFieldUnionMember12ValueUnionMember4Value",
    "DataFieldUnionMember12ValueUnionMember5",
    "DataFieldUnionMember12ValueUnionMember5InvalidData",
    "DataFieldUnionMember12ValueUnionMember5Value",
    "DataFieldUnionMember12ValueUnionMember6",
    "DataFieldUnionMember12ValueUnionMember6InvalidData",
    "DataFieldUnionMember12ValueUnionMember6Value",
    "DataFieldUnionMember12ValueUnionMember7",
    "DataFieldUnionMember12ValueUnionMember7InvalidData",
    "DataFieldUnionMember12ValueUnionMember7Value",
    "DataFieldUnionMember12ValueUnionMember8",
    "DataFieldUnionMember12ValueUnionMember8InvalidData",
    "DataFieldUnionMember12ValueUnionMember9",
    "DataFieldUnionMember12ValueUnionMember9InvalidData",
    "DataFieldUnionMember12ValueUnionMember9Value",
    "DataFieldUnionMember12ValueUnionMember9ValueURL",
    "DataFieldUnionMember12ValueUnionMember10",
    "DataFieldUnionMember12ValueUnionMember10InvalidData",
    "DataFieldUnionMember12ValueUnionMember10Value",
]


class DataFieldUnionMember0InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember0(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["text"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember0InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[str] = None


class DataFieldUnionMember1InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember1(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["integer"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember1InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[float] = None


class DataFieldUnionMember2InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember2(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["float"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember2InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[float] = None


class DataFieldUnionMember3InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember3(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["boolean"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember3InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[bool] = None


class DataFieldUnionMember4InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember4Value(BaseModel):
    row_ids: List[str] = FieldInfo(alias="rowIds")


class DataFieldUnionMember4(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["reference"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember4InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember4Value] = None


class DataFieldUnionMember5InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember5Value(BaseModel):
    top_level_blocks: List[object] = FieldInfo(alias="topLevelBlocks")


class DataFieldUnionMember5(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["editor"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember5InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember5Value] = None


class DataFieldUnionMember6InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember6Value(BaseModel):
    file_ids: List[str] = FieldInfo(alias="fileIds")


class DataFieldUnionMember6(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["file"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember6InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember6Value] = None


class DataFieldUnionMember7InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember7Value(BaseModel):
    selected_options: List[str] = FieldInfo(alias="selectedOptions")


class DataFieldUnionMember7(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["select"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember7InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember7Value] = None


class DataFieldUnionMember8InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember8(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["date"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember8InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[str] = None


class DataFieldUnionMember9InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember9ValueURL(BaseModel):
    url: str

    title: Optional[str] = None


class DataFieldUnionMember9Value(BaseModel):
    urls: List[DataFieldUnionMember9ValueURL]


class DataFieldUnionMember9(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["url"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember9InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember9Value] = None


class DataFieldUnionMember10InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember10Value(BaseModel):
    user_drns: List[str] = FieldInfo(alias="userDrns")


class DataFieldUnionMember10(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["user"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember10InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember10Value] = None


class DataFieldUnionMember11Value(BaseModel):
    error_message: Optional[str] = FieldInfo(alias="errorMessage", default=None)
    """The error message from executing the expression, if one occurred."""

    invalid_result: Optional[object] = FieldInfo(alias="invalidResult", default=None)
    """Expression result that is not a valid return value type."""

    result: Union[str, float, None] = None
    """The return value from executing the expression."""


class DataFieldUnionMember11InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember11(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["expression"]

    validation_status: str = FieldInfo(alias="validationStatus")

    value: DataFieldUnionMember11Value

    invalid_data: Optional[DataFieldUnionMember11InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )


class DataFieldUnionMember12InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember0InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember0(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["text"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember12ValueUnionMember0InvalidData] = (
        FieldInfo(alias="invalidData", default=None)
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember1InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember1(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["integer"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember12ValueUnionMember1InvalidData] = (
        FieldInfo(alias="invalidData", default=None)
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[float] = None


class DataFieldUnionMember12ValueUnionMember2InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember2(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["float"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember12ValueUnionMember2InvalidData] = (
        FieldInfo(alias="invalidData", default=None)
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[float] = None


class DataFieldUnionMember12ValueUnionMember3InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember3(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["boolean"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember12ValueUnionMember3InvalidData] = (
        FieldInfo(alias="invalidData", default=None)
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[bool] = None


class DataFieldUnionMember12ValueUnionMember4InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember4Value(BaseModel):
    row_ids: List[str] = FieldInfo(alias="rowIds")


class DataFieldUnionMember12ValueUnionMember4(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["reference"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember12ValueUnionMember4InvalidData] = (
        FieldInfo(alias="invalidData", default=None)
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember12ValueUnionMember4Value] = None


class DataFieldUnionMember12ValueUnionMember5InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember5Value(BaseModel):
    top_level_blocks: List[object] = FieldInfo(alias="topLevelBlocks")


class DataFieldUnionMember12ValueUnionMember5(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["editor"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember12ValueUnionMember5InvalidData] = (
        FieldInfo(alias="invalidData", default=None)
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember12ValueUnionMember5Value] = None


class DataFieldUnionMember12ValueUnionMember6InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember6Value(BaseModel):
    file_ids: List[str] = FieldInfo(alias="fileIds")


class DataFieldUnionMember12ValueUnionMember6(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["file"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember12ValueUnionMember6InvalidData] = (
        FieldInfo(alias="invalidData", default=None)
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember12ValueUnionMember6Value] = None


class DataFieldUnionMember12ValueUnionMember7InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember7Value(BaseModel):
    selected_options: List[str] = FieldInfo(alias="selectedOptions")


class DataFieldUnionMember12ValueUnionMember7(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["select"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember12ValueUnionMember7InvalidData] = (
        FieldInfo(alias="invalidData", default=None)
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember12ValueUnionMember7Value] = None


class DataFieldUnionMember12ValueUnionMember8InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember8(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["date"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember12ValueUnionMember8InvalidData] = (
        FieldInfo(alias="invalidData", default=None)
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember9InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember9ValueURL(BaseModel):
    url: str

    title: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember9Value(BaseModel):
    urls: List[DataFieldUnionMember12ValueUnionMember9ValueURL]


class DataFieldUnionMember12ValueUnionMember9(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["url"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember12ValueUnionMember9InvalidData] = (
        FieldInfo(alias="invalidData", default=None)
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember12ValueUnionMember9Value] = None


class DataFieldUnionMember12ValueUnionMember10InvalidData(BaseModel):
    message: Optional[str] = None


class DataFieldUnionMember12ValueUnionMember10Value(BaseModel):
    user_drns: List[str] = FieldInfo(alias="userDrns")


class DataFieldUnionMember12ValueUnionMember10(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["user"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember12ValueUnionMember10InvalidData] = (
        FieldInfo(alias="invalidData", default=None)
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember12ValueUnionMember10Value] = None


DataFieldUnionMember12Value = Union[
    DataFieldUnionMember12ValueUnionMember0,
    DataFieldUnionMember12ValueUnionMember1,
    DataFieldUnionMember12ValueUnionMember2,
    DataFieldUnionMember12ValueUnionMember3,
    DataFieldUnionMember12ValueUnionMember4,
    DataFieldUnionMember12ValueUnionMember5,
    DataFieldUnionMember12ValueUnionMember6,
    DataFieldUnionMember12ValueUnionMember7,
    DataFieldUnionMember12ValueUnionMember8,
    DataFieldUnionMember12ValueUnionMember9,
    DataFieldUnionMember12ValueUnionMember10,
]


class DataFieldUnionMember12(BaseModel):
    cell_id: str = FieldInfo(alias="cellId")

    column_id: str = FieldInfo(alias="columnId")

    type: Literal["lookup"]

    validation_status: str = FieldInfo(alias="validationStatus")

    invalid_data: Optional[DataFieldUnionMember12InvalidData] = FieldInfo(
        alias="invalidData", default=None
    )

    system_type: Optional[Literal["name", "bodyDocument"]] = FieldInfo(
        alias="systemType", default=None
    )

    value: Optional[DataFieldUnionMember12Value] = None


DataField = Union[
    DataFieldUnionMember0,
    DataFieldUnionMember1,
    DataFieldUnionMember2,
    DataFieldUnionMember3,
    DataFieldUnionMember4,
    DataFieldUnionMember5,
    DataFieldUnionMember6,
    DataFieldUnionMember7,
    DataFieldUnionMember8,
    DataFieldUnionMember9,
    DataFieldUnionMember10,
    DataFieldUnionMember11,
    DataFieldUnionMember12,
]


class Data(BaseModel):
    id: str

    date_created: str = FieldInfo(alias="dateCreated")

    hid: str

    type: Literal["row"]

    created_by_user_drn: Optional[str] = FieldInfo(
        alias="createdByUserDrn", default=None
    )

    creation_block_id: Optional[str] = FieldInfo(alias="creationBlockId", default=None)

    creation_parent_id: Optional[str] = FieldInfo(
        alias="creationParentId", default=None
    )

    date_updated: Optional[str] = FieldInfo(alias="dateUpdated", default=None)

    edited_by_user_drn: Optional[str] = FieldInfo(alias="editedByUserDrn", default=None)

    fields: Optional[List[DataField]] = None

    is_template: Optional[bool] = FieldInfo(alias="isTemplate", default=None)

    name: Optional[str] = None

    parent_id: Optional[str] = FieldInfo(alias="parentId", default=None)

    submission_status: Optional[Literal["draft", "final"]] = FieldInfo(
        alias="submissionStatus", default=None
    )

    validation_status: Optional[str] = FieldInfo(alias="validationStatus", default=None)


class RowListResponse(BaseModel):
    data: List[Data]
