# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import List, Union, Iterable
from typing_extensions import Literal, Required, Annotated, TypedDict

from ..._utils import PropertyInfo

__all__ = [
    "RowListParams",
    "ColumnSelection",
    "Filter",
    "FilterUnionMember0",
    "FilterUnionMember1",
    "FilterRowFilterJoin",
    "FilterRowFilterJoinCondition",
    "FilterRowFilterJoinConditionRowFilterJoin",
    "FilterRowFilterJoinConditionRowFilterJoinCondition",
    "FilterRowFilterJoinConditionRowFilterJoinConditionUnionMember0",
    "FilterRowFilterJoinConditionRowFilterJoinConditionUnionMember1",
    "FilterRowFilterJoinConditionUnionMember1",
    "FilterRowFilterJoinConditionUnionMember2",
]


class RowListParams(TypedDict, total=False):
    database_row_id: Required[Annotated[str, PropertyInfo(alias="databaseRowId")]]

    column_selection: Annotated[ColumnSelection, PropertyInfo(alias="columnSelection")]

    creation_block_id: Annotated[str, PropertyInfo(alias="creationBlockId")]

    creation_parent_id: Annotated[str, PropertyInfo(alias="creationParentId")]

    filter: Filter


class ColumnSelection(TypedDict, total=False):
    exclude: List[str]

    include: List[str]


class FilterUnionMember0(TypedDict, total=False):
    column_id: Required[Annotated[str, PropertyInfo(alias="columnId")]]

    filter_type: Required[Annotated[Literal["text"], PropertyInfo(alias="filterType")]]

    filter_value: Required[Annotated[str, PropertyInfo(alias="filterValue")]]

    operator: Required[Literal["equals", "notEqual", "blank", "notBlank"]]


class FilterUnionMember1(TypedDict, total=False):
    column_id: Required[Annotated[str, PropertyInfo(alias="columnId")]]

    filter_type: Required[Annotated[Literal["number"], PropertyInfo(alias="filterType")]]

    filter_value: Required[Annotated[float, PropertyInfo(alias="filterValue")]]

    operator: Required[
        Literal[
            "equals",
            "notEqual",
            "lessThan",
            "lessThanOrEqual",
            "greaterThan",
            "greaterThanOrEqual",
            "blank",
            "notBlank",
        ]
    ]


class FilterRowFilterJoinConditionRowFilterJoinConditionUnionMember0(TypedDict, total=False):
    column_id: Required[Annotated[str, PropertyInfo(alias="columnId")]]

    filter_type: Required[Annotated[Literal["text"], PropertyInfo(alias="filterType")]]

    filter_value: Required[Annotated[str, PropertyInfo(alias="filterValue")]]

    operator: Required[Literal["equals", "notEqual", "blank", "notBlank"]]


class FilterRowFilterJoinConditionRowFilterJoinConditionUnionMember1(TypedDict, total=False):
    column_id: Required[Annotated[str, PropertyInfo(alias="columnId")]]

    filter_type: Required[Annotated[Literal["number"], PropertyInfo(alias="filterType")]]

    filter_value: Required[Annotated[float, PropertyInfo(alias="filterValue")]]

    operator: Required[
        Literal[
            "equals",
            "notEqual",
            "lessThan",
            "lessThanOrEqual",
            "greaterThan",
            "greaterThanOrEqual",
            "blank",
            "notBlank",
        ]
    ]


FilterRowFilterJoinConditionRowFilterJoinCondition = Union[
    FilterRowFilterJoinConditionRowFilterJoinConditionUnionMember0,
    FilterRowFilterJoinConditionRowFilterJoinConditionUnionMember1,
    object,
]


class FilterRowFilterJoinConditionRowFilterJoin(TypedDict, total=False):
    conditions: Required[Iterable[FilterRowFilterJoinConditionRowFilterJoinCondition]]

    filter_type: Required[Annotated[Literal["join"], PropertyInfo(alias="filterType")]]

    join_type: Required[Annotated[Literal["and", "or"], PropertyInfo(alias="joinType")]]


class FilterRowFilterJoinConditionUnionMember1(TypedDict, total=False):
    column_id: Required[Annotated[str, PropertyInfo(alias="columnId")]]

    filter_type: Required[Annotated[Literal["text"], PropertyInfo(alias="filterType")]]

    filter_value: Required[Annotated[str, PropertyInfo(alias="filterValue")]]

    operator: Required[Literal["equals", "notEqual", "blank", "notBlank"]]


class FilterRowFilterJoinConditionUnionMember2(TypedDict, total=False):
    column_id: Required[Annotated[str, PropertyInfo(alias="columnId")]]

    filter_type: Required[Annotated[Literal["number"], PropertyInfo(alias="filterType")]]

    filter_value: Required[Annotated[float, PropertyInfo(alias="filterValue")]]

    operator: Required[
        Literal[
            "equals",
            "notEqual",
            "lessThan",
            "lessThanOrEqual",
            "greaterThan",
            "greaterThanOrEqual",
            "blank",
            "notBlank",
        ]
    ]


FilterRowFilterJoinCondition = Union[
    FilterRowFilterJoinConditionRowFilterJoin,
    FilterRowFilterJoinConditionUnionMember1,
    FilterRowFilterJoinConditionUnionMember2,
]


class FilterRowFilterJoin(TypedDict, total=False):
    conditions: Required[Iterable[FilterRowFilterJoinCondition]]

    filter_type: Required[Annotated[Literal["join"], PropertyInfo(alias="filterType")]]

    join_type: Required[Annotated[Literal["and", "or"], PropertyInfo(alias="joinType")]]


Filter = Union[FilterUnionMember0, FilterUnionMember1, FilterRowFilterJoin]
