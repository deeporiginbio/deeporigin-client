# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import List
from typing_extensions import Required, Annotated, TypedDict

from .._utils import PropertyInfo

__all__ = ["DescribeRowRetrieveParams", "ColumnSelection"]


class DescribeRowRetrieveParams(TypedDict, total=False):
    row_id: Required[Annotated[str, PropertyInfo(alias="rowId")]]

    column_selection: Annotated[ColumnSelection, PropertyInfo(alias="columnSelection")]

    fields: bool


class ColumnSelection(TypedDict, total=False):
    exclude: List[str]

    include: List[str]
