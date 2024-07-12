# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Union, Iterable
from typing_extensions import Required, Annotated, TypedDict

from ..._utils import PropertyInfo

__all__ = [
    "HierarchyListParams",
    "Filter",
    "FilterParent",
    "FilterParentParent",
    "FilterParentParentID",
    "FilterParentParentIsRoot",
    "FilterRowType",
]


class HierarchyListParams(TypedDict, total=False):
    filters: Required[Iterable[Filter]]


class FilterParentParentID(TypedDict, total=False):
    id: Required[str]


class FilterParentParentIsRoot(TypedDict, total=False):
    is_root: Required[Annotated[bool, PropertyInfo(alias="isRoot")]]


FilterParentParent = Union[FilterParentParentID, FilterParentParentIsRoot]


class FilterParent(TypedDict, total=False):
    parent: Required[FilterParentParent]


class FilterRowType(TypedDict, total=False):
    row_type: Required[Annotated[str, PropertyInfo(alias="rowType")]]


Filter = Union[FilterParent, FilterRowType]
