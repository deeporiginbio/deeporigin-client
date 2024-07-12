# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import List, Iterable
from typing_extensions import Annotated, TypedDict

from .._utils import PropertyInfo

__all__ = ["FileListParams", "Filter"]


class FileListParams(TypedDict, total=False):
    filters: Iterable[Filter]


class Filter(TypedDict, total=False):
    assigned_row_ids: Annotated[List[str], PropertyInfo(alias="assignedRowIds")]

    is_unassigned: Annotated[bool, PropertyInfo(alias="isUnassigned")]
