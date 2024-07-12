# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Iterable, Optional
from typing_extensions import Literal, Required, Annotated, TypedDict

from .._utils import PropertyInfo

__all__ = ["EnsureRowEnsureParams", "Row", "RowCell", "RowRow"]


class EnsureRowEnsureParams(TypedDict, total=False):
    database_id: Required[Annotated[str, PropertyInfo(alias="databaseId")]]

    rows: Required[Iterable[Row]]


class RowCell(TypedDict, total=False):
    column_id: Required[Annotated[str, PropertyInfo(alias="columnId")]]
    """The programmatic key or system ID for the column."""

    cell_id: Annotated[str, PropertyInfo(alias="cellId")]

    value: object


class RowRow(TypedDict, total=False):
    creation_block_id: Annotated[Optional[str], PropertyInfo(alias="creationBlockId")]

    creation_parent_id: Annotated[Optional[str], PropertyInfo(alias="creationParentId")]

    is_template: Annotated[Optional[bool], PropertyInfo(alias="isTemplate")]

    submission_status: Annotated[
        Literal["draft", "final"], PropertyInfo(alias="submissionStatus")
    ]


class Row(TypedDict, total=False):
    cells: Iterable[RowCell]

    row: RowRow

    row_id: Annotated[str, PropertyInfo(alias="rowId")]
