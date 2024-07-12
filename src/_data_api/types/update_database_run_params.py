# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Optional
from typing_extensions import Required, Annotated, TypedDict

from .._utils import PropertyInfo

__all__ = ["UpdateDatabaseRunParams", "Database"]


class UpdateDatabaseRunParams(TypedDict, total=False):
    id: Required[str]

    database: Required[Database]


class Database(TypedDict, total=False):
    editor: object

    hid: str

    hid_prefix: Annotated[Optional[str], PropertyInfo(alias="hidPrefix")]

    name: Optional[str]

    parent_id: Annotated[Optional[str], PropertyInfo(alias="parentId")]
