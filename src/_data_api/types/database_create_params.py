# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Iterable, Optional
from typing_extensions import Required, Annotated, TypedDict

from .._utils import PropertyInfo

__all__ = ["DatabaseCreateParams", "Database"]


class DatabaseCreateParams(TypedDict, total=False):
    database: Required[Database]


class Database(TypedDict, total=False):
    hid: Required[str]

    hid_prefix: Required[Annotated[Optional[str], PropertyInfo(alias="hidPrefix")]]

    name: Required[Optional[str]]

    cols: Iterable[object]

    parent_id: Annotated[Optional[str], PropertyInfo(alias="parentId")]
