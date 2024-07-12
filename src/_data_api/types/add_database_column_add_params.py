# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing_extensions import Required, Annotated, TypedDict

from ..types import add_database_column_add_params
from .._utils import PropertyInfo

__all__ = ["AddDatabaseColumnAddParams"]


class AddDatabaseColumnAddParams(TypedDict, total=False):
    column: Required[add_database_column_add_params.Column]

    database_id: Required[Annotated[str, PropertyInfo(alias="databaseId")]]
