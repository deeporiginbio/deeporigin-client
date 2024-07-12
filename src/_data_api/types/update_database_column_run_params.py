# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing_extensions import Required, Annotated, TypedDict

from ..types import update_database_column_run_params
from .._utils import PropertyInfo

__all__ = ["UpdateDatabaseColumnRunParams"]


class UpdateDatabaseColumnRunParams(TypedDict, total=False):
    column: Required[update_database_column_run_params.Column]

    column_id: Required[Annotated[str, PropertyInfo(alias="columnId")]]
