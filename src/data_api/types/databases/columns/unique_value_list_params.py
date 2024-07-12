# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Union
from typing_extensions import Literal, Required, Annotated, TypedDict

from ...._utils import PropertyInfo

__all__ = ["UniqueValueListParams", "Variant0", "Variant1"]


class Variant0(TypedDict, total=False):
    column_id: Required[Annotated[str, PropertyInfo(alias="columnId")]]


class Variant1(TypedDict, total=False):
    database_row_id: Required[Annotated[str, PropertyInfo(alias="databaseRowId")]]

    system_column_name: Required[
        Annotated[Literal["creationParentId"], PropertyInfo(alias="systemColumnName")]
    ]


UniqueValueListParams = Union[Variant0, Variant1]
