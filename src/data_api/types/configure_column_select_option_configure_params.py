# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Iterable
from typing_extensions import Literal, Required, Annotated, TypedDict

from .._utils import PropertyInfo

__all__ = ["ConfigureColumnSelectOptionConfigureParams", "OptionConfiguration"]


class ConfigureColumnSelectOptionConfigureParams(TypedDict, total=False):
    column_id: Required[Annotated[str, PropertyInfo(alias="columnId")]]

    option_configuration: Required[Annotated[Iterable[OptionConfiguration], PropertyInfo(alias="optionConfiguration")]]


class OptionConfiguration(TypedDict, total=False):
    op: Required[Literal["add", "remove"]]

    option: Required[str]
