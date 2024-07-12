# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing_extensions import Required, Annotated, TypedDict

from .._utils import PropertyInfo

__all__ = ["SequenceParseParams"]


class SequenceParseParams(TypedDict, total=False):
    file_id: Required[Annotated[str, PropertyInfo(alias="fileId")]]
