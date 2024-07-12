# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Optional
from typing_extensions import Required, Annotated, TypedDict

from .._utils import PropertyInfo

__all__ = ["CreateFileUploadCreateParams"]


class CreateFileUploadCreateParams(TypedDict, total=False):
    content_length: Required[Annotated[str, PropertyInfo(alias="contentLength")]]

    content_type: Required[Annotated[Optional[str], PropertyInfo(alias="contentType")]]

    name: Required[str]
