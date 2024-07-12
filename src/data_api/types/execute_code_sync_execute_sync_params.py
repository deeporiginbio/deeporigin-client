# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing_extensions import Literal, Required, Annotated, TypedDict

from .._utils import PropertyInfo

__all__ = ["ExecuteCodeSyncExecuteSyncParams"]


class ExecuteCodeSyncExecuteSyncParams(TypedDict, total=False):
    code: Required[str]

    code_language: Required[Annotated[Literal["python"], PropertyInfo(alias="codeLanguage")]]
