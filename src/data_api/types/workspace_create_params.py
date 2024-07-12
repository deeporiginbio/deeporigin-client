# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Optional
from typing_extensions import Required, Annotated, TypedDict

from .._utils import PropertyInfo

__all__ = ["WorkspaceCreateParams", "Workspace"]


class WorkspaceCreateParams(TypedDict, total=False):
    workspace: Required[Workspace]


class Workspace(TypedDict, total=False):
    hid: Required[str]

    name: Required[Optional[str]]

    parent_id: Annotated[Optional[str], PropertyInfo(alias="parentId")]
