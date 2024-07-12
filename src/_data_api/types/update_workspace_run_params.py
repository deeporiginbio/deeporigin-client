# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Optional
from typing_extensions import Required, Annotated, TypedDict

from .._utils import PropertyInfo

__all__ = ["UpdateWorkspaceRunParams", "Workspace"]


class UpdateWorkspaceRunParams(TypedDict, total=False):
    id: Required[str]

    workspace: Required[Workspace]


class Workspace(TypedDict, total=False):
    editor: object

    hid: str

    name: Optional[str]

    parent_id: Annotated[Optional[str], PropertyInfo(alias="parentId")]
