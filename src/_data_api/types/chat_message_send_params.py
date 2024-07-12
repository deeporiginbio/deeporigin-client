# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Iterable
from typing_extensions import Literal, Required, Annotated, TypedDict

from .._utils import PropertyInfo

__all__ = [
    "ChatMessageSendParams",
    "Message",
    "Context",
    "ContextDatabase",
    "ContextDatabaseDatabase",
    "ContextDatabaseColumn",
    "ContextDatabaseRow",
]


class ChatMessageSendParams(TypedDict, total=False):
    messages: Required[Iterable[Message]]

    thread_id: Required[Annotated[str, PropertyInfo(alias="threadId")]]

    context: Context


class Message(TypedDict, total=False):
    content: Required[str]

    role: Required[Literal["user", "assistant"]]


class ContextDatabaseDatabase(TypedDict, total=False):
    hid: Required[str]

    hid_prefix: Required[Annotated[str, PropertyInfo(alias="hidPrefix")]]

    name: Required[str]


class ContextDatabaseColumn(TypedDict, total=False):
    name: Required[str]


class ContextDatabaseRow(TypedDict, total=False):
    hid: Required[str]

    name: str


class ContextDatabase(TypedDict, total=False):
    database: Required[ContextDatabaseDatabase]

    columns: Iterable[ContextDatabaseColumn]

    rows: Iterable[ContextDatabaseRow]


class Context(TypedDict, total=False):
    databases: Iterable[ContextDatabase]
