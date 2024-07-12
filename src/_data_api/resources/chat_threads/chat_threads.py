# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from .messages import (
    MessagesResource,
    AsyncMessagesResource,
    MessagesResourceWithRawResponse,
    AsyncMessagesResourceWithRawResponse,
    MessagesResourceWithStreamingResponse,
    AsyncMessagesResourceWithStreamingResponse,
)
from ..._compat import cached_property
from ..._resource import SyncAPIResource, AsyncAPIResource

__all__ = ["ChatThreadsResource", "AsyncChatThreadsResource"]


class ChatThreadsResource(SyncAPIResource):
    @cached_property
    def messages(self) -> MessagesResource:
        return MessagesResource(self._client)

    @cached_property
    def with_raw_response(self) -> ChatThreadsResourceWithRawResponse:
        return ChatThreadsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> ChatThreadsResourceWithStreamingResponse:
        return ChatThreadsResourceWithStreamingResponse(self)


class AsyncChatThreadsResource(AsyncAPIResource):
    @cached_property
    def messages(self) -> AsyncMessagesResource:
        return AsyncMessagesResource(self._client)

    @cached_property
    def with_raw_response(self) -> AsyncChatThreadsResourceWithRawResponse:
        return AsyncChatThreadsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncChatThreadsResourceWithStreamingResponse:
        return AsyncChatThreadsResourceWithStreamingResponse(self)


class ChatThreadsResourceWithRawResponse:
    def __init__(self, chat_threads: ChatThreadsResource) -> None:
        self._chat_threads = chat_threads

    @cached_property
    def messages(self) -> MessagesResourceWithRawResponse:
        return MessagesResourceWithRawResponse(self._chat_threads.messages)


class AsyncChatThreadsResourceWithRawResponse:
    def __init__(self, chat_threads: AsyncChatThreadsResource) -> None:
        self._chat_threads = chat_threads

    @cached_property
    def messages(self) -> AsyncMessagesResourceWithRawResponse:
        return AsyncMessagesResourceWithRawResponse(self._chat_threads.messages)


class ChatThreadsResourceWithStreamingResponse:
    def __init__(self, chat_threads: ChatThreadsResource) -> None:
        self._chat_threads = chat_threads

    @cached_property
    def messages(self) -> MessagesResourceWithStreamingResponse:
        return MessagesResourceWithStreamingResponse(self._chat_threads.messages)


class AsyncChatThreadsResourceWithStreamingResponse:
    def __init__(self, chat_threads: AsyncChatThreadsResource) -> None:
        self._chat_threads = chat_threads

    @cached_property
    def messages(self) -> AsyncMessagesResourceWithStreamingResponse:
        return AsyncMessagesResourceWithStreamingResponse(self._chat_threads.messages)
