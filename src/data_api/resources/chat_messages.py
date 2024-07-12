# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Iterable

import httpx

from ..types import chat_message_send_params
from .._types import NOT_GIVEN, Body, Query, Headers, NoneType, NotGiven
from .._utils import (
    maybe_transform,
    async_maybe_transform,
)
from .._compat import cached_property
from .._resource import SyncAPIResource, AsyncAPIResource
from .._response import (
    to_raw_response_wrapper,
    to_streamed_response_wrapper,
    async_to_raw_response_wrapper,
    async_to_streamed_response_wrapper,
)
from .._base_client import (
    make_request_options,
)

__all__ = ["ChatMessagesResource", "AsyncChatMessagesResource"]


class ChatMessagesResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> ChatMessagesResourceWithRawResponse:
        return ChatMessagesResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> ChatMessagesResourceWithStreamingResponse:
        return ChatMessagesResourceWithStreamingResponse(self)

    def send(
        self,
        *,
        messages: Iterable[chat_message_send_params.Message],
        thread_id: str,
        context: chat_message_send_params.Context | NotGiven = NOT_GIVEN,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> None:
        """
        Send a chat message to the Deep Origin assistant and streams results via SSE.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        extra_headers = {"Accept": "*/*", **(extra_headers or {})}
        return self._post(
            "/SendChatMessage",
            body=maybe_transform(
                {
                    "messages": messages,
                    "thread_id": thread_id,
                    "context": context,
                },
                chat_message_send_params.ChatMessageSendParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=NoneType,
        )


class AsyncChatMessagesResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncChatMessagesResourceWithRawResponse:
        return AsyncChatMessagesResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncChatMessagesResourceWithStreamingResponse:
        return AsyncChatMessagesResourceWithStreamingResponse(self)

    async def send(
        self,
        *,
        messages: Iterable[chat_message_send_params.Message],
        thread_id: str,
        context: chat_message_send_params.Context | NotGiven = NOT_GIVEN,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> None:
        """
        Send a chat message to the Deep Origin assistant and streams results via SSE.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        extra_headers = {"Accept": "*/*", **(extra_headers or {})}
        return await self._post(
            "/SendChatMessage",
            body=await async_maybe_transform(
                {
                    "messages": messages,
                    "thread_id": thread_id,
                    "context": context,
                },
                chat_message_send_params.ChatMessageSendParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=NoneType,
        )


class ChatMessagesResourceWithRawResponse:
    def __init__(self, chat_messages: ChatMessagesResource) -> None:
        self._chat_messages = chat_messages

        self.send = to_raw_response_wrapper(
            chat_messages.send,
        )


class AsyncChatMessagesResourceWithRawResponse:
    def __init__(self, chat_messages: AsyncChatMessagesResource) -> None:
        self._chat_messages = chat_messages

        self.send = async_to_raw_response_wrapper(
            chat_messages.send,
        )


class ChatMessagesResourceWithStreamingResponse:
    def __init__(self, chat_messages: ChatMessagesResource) -> None:
        self._chat_messages = chat_messages

        self.send = to_streamed_response_wrapper(
            chat_messages.send,
        )


class AsyncChatMessagesResourceWithStreamingResponse:
    def __init__(self, chat_messages: AsyncChatMessagesResource) -> None:
        self._chat_messages = chat_messages

        self.send = async_to_streamed_response_wrapper(
            chat_messages.send,
        )
