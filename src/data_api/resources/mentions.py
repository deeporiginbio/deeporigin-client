# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..types import mention_list_params
from .._types import NOT_GIVEN, Body, Query, Headers, NotGiven
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
from ..types.mention_list_response import MentionListResponse

__all__ = ["MentionsResource", "AsyncMentionsResource"]


class MentionsResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> MentionsResourceWithRawResponse:
        return MentionsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> MentionsResourceWithStreamingResponse:
        return MentionsResourceWithStreamingResponse(self)

    def list(
        self,
        *,
        query: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> MentionListResponse:
        """
        Returns entities that can be mentioned in a body document.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/ListMentions",
            body=maybe_transform({"query": query}, mention_list_params.MentionListParams),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=MentionListResponse,
        )


class AsyncMentionsResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncMentionsResourceWithRawResponse:
        return AsyncMentionsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncMentionsResourceWithStreamingResponse:
        return AsyncMentionsResourceWithStreamingResponse(self)

    async def list(
        self,
        *,
        query: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> MentionListResponse:
        """
        Returns entities that can be mentioned in a body document.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/ListMentions",
            body=await async_maybe_transform({"query": query}, mention_list_params.MentionListParams),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=MentionListResponse,
        )


class MentionsResourceWithRawResponse:
    def __init__(self, mentions: MentionsResource) -> None:
        self._mentions = mentions

        self.list = to_raw_response_wrapper(
            mentions.list,
        )


class AsyncMentionsResourceWithRawResponse:
    def __init__(self, mentions: AsyncMentionsResource) -> None:
        self._mentions = mentions

        self.list = async_to_raw_response_wrapper(
            mentions.list,
        )


class MentionsResourceWithStreamingResponse:
    def __init__(self, mentions: MentionsResource) -> None:
        self._mentions = mentions

        self.list = to_streamed_response_wrapper(
            mentions.list,
        )


class AsyncMentionsResourceWithStreamingResponse:
    def __init__(self, mentions: AsyncMentionsResource) -> None:
        self._mentions = mentions

        self.list = async_to_streamed_response_wrapper(
            mentions.list,
        )
