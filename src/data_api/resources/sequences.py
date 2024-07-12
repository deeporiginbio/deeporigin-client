# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..types import sequence_parse_params
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
from ..types.sequence_parse_response import SequenceParseResponse

__all__ = ["SequencesResource", "AsyncSequencesResource"]


class SequencesResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> SequencesResourceWithRawResponse:
        return SequencesResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> SequencesResourceWithStreamingResponse:
        return SequencesResourceWithStreamingResponse(self)

    def parse(
        self,
        *,
        file_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> SequenceParseResponse:
        """
        Parses a base sequence file and returns the parsed result.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/ParseBaseSequenceData",
            body=maybe_transform(
                {"file_id": file_id}, sequence_parse_params.SequenceParseParams
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=SequenceParseResponse,
        )


class AsyncSequencesResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncSequencesResourceWithRawResponse:
        return AsyncSequencesResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncSequencesResourceWithStreamingResponse:
        return AsyncSequencesResourceWithStreamingResponse(self)

    async def parse(
        self,
        *,
        file_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> SequenceParseResponse:
        """
        Parses a base sequence file and returns the parsed result.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/ParseBaseSequenceData",
            body=await async_maybe_transform(
                {"file_id": file_id}, sequence_parse_params.SequenceParseParams
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=SequenceParseResponse,
        )


class SequencesResourceWithRawResponse:
    def __init__(self, sequences: SequencesResource) -> None:
        self._sequences = sequences

        self.parse = to_raw_response_wrapper(
            sequences.parse,
        )


class AsyncSequencesResourceWithRawResponse:
    def __init__(self, sequences: AsyncSequencesResource) -> None:
        self._sequences = sequences

        self.parse = async_to_raw_response_wrapper(
            sequences.parse,
        )


class SequencesResourceWithStreamingResponse:
    def __init__(self, sequences: SequencesResource) -> None:
        self._sequences = sequences

        self.parse = to_streamed_response_wrapper(
            sequences.parse,
        )


class AsyncSequencesResourceWithStreamingResponse:
    def __init__(self, sequences: AsyncSequencesResource) -> None:
        self._sequences = sequences

        self.parse = async_to_streamed_response_wrapper(
            sequences.parse,
        )
