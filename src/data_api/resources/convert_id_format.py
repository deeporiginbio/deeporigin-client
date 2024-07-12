# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Iterable

import httpx

from ..types import convert_id_format_convert_params
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
from ..types.convert_id_format_convert_response import ConvertIDFormatConvertResponse

__all__ = ["ConvertIDFormatResource", "AsyncConvertIDFormatResource"]


class ConvertIDFormatResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> ConvertIDFormatResourceWithRawResponse:
        return ConvertIDFormatResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> ConvertIDFormatResourceWithStreamingResponse:
        return ConvertIDFormatResourceWithStreamingResponse(self)

    def convert(
        self,
        *,
        conversions: Iterable[convert_id_format_convert_params.Conversion],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> ConvertIDFormatConvertResponse:
        """
        Converts between system IDs and human IDs (HIDs).

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/ConvertIdFormat",
            body=maybe_transform(
                {"conversions": conversions}, convert_id_format_convert_params.ConvertIDFormatConvertParams
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=ConvertIDFormatConvertResponse,
        )


class AsyncConvertIDFormatResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncConvertIDFormatResourceWithRawResponse:
        return AsyncConvertIDFormatResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncConvertIDFormatResourceWithStreamingResponse:
        return AsyncConvertIDFormatResourceWithStreamingResponse(self)

    async def convert(
        self,
        *,
        conversions: Iterable[convert_id_format_convert_params.Conversion],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> ConvertIDFormatConvertResponse:
        """
        Converts between system IDs and human IDs (HIDs).

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/ConvertIdFormat",
            body=await async_maybe_transform(
                {"conversions": conversions}, convert_id_format_convert_params.ConvertIDFormatConvertParams
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=ConvertIDFormatConvertResponse,
        )


class ConvertIDFormatResourceWithRawResponse:
    def __init__(self, convert_id_format: ConvertIDFormatResource) -> None:
        self._convert_id_format = convert_id_format

        self.convert = to_raw_response_wrapper(
            convert_id_format.convert,
        )


class AsyncConvertIDFormatResourceWithRawResponse:
    def __init__(self, convert_id_format: AsyncConvertIDFormatResource) -> None:
        self._convert_id_format = convert_id_format

        self.convert = async_to_raw_response_wrapper(
            convert_id_format.convert,
        )


class ConvertIDFormatResourceWithStreamingResponse:
    def __init__(self, convert_id_format: ConvertIDFormatResource) -> None:
        self._convert_id_format = convert_id_format

        self.convert = to_streamed_response_wrapper(
            convert_id_format.convert,
        )


class AsyncConvertIDFormatResourceWithStreamingResponse:
    def __init__(self, convert_id_format: AsyncConvertIDFormatResource) -> None:
        self._convert_id_format = convert_id_format

        self.convert = async_to_streamed_response_wrapper(
            convert_id_format.convert,
        )
