# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Iterable

import httpx

from ..types import ensure_row_ensure_params
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
from ..types.ensure_row_ensure_response import EnsureRowEnsureResponse

__all__ = ["EnsureRowsResource", "AsyncEnsureRowsResource"]


class EnsureRowsResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> EnsureRowsResourceWithRawResponse:
        return EnsureRowsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> EnsureRowsResourceWithStreamingResponse:
        return EnsureRowsResourceWithStreamingResponse(self)

    def ensure(
        self,
        *,
        database_id: str,
        rows: Iterable[ensure_row_ensure_params.Row],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> EnsureRowEnsureResponse:
        """Either creates or updates an existing row.

        Supports updates to both system and
        user defined columns.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/EnsureRows",
            body=maybe_transform(
                {
                    "database_id": database_id,
                    "rows": rows,
                },
                ensure_row_ensure_params.EnsureRowEnsureParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=EnsureRowEnsureResponse,
        )


class AsyncEnsureRowsResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncEnsureRowsResourceWithRawResponse:
        return AsyncEnsureRowsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncEnsureRowsResourceWithStreamingResponse:
        return AsyncEnsureRowsResourceWithStreamingResponse(self)

    async def ensure(
        self,
        *,
        database_id: str,
        rows: Iterable[ensure_row_ensure_params.Row],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> EnsureRowEnsureResponse:
        """Either creates or updates an existing row.

        Supports updates to both system and
        user defined columns.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/EnsureRows",
            body=await async_maybe_transform(
                {
                    "database_id": database_id,
                    "rows": rows,
                },
                ensure_row_ensure_params.EnsureRowEnsureParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=EnsureRowEnsureResponse,
        )


class EnsureRowsResourceWithRawResponse:
    def __init__(self, ensure_rows: EnsureRowsResource) -> None:
        self._ensure_rows = ensure_rows

        self.ensure = to_raw_response_wrapper(
            ensure_rows.ensure,
        )


class AsyncEnsureRowsResourceWithRawResponse:
    def __init__(self, ensure_rows: AsyncEnsureRowsResource) -> None:
        self._ensure_rows = ensure_rows

        self.ensure = async_to_raw_response_wrapper(
            ensure_rows.ensure,
        )


class EnsureRowsResourceWithStreamingResponse:
    def __init__(self, ensure_rows: EnsureRowsResource) -> None:
        self._ensure_rows = ensure_rows

        self.ensure = to_streamed_response_wrapper(
            ensure_rows.ensure,
        )


class AsyncEnsureRowsResourceWithStreamingResponse:
    def __init__(self, ensure_rows: AsyncEnsureRowsResource) -> None:
        self._ensure_rows = ensure_rows

        self.ensure = async_to_streamed_response_wrapper(
            ensure_rows.ensure,
        )
