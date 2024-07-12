# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import List

import httpx

from ..types import delete_row_delete_params
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
from ..types.delete_row_delete_response import DeleteRowDeleteResponse

__all__ = ["DeleteRowsResource", "AsyncDeleteRowsResource"]


class DeleteRowsResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> DeleteRowsResourceWithRawResponse:
        return DeleteRowsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> DeleteRowsResourceWithStreamingResponse:
        return DeleteRowsResourceWithStreamingResponse(self)

    def delete(
        self,
        *,
        row_ids: List[str],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> DeleteRowDeleteResponse:
        """
        Delete rows by their ids.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/DeleteRows",
            body=maybe_transform({"row_ids": row_ids}, delete_row_delete_params.DeleteRowDeleteParams),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=DeleteRowDeleteResponse,
        )


class AsyncDeleteRowsResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncDeleteRowsResourceWithRawResponse:
        return AsyncDeleteRowsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncDeleteRowsResourceWithStreamingResponse:
        return AsyncDeleteRowsResourceWithStreamingResponse(self)

    async def delete(
        self,
        *,
        row_ids: List[str],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> DeleteRowDeleteResponse:
        """
        Delete rows by their ids.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/DeleteRows",
            body=await async_maybe_transform({"row_ids": row_ids}, delete_row_delete_params.DeleteRowDeleteParams),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=DeleteRowDeleteResponse,
        )


class DeleteRowsResourceWithRawResponse:
    def __init__(self, delete_rows: DeleteRowsResource) -> None:
        self._delete_rows = delete_rows

        self.delete = to_raw_response_wrapper(
            delete_rows.delete,
        )


class AsyncDeleteRowsResourceWithRawResponse:
    def __init__(self, delete_rows: AsyncDeleteRowsResource) -> None:
        self._delete_rows = delete_rows

        self.delete = async_to_raw_response_wrapper(
            delete_rows.delete,
        )


class DeleteRowsResourceWithStreamingResponse:
    def __init__(self, delete_rows: DeleteRowsResource) -> None:
        self._delete_rows = delete_rows

        self.delete = to_streamed_response_wrapper(
            delete_rows.delete,
        )


class AsyncDeleteRowsResourceWithStreamingResponse:
    def __init__(self, delete_rows: AsyncDeleteRowsResource) -> None:
        self._delete_rows = delete_rows

        self.delete = async_to_streamed_response_wrapper(
            delete_rows.delete,
        )
