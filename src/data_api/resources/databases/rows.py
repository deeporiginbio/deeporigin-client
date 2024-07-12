# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..._types import NOT_GIVEN, Body, Query, Headers, NotGiven
from ..._utils import (
    maybe_transform,
    async_maybe_transform,
)
from ..._compat import cached_property
from ..._resource import SyncAPIResource, AsyncAPIResource
from ..._response import (
    to_raw_response_wrapper,
    to_streamed_response_wrapper,
    async_to_raw_response_wrapper,
    async_to_streamed_response_wrapper,
)
from ..._base_client import (
    make_request_options,
)
from ...types.databases import row_list_params
from ...types.databases.row_list_response import RowListResponse

__all__ = ["RowsResource", "AsyncRowsResource"]


class RowsResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> RowsResourceWithRawResponse:
        return RowsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> RowsResourceWithStreamingResponse:
        return RowsResourceWithStreamingResponse(self)

    def list(
        self,
        *,
        database_row_id: str,
        column_selection: row_list_params.ColumnSelection | NotGiven = NOT_GIVEN,
        creation_block_id: str | NotGiven = NOT_GIVEN,
        creation_parent_id: str | NotGiven = NOT_GIVEN,
        filter: row_list_params.Filter | NotGiven = NOT_GIVEN,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> RowListResponse:
        """
        List database rows with full row data.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/ListDatabaseRows",
            body=maybe_transform(
                {
                    "database_row_id": database_row_id,
                    "column_selection": column_selection,
                    "creation_block_id": creation_block_id,
                    "creation_parent_id": creation_parent_id,
                    "filter": filter,
                },
                row_list_params.RowListParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=RowListResponse,
        )


class AsyncRowsResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncRowsResourceWithRawResponse:
        return AsyncRowsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncRowsResourceWithStreamingResponse:
        return AsyncRowsResourceWithStreamingResponse(self)

    async def list(
        self,
        *,
        database_row_id: str,
        column_selection: row_list_params.ColumnSelection | NotGiven = NOT_GIVEN,
        creation_block_id: str | NotGiven = NOT_GIVEN,
        creation_parent_id: str | NotGiven = NOT_GIVEN,
        filter: row_list_params.Filter | NotGiven = NOT_GIVEN,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> RowListResponse:
        """
        List database rows with full row data.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/ListDatabaseRows",
            body=await async_maybe_transform(
                {
                    "database_row_id": database_row_id,
                    "column_selection": column_selection,
                    "creation_block_id": creation_block_id,
                    "creation_parent_id": creation_parent_id,
                    "filter": filter,
                },
                row_list_params.RowListParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=RowListResponse,
        )


class RowsResourceWithRawResponse:
    def __init__(self, rows: RowsResource) -> None:
        self._rows = rows

        self.list = to_raw_response_wrapper(
            rows.list,
        )


class AsyncRowsResourceWithRawResponse:
    def __init__(self, rows: AsyncRowsResource) -> None:
        self._rows = rows

        self.list = async_to_raw_response_wrapper(
            rows.list,
        )


class RowsResourceWithStreamingResponse:
    def __init__(self, rows: RowsResource) -> None:
        self._rows = rows

        self.list = to_streamed_response_wrapper(
            rows.list,
        )


class AsyncRowsResourceWithStreamingResponse:
    def __init__(self, rows: AsyncRowsResource) -> None:
        self._rows = rows

        self.list = async_to_streamed_response_wrapper(
            rows.list,
        )
