# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..types import delete_database_column_delete_params
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
from ..types.delete_database_column_delete_response import DeleteDatabaseColumnDeleteResponse

__all__ = ["DeleteDatabaseColumnResource", "AsyncDeleteDatabaseColumnResource"]


class DeleteDatabaseColumnResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> DeleteDatabaseColumnResourceWithRawResponse:
        return DeleteDatabaseColumnResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> DeleteDatabaseColumnResourceWithStreamingResponse:
        return DeleteDatabaseColumnResourceWithStreamingResponse(self)

    def delete(
        self,
        *,
        column_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> DeleteDatabaseColumnDeleteResponse:
        """
        Delete a column from a database.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/DeleteDatabaseColumn",
            body=maybe_transform(
                {"column_id": column_id}, delete_database_column_delete_params.DeleteDatabaseColumnDeleteParams
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=DeleteDatabaseColumnDeleteResponse,
        )


class AsyncDeleteDatabaseColumnResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncDeleteDatabaseColumnResourceWithRawResponse:
        return AsyncDeleteDatabaseColumnResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncDeleteDatabaseColumnResourceWithStreamingResponse:
        return AsyncDeleteDatabaseColumnResourceWithStreamingResponse(self)

    async def delete(
        self,
        *,
        column_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> DeleteDatabaseColumnDeleteResponse:
        """
        Delete a column from a database.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/DeleteDatabaseColumn",
            body=await async_maybe_transform(
                {"column_id": column_id}, delete_database_column_delete_params.DeleteDatabaseColumnDeleteParams
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=DeleteDatabaseColumnDeleteResponse,
        )


class DeleteDatabaseColumnResourceWithRawResponse:
    def __init__(self, delete_database_column: DeleteDatabaseColumnResource) -> None:
        self._delete_database_column = delete_database_column

        self.delete = to_raw_response_wrapper(
            delete_database_column.delete,
        )


class AsyncDeleteDatabaseColumnResourceWithRawResponse:
    def __init__(self, delete_database_column: AsyncDeleteDatabaseColumnResource) -> None:
        self._delete_database_column = delete_database_column

        self.delete = async_to_raw_response_wrapper(
            delete_database_column.delete,
        )


class DeleteDatabaseColumnResourceWithStreamingResponse:
    def __init__(self, delete_database_column: DeleteDatabaseColumnResource) -> None:
        self._delete_database_column = delete_database_column

        self.delete = to_streamed_response_wrapper(
            delete_database_column.delete,
        )


class AsyncDeleteDatabaseColumnResourceWithStreamingResponse:
    def __init__(self, delete_database_column: AsyncDeleteDatabaseColumnResource) -> None:
        self._delete_database_column = delete_database_column

        self.delete = async_to_streamed_response_wrapper(
            delete_database_column.delete,
        )
