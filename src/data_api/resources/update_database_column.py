# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..types import update_database_column_run_params
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
from ..types.update_database_column_run_response import UpdateDatabaseColumnRunResponse

__all__ = ["UpdateDatabaseColumnResource", "AsyncUpdateDatabaseColumnResource"]


class UpdateDatabaseColumnResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> UpdateDatabaseColumnResourceWithRawResponse:
        return UpdateDatabaseColumnResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(
        self,
    ) -> UpdateDatabaseColumnResourceWithStreamingResponse:
        return UpdateDatabaseColumnResourceWithStreamingResponse(self)

    def run(
        self,
        *,
        column: update_database_column_run_params.Column,
        column_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> UpdateDatabaseColumnRunResponse:
        """
        Update a column in a database.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/UpdateDatabaseColumn",
            body=maybe_transform(
                {
                    "column": column,
                    "column_id": column_id,
                },
                update_database_column_run_params.UpdateDatabaseColumnRunParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=UpdateDatabaseColumnRunResponse,
        )


class AsyncUpdateDatabaseColumnResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncUpdateDatabaseColumnResourceWithRawResponse:
        return AsyncUpdateDatabaseColumnResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(
        self,
    ) -> AsyncUpdateDatabaseColumnResourceWithStreamingResponse:
        return AsyncUpdateDatabaseColumnResourceWithStreamingResponse(self)

    async def run(
        self,
        *,
        column: update_database_column_run_params.Column,
        column_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> UpdateDatabaseColumnRunResponse:
        """
        Update a column in a database.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/UpdateDatabaseColumn",
            body=await async_maybe_transform(
                {
                    "column": column,
                    "column_id": column_id,
                },
                update_database_column_run_params.UpdateDatabaseColumnRunParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=UpdateDatabaseColumnRunResponse,
        )


class UpdateDatabaseColumnResourceWithRawResponse:
    def __init__(self, update_database_column: UpdateDatabaseColumnResource) -> None:
        self._update_database_column = update_database_column

        self.run = to_raw_response_wrapper(
            update_database_column.run,
        )


class AsyncUpdateDatabaseColumnResourceWithRawResponse:
    def __init__(
        self, update_database_column: AsyncUpdateDatabaseColumnResource
    ) -> None:
        self._update_database_column = update_database_column

        self.run = async_to_raw_response_wrapper(
            update_database_column.run,
        )


class UpdateDatabaseColumnResourceWithStreamingResponse:
    def __init__(self, update_database_column: UpdateDatabaseColumnResource) -> None:
        self._update_database_column = update_database_column

        self.run = to_streamed_response_wrapper(
            update_database_column.run,
        )


class AsyncUpdateDatabaseColumnResourceWithStreamingResponse:
    def __init__(
        self, update_database_column: AsyncUpdateDatabaseColumnResource
    ) -> None:
        self._update_database_column = update_database_column

        self.run = async_to_streamed_response_wrapper(
            update_database_column.run,
        )
