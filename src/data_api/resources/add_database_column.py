# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..types import add_database_column_add_params
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
from ..types.add_database_column_add_response import AddDatabaseColumnAddResponse

__all__ = ["AddDatabaseColumnResource", "AsyncAddDatabaseColumnResource"]


class AddDatabaseColumnResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AddDatabaseColumnResourceWithRawResponse:
        return AddDatabaseColumnResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AddDatabaseColumnResourceWithStreamingResponse:
        return AddDatabaseColumnResourceWithStreamingResponse(self)

    def add(
        self,
        *,
        column: add_database_column_add_params.Column,
        database_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> AddDatabaseColumnAddResponse:
        """
        Add a column to a database.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/AddDatabaseColumn",
            body=maybe_transform(
                {
                    "column": column,
                    "database_id": database_id,
                },
                add_database_column_add_params.AddDatabaseColumnAddParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=AddDatabaseColumnAddResponse,
        )


class AsyncAddDatabaseColumnResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncAddDatabaseColumnResourceWithRawResponse:
        return AsyncAddDatabaseColumnResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncAddDatabaseColumnResourceWithStreamingResponse:
        return AsyncAddDatabaseColumnResourceWithStreamingResponse(self)

    async def add(
        self,
        *,
        column: add_database_column_add_params.Column,
        database_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> AddDatabaseColumnAddResponse:
        """
        Add a column to a database.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/AddDatabaseColumn",
            body=await async_maybe_transform(
                {
                    "column": column,
                    "database_id": database_id,
                },
                add_database_column_add_params.AddDatabaseColumnAddParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=AddDatabaseColumnAddResponse,
        )


class AddDatabaseColumnResourceWithRawResponse:
    def __init__(self, add_database_column: AddDatabaseColumnResource) -> None:
        self._add_database_column = add_database_column

        self.add = to_raw_response_wrapper(
            add_database_column.add,
        )


class AsyncAddDatabaseColumnResourceWithRawResponse:
    def __init__(self, add_database_column: AsyncAddDatabaseColumnResource) -> None:
        self._add_database_column = add_database_column

        self.add = async_to_raw_response_wrapper(
            add_database_column.add,
        )


class AddDatabaseColumnResourceWithStreamingResponse:
    def __init__(self, add_database_column: AddDatabaseColumnResource) -> None:
        self._add_database_column = add_database_column

        self.add = to_streamed_response_wrapper(
            add_database_column.add,
        )


class AsyncAddDatabaseColumnResourceWithStreamingResponse:
    def __init__(self, add_database_column: AsyncAddDatabaseColumnResource) -> None:
        self._add_database_column = add_database_column

        self.add = async_to_streamed_response_wrapper(
            add_database_column.add,
        )
