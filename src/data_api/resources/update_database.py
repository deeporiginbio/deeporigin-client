# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..types import update_database_run_params
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
from ..types.update_database_run_response import UpdateDatabaseRunResponse

__all__ = ["UpdateDatabaseResource", "AsyncUpdateDatabaseResource"]


class UpdateDatabaseResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> UpdateDatabaseResourceWithRawResponse:
        return UpdateDatabaseResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> UpdateDatabaseResourceWithStreamingResponse:
        return UpdateDatabaseResourceWithStreamingResponse(self)

    def run(
        self,
        *,
        id: str,
        database: update_database_run_params.Database,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> UpdateDatabaseRunResponse:
        """
        Update a database.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/UpdateDatabase",
            body=maybe_transform(
                {
                    "id": id,
                    "database": database,
                },
                update_database_run_params.UpdateDatabaseRunParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=UpdateDatabaseRunResponse,
        )


class AsyncUpdateDatabaseResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncUpdateDatabaseResourceWithRawResponse:
        return AsyncUpdateDatabaseResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(
        self,
    ) -> AsyncUpdateDatabaseResourceWithStreamingResponse:
        return AsyncUpdateDatabaseResourceWithStreamingResponse(self)

    async def run(
        self,
        *,
        id: str,
        database: update_database_run_params.Database,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> UpdateDatabaseRunResponse:
        """
        Update a database.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/UpdateDatabase",
            body=await async_maybe_transform(
                {
                    "id": id,
                    "database": database,
                },
                update_database_run_params.UpdateDatabaseRunParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=UpdateDatabaseRunResponse,
        )


class UpdateDatabaseResourceWithRawResponse:
    def __init__(self, update_database: UpdateDatabaseResource) -> None:
        self._update_database = update_database

        self.run = to_raw_response_wrapper(
            update_database.run,
        )


class AsyncUpdateDatabaseResourceWithRawResponse:
    def __init__(self, update_database: AsyncUpdateDatabaseResource) -> None:
        self._update_database = update_database

        self.run = async_to_raw_response_wrapper(
            update_database.run,
        )


class UpdateDatabaseResourceWithStreamingResponse:
    def __init__(self, update_database: UpdateDatabaseResource) -> None:
        self._update_database = update_database

        self.run = to_streamed_response_wrapper(
            update_database.run,
        )


class AsyncUpdateDatabaseResourceWithStreamingResponse:
    def __init__(self, update_database: AsyncUpdateDatabaseResource) -> None:
        self._update_database = update_database

        self.run = async_to_streamed_response_wrapper(
            update_database.run,
        )
