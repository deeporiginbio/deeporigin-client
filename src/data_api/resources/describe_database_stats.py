# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..types import describe_database_stat_retrieve_params
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
from ..types.describe_database_stat_retrieve_response import (
    DescribeDatabaseStatRetrieveResponse,
)

__all__ = ["DescribeDatabaseStatsResource", "AsyncDescribeDatabaseStatsResource"]


class DescribeDatabaseStatsResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> DescribeDatabaseStatsResourceWithRawResponse:
        return DescribeDatabaseStatsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(
        self,
    ) -> DescribeDatabaseStatsResourceWithStreamingResponse:
        return DescribeDatabaseStatsResourceWithStreamingResponse(self)

    def retrieve(
        self,
        *,
        database_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> DescribeDatabaseStatRetrieveResponse:
        """
        Returns aggregation information about a particular database.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/DescribeDatabaseStats",
            body=maybe_transform(
                {"database_id": database_id},
                describe_database_stat_retrieve_params.DescribeDatabaseStatRetrieveParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=DescribeDatabaseStatRetrieveResponse,
        )


class AsyncDescribeDatabaseStatsResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncDescribeDatabaseStatsResourceWithRawResponse:
        return AsyncDescribeDatabaseStatsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(
        self,
    ) -> AsyncDescribeDatabaseStatsResourceWithStreamingResponse:
        return AsyncDescribeDatabaseStatsResourceWithStreamingResponse(self)

    async def retrieve(
        self,
        *,
        database_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> DescribeDatabaseStatRetrieveResponse:
        """
        Returns aggregation information about a particular database.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/DescribeDatabaseStats",
            body=await async_maybe_transform(
                {"database_id": database_id},
                describe_database_stat_retrieve_params.DescribeDatabaseStatRetrieveParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=DescribeDatabaseStatRetrieveResponse,
        )


class DescribeDatabaseStatsResourceWithRawResponse:
    def __init__(self, describe_database_stats: DescribeDatabaseStatsResource) -> None:
        self._describe_database_stats = describe_database_stats

        self.retrieve = to_raw_response_wrapper(
            describe_database_stats.retrieve,
        )


class AsyncDescribeDatabaseStatsResourceWithRawResponse:
    def __init__(
        self, describe_database_stats: AsyncDescribeDatabaseStatsResource
    ) -> None:
        self._describe_database_stats = describe_database_stats

        self.retrieve = async_to_raw_response_wrapper(
            describe_database_stats.retrieve,
        )


class DescribeDatabaseStatsResourceWithStreamingResponse:
    def __init__(self, describe_database_stats: DescribeDatabaseStatsResource) -> None:
        self._describe_database_stats = describe_database_stats

        self.retrieve = to_streamed_response_wrapper(
            describe_database_stats.retrieve,
        )


class AsyncDescribeDatabaseStatsResourceWithStreamingResponse:
    def __init__(
        self, describe_database_stats: AsyncDescribeDatabaseStatsResource
    ) -> None:
        self._describe_database_stats = describe_database_stats

        self.retrieve = async_to_streamed_response_wrapper(
            describe_database_stats.retrieve,
        )
