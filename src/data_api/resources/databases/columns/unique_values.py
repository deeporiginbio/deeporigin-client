# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import overload
from typing_extensions import Literal

import httpx

from ...._types import NOT_GIVEN, Body, Query, Headers, NotGiven
from ...._utils import (
    required_args,
    maybe_transform,
    async_maybe_transform,
)
from ...._compat import cached_property
from ...._resource import SyncAPIResource, AsyncAPIResource
from ...._response import (
    to_raw_response_wrapper,
    to_streamed_response_wrapper,
    async_to_raw_response_wrapper,
    async_to_streamed_response_wrapper,
)
from ...._base_client import (
    make_request_options,
)
from ....types.databases.columns import unique_value_list_params
from ....types.databases.columns.unique_value_list_response import UniqueValueListResponse

__all__ = ["UniqueValuesResource", "AsyncUniqueValuesResource"]


class UniqueValuesResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> UniqueValuesResourceWithRawResponse:
        return UniqueValuesResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> UniqueValuesResourceWithStreamingResponse:
        return UniqueValuesResourceWithStreamingResponse(self)

    @overload
    def list(
        self,
        *,
        column_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> UniqueValueListResponse:
        """
        Returns the unique values for every cell within the column.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        ...

    @overload
    def list(
        self,
        *,
        database_row_id: str,
        system_column_name: Literal["creationParentId"],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> UniqueValueListResponse:
        """
        Returns the unique values for every cell within the column.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        ...

    @required_args(["column_id"], ["database_row_id", "system_column_name"])
    def list(
        self,
        *,
        column_id: str | NotGiven = NOT_GIVEN,
        database_row_id: str | NotGiven = NOT_GIVEN,
        system_column_name: Literal["creationParentId"] | NotGiven = NOT_GIVEN,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> UniqueValueListResponse:
        return self._post(
            "/ListDatabaseColumnUniqueValues",
            body=maybe_transform(
                {
                    "column_id": column_id,
                    "database_row_id": database_row_id,
                    "system_column_name": system_column_name,
                },
                unique_value_list_params.UniqueValueListParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=UniqueValueListResponse,
        )


class AsyncUniqueValuesResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncUniqueValuesResourceWithRawResponse:
        return AsyncUniqueValuesResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncUniqueValuesResourceWithStreamingResponse:
        return AsyncUniqueValuesResourceWithStreamingResponse(self)

    @overload
    async def list(
        self,
        *,
        column_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> UniqueValueListResponse:
        """
        Returns the unique values for every cell within the column.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        ...

    @overload
    async def list(
        self,
        *,
        database_row_id: str,
        system_column_name: Literal["creationParentId"],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> UniqueValueListResponse:
        """
        Returns the unique values for every cell within the column.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        ...

    @required_args(["column_id"], ["database_row_id", "system_column_name"])
    async def list(
        self,
        *,
        column_id: str | NotGiven = NOT_GIVEN,
        database_row_id: str | NotGiven = NOT_GIVEN,
        system_column_name: Literal["creationParentId"] | NotGiven = NOT_GIVEN,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> UniqueValueListResponse:
        return await self._post(
            "/ListDatabaseColumnUniqueValues",
            body=await async_maybe_transform(
                {
                    "column_id": column_id,
                    "database_row_id": database_row_id,
                    "system_column_name": system_column_name,
                },
                unique_value_list_params.UniqueValueListParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=UniqueValueListResponse,
        )


class UniqueValuesResourceWithRawResponse:
    def __init__(self, unique_values: UniqueValuesResource) -> None:
        self._unique_values = unique_values

        self.list = to_raw_response_wrapper(
            unique_values.list,
        )


class AsyncUniqueValuesResourceWithRawResponse:
    def __init__(self, unique_values: AsyncUniqueValuesResource) -> None:
        self._unique_values = unique_values

        self.list = async_to_raw_response_wrapper(
            unique_values.list,
        )


class UniqueValuesResourceWithStreamingResponse:
    def __init__(self, unique_values: UniqueValuesResource) -> None:
        self._unique_values = unique_values

        self.list = to_streamed_response_wrapper(
            unique_values.list,
        )


class AsyncUniqueValuesResourceWithStreamingResponse:
    def __init__(self, unique_values: AsyncUniqueValuesResource) -> None:
        self._unique_values = unique_values

        self.list = async_to_streamed_response_wrapper(
            unique_values.list,
        )
