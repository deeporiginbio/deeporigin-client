# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..types import describe_code_execution_retrieve_params
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
from ..types.describe_code_execution_retrieve_response import DescribeCodeExecutionRetrieveResponse

__all__ = ["DescribeCodeExecutionResource", "AsyncDescribeCodeExecutionResource"]


class DescribeCodeExecutionResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> DescribeCodeExecutionResourceWithRawResponse:
        return DescribeCodeExecutionResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> DescribeCodeExecutionResourceWithStreamingResponse:
        return DescribeCodeExecutionResourceWithStreamingResponse(self)

    def retrieve(
        self,
        *,
        id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> DescribeCodeExecutionRetrieveResponse:
        """
        Returns information about a particular code execution.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/DescribeCodeExecution",
            body=maybe_transform(
                {"id": id}, describe_code_execution_retrieve_params.DescribeCodeExecutionRetrieveParams
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=DescribeCodeExecutionRetrieveResponse,
        )


class AsyncDescribeCodeExecutionResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncDescribeCodeExecutionResourceWithRawResponse:
        return AsyncDescribeCodeExecutionResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncDescribeCodeExecutionResourceWithStreamingResponse:
        return AsyncDescribeCodeExecutionResourceWithStreamingResponse(self)

    async def retrieve(
        self,
        *,
        id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> DescribeCodeExecutionRetrieveResponse:
        """
        Returns information about a particular code execution.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/DescribeCodeExecution",
            body=await async_maybe_transform(
                {"id": id}, describe_code_execution_retrieve_params.DescribeCodeExecutionRetrieveParams
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=DescribeCodeExecutionRetrieveResponse,
        )


class DescribeCodeExecutionResourceWithRawResponse:
    def __init__(self, describe_code_execution: DescribeCodeExecutionResource) -> None:
        self._describe_code_execution = describe_code_execution

        self.retrieve = to_raw_response_wrapper(
            describe_code_execution.retrieve,
        )


class AsyncDescribeCodeExecutionResourceWithRawResponse:
    def __init__(self, describe_code_execution: AsyncDescribeCodeExecutionResource) -> None:
        self._describe_code_execution = describe_code_execution

        self.retrieve = async_to_raw_response_wrapper(
            describe_code_execution.retrieve,
        )


class DescribeCodeExecutionResourceWithStreamingResponse:
    def __init__(self, describe_code_execution: DescribeCodeExecutionResource) -> None:
        self._describe_code_execution = describe_code_execution

        self.retrieve = to_streamed_response_wrapper(
            describe_code_execution.retrieve,
        )


class AsyncDescribeCodeExecutionResourceWithStreamingResponse:
    def __init__(self, describe_code_execution: AsyncDescribeCodeExecutionResource) -> None:
        self._describe_code_execution = describe_code_execution

        self.retrieve = async_to_streamed_response_wrapper(
            describe_code_execution.retrieve,
        )
