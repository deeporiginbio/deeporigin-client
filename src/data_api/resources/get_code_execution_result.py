# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..types import get_code_execution_result_retrieve_params
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

__all__ = ["GetCodeExecutionResultResource", "AsyncGetCodeExecutionResultResource"]


class GetCodeExecutionResultResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> GetCodeExecutionResultResourceWithRawResponse:
        return GetCodeExecutionResultResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> GetCodeExecutionResultResourceWithStreamingResponse:
        return GetCodeExecutionResultResourceWithStreamingResponse(self)

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
    ) -> object:
        """
        Returns the result of a code execution.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/GetCodeExecutionResult",
            body=maybe_transform(
                {"id": id}, get_code_execution_result_retrieve_params.GetCodeExecutionResultRetrieveParams
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=object,
        )


class AsyncGetCodeExecutionResultResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncGetCodeExecutionResultResourceWithRawResponse:
        return AsyncGetCodeExecutionResultResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncGetCodeExecutionResultResourceWithStreamingResponse:
        return AsyncGetCodeExecutionResultResourceWithStreamingResponse(self)

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
    ) -> object:
        """
        Returns the result of a code execution.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/GetCodeExecutionResult",
            body=await async_maybe_transform(
                {"id": id}, get_code_execution_result_retrieve_params.GetCodeExecutionResultRetrieveParams
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=object,
        )


class GetCodeExecutionResultResourceWithRawResponse:
    def __init__(self, get_code_execution_result: GetCodeExecutionResultResource) -> None:
        self._get_code_execution_result = get_code_execution_result

        self.retrieve = to_raw_response_wrapper(
            get_code_execution_result.retrieve,
        )


class AsyncGetCodeExecutionResultResourceWithRawResponse:
    def __init__(self, get_code_execution_result: AsyncGetCodeExecutionResultResource) -> None:
        self._get_code_execution_result = get_code_execution_result

        self.retrieve = async_to_raw_response_wrapper(
            get_code_execution_result.retrieve,
        )


class GetCodeExecutionResultResourceWithStreamingResponse:
    def __init__(self, get_code_execution_result: GetCodeExecutionResultResource) -> None:
        self._get_code_execution_result = get_code_execution_result

        self.retrieve = to_streamed_response_wrapper(
            get_code_execution_result.retrieve,
        )


class AsyncGetCodeExecutionResultResourceWithStreamingResponse:
    def __init__(self, get_code_execution_result: AsyncGetCodeExecutionResultResource) -> None:
        self._get_code_execution_result = get_code_execution_result

        self.retrieve = async_to_streamed_response_wrapper(
            get_code_execution_result.retrieve,
        )
