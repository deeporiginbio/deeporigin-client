# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing_extensions import Literal

import httpx

from ..types import execute_code_execute_params
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
from ..types.execute_code_execute_response import ExecuteCodeExecuteResponse

__all__ = ["ExecuteCodeResource", "AsyncExecuteCodeResource"]


class ExecuteCodeResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> ExecuteCodeResourceWithRawResponse:
        return ExecuteCodeResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> ExecuteCodeResourceWithStreamingResponse:
        return ExecuteCodeResourceWithStreamingResponse(self)

    def execute(
        self,
        *,
        code: str,
        code_language: Literal["python"],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> ExecuteCodeExecuteResponse:
        """
        Execute code asynchronously.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/ExecuteCode",
            body=maybe_transform(
                {
                    "code": code,
                    "code_language": code_language,
                },
                execute_code_execute_params.ExecuteCodeExecuteParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=ExecuteCodeExecuteResponse,
        )


class AsyncExecuteCodeResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncExecuteCodeResourceWithRawResponse:
        return AsyncExecuteCodeResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncExecuteCodeResourceWithStreamingResponse:
        return AsyncExecuteCodeResourceWithStreamingResponse(self)

    async def execute(
        self,
        *,
        code: str,
        code_language: Literal["python"],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> ExecuteCodeExecuteResponse:
        """
        Execute code asynchronously.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/ExecuteCode",
            body=await async_maybe_transform(
                {
                    "code": code,
                    "code_language": code_language,
                },
                execute_code_execute_params.ExecuteCodeExecuteParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=ExecuteCodeExecuteResponse,
        )


class ExecuteCodeResourceWithRawResponse:
    def __init__(self, execute_code: ExecuteCodeResource) -> None:
        self._execute_code = execute_code

        self.execute = to_raw_response_wrapper(
            execute_code.execute,
        )


class AsyncExecuteCodeResourceWithRawResponse:
    def __init__(self, execute_code: AsyncExecuteCodeResource) -> None:
        self._execute_code = execute_code

        self.execute = async_to_raw_response_wrapper(
            execute_code.execute,
        )


class ExecuteCodeResourceWithStreamingResponse:
    def __init__(self, execute_code: ExecuteCodeResource) -> None:
        self._execute_code = execute_code

        self.execute = to_streamed_response_wrapper(
            execute_code.execute,
        )


class AsyncExecuteCodeResourceWithStreamingResponse:
    def __init__(self, execute_code: AsyncExecuteCodeResource) -> None:
        self._execute_code = execute_code

        self.execute = async_to_streamed_response_wrapper(
            execute_code.execute,
        )
