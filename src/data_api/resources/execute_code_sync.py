# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing_extensions import Literal

import httpx

from ..types import execute_code_sync_execute_sync_params
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
from ..types.execute_code_sync_execute_sync_response import ExecuteCodeSyncExecuteSyncResponse

__all__ = ["ExecuteCodeSyncResource", "AsyncExecuteCodeSyncResource"]


class ExecuteCodeSyncResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> ExecuteCodeSyncResourceWithRawResponse:
        return ExecuteCodeSyncResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> ExecuteCodeSyncResourceWithStreamingResponse:
        return ExecuteCodeSyncResourceWithStreamingResponse(self)

    def execute_sync(
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
    ) -> ExecuteCodeSyncExecuteSyncResponse:
        """
        Execute code synchronously.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/ExecuteCodeSync",
            body=maybe_transform(
                {
                    "code": code,
                    "code_language": code_language,
                },
                execute_code_sync_execute_sync_params.ExecuteCodeSyncExecuteSyncParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=ExecuteCodeSyncExecuteSyncResponse,
        )


class AsyncExecuteCodeSyncResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncExecuteCodeSyncResourceWithRawResponse:
        return AsyncExecuteCodeSyncResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncExecuteCodeSyncResourceWithStreamingResponse:
        return AsyncExecuteCodeSyncResourceWithStreamingResponse(self)

    async def execute_sync(
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
    ) -> ExecuteCodeSyncExecuteSyncResponse:
        """
        Execute code synchronously.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/ExecuteCodeSync",
            body=await async_maybe_transform(
                {
                    "code": code,
                    "code_language": code_language,
                },
                execute_code_sync_execute_sync_params.ExecuteCodeSyncExecuteSyncParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=ExecuteCodeSyncExecuteSyncResponse,
        )


class ExecuteCodeSyncResourceWithRawResponse:
    def __init__(self, execute_code_sync: ExecuteCodeSyncResource) -> None:
        self._execute_code_sync = execute_code_sync

        self.execute_sync = to_raw_response_wrapper(
            execute_code_sync.execute_sync,
        )


class AsyncExecuteCodeSyncResourceWithRawResponse:
    def __init__(self, execute_code_sync: AsyncExecuteCodeSyncResource) -> None:
        self._execute_code_sync = execute_code_sync

        self.execute_sync = async_to_raw_response_wrapper(
            execute_code_sync.execute_sync,
        )


class ExecuteCodeSyncResourceWithStreamingResponse:
    def __init__(self, execute_code_sync: ExecuteCodeSyncResource) -> None:
        self._execute_code_sync = execute_code_sync

        self.execute_sync = to_streamed_response_wrapper(
            execute_code_sync.execute_sync,
        )


class AsyncExecuteCodeSyncResourceWithStreamingResponse:
    def __init__(self, execute_code_sync: AsyncExecuteCodeSyncResource) -> None:
        self._execute_code_sync = execute_code_sync

        self.execute_sync = async_to_streamed_response_wrapper(
            execute_code_sync.execute_sync,
        )
