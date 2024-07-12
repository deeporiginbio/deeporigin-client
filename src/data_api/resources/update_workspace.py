# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..types import update_workspace_run_params
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
from ..types.update_workspace_run_response import UpdateWorkspaceRunResponse

__all__ = ["UpdateWorkspaceResource", "AsyncUpdateWorkspaceResource"]


class UpdateWorkspaceResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> UpdateWorkspaceResourceWithRawResponse:
        return UpdateWorkspaceResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> UpdateWorkspaceResourceWithStreamingResponse:
        return UpdateWorkspaceResourceWithStreamingResponse(self)

    def run(
        self,
        *,
        id: str,
        workspace: update_workspace_run_params.Workspace,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> UpdateWorkspaceRunResponse:
        """
        Update a workspace.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/UpdateWorkspace",
            body=maybe_transform(
                {
                    "id": id,
                    "workspace": workspace,
                },
                update_workspace_run_params.UpdateWorkspaceRunParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=UpdateWorkspaceRunResponse,
        )


class AsyncUpdateWorkspaceResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncUpdateWorkspaceResourceWithRawResponse:
        return AsyncUpdateWorkspaceResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncUpdateWorkspaceResourceWithStreamingResponse:
        return AsyncUpdateWorkspaceResourceWithStreamingResponse(self)

    async def run(
        self,
        *,
        id: str,
        workspace: update_workspace_run_params.Workspace,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> UpdateWorkspaceRunResponse:
        """
        Update a workspace.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/UpdateWorkspace",
            body=await async_maybe_transform(
                {
                    "id": id,
                    "workspace": workspace,
                },
                update_workspace_run_params.UpdateWorkspaceRunParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=UpdateWorkspaceRunResponse,
        )


class UpdateWorkspaceResourceWithRawResponse:
    def __init__(self, update_workspace: UpdateWorkspaceResource) -> None:
        self._update_workspace = update_workspace

        self.run = to_raw_response_wrapper(
            update_workspace.run,
        )


class AsyncUpdateWorkspaceResourceWithRawResponse:
    def __init__(self, update_workspace: AsyncUpdateWorkspaceResource) -> None:
        self._update_workspace = update_workspace

        self.run = async_to_raw_response_wrapper(
            update_workspace.run,
        )


class UpdateWorkspaceResourceWithStreamingResponse:
    def __init__(self, update_workspace: UpdateWorkspaceResource) -> None:
        self._update_workspace = update_workspace

        self.run = to_streamed_response_wrapper(
            update_workspace.run,
        )


class AsyncUpdateWorkspaceResourceWithStreamingResponse:
    def __init__(self, update_workspace: AsyncUpdateWorkspaceResource) -> None:
        self._update_workspace = update_workspace

        self.run = async_to_streamed_response_wrapper(
            update_workspace.run,
        )
