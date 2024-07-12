# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..types import create_file_download_url_create_params
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
from ..types.create_file_download_url_create_response import CreateFileDownloadURLCreateResponse

__all__ = ["CreateFileDownloadURLResource", "AsyncCreateFileDownloadURLResource"]


class CreateFileDownloadURLResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> CreateFileDownloadURLResourceWithRawResponse:
        return CreateFileDownloadURLResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> CreateFileDownloadURLResourceWithStreamingResponse:
        return CreateFileDownloadURLResourceWithStreamingResponse(self)

    def create(
        self,
        *,
        file_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> CreateFileDownloadURLCreateResponse:
        """
        Returns a pre-signed S3 URL.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/CreateFileDownloadUrl",
            body=maybe_transform(
                {"file_id": file_id}, create_file_download_url_create_params.CreateFileDownloadURLCreateParams
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=CreateFileDownloadURLCreateResponse,
        )


class AsyncCreateFileDownloadURLResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncCreateFileDownloadURLResourceWithRawResponse:
        return AsyncCreateFileDownloadURLResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncCreateFileDownloadURLResourceWithStreamingResponse:
        return AsyncCreateFileDownloadURLResourceWithStreamingResponse(self)

    async def create(
        self,
        *,
        file_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> CreateFileDownloadURLCreateResponse:
        """
        Returns a pre-signed S3 URL.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/CreateFileDownloadUrl",
            body=await async_maybe_transform(
                {"file_id": file_id}, create_file_download_url_create_params.CreateFileDownloadURLCreateParams
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=CreateFileDownloadURLCreateResponse,
        )


class CreateFileDownloadURLResourceWithRawResponse:
    def __init__(self, create_file_download_url: CreateFileDownloadURLResource) -> None:
        self._create_file_download_url = create_file_download_url

        self.create = to_raw_response_wrapper(
            create_file_download_url.create,
        )


class AsyncCreateFileDownloadURLResourceWithRawResponse:
    def __init__(self, create_file_download_url: AsyncCreateFileDownloadURLResource) -> None:
        self._create_file_download_url = create_file_download_url

        self.create = async_to_raw_response_wrapper(
            create_file_download_url.create,
        )


class CreateFileDownloadURLResourceWithStreamingResponse:
    def __init__(self, create_file_download_url: CreateFileDownloadURLResource) -> None:
        self._create_file_download_url = create_file_download_url

        self.create = to_streamed_response_wrapper(
            create_file_download_url.create,
        )


class AsyncCreateFileDownloadURLResourceWithStreamingResponse:
    def __init__(self, create_file_download_url: AsyncCreateFileDownloadURLResource) -> None:
        self._create_file_download_url = create_file_download_url

        self.create = async_to_streamed_response_wrapper(
            create_file_download_url.create,
        )
