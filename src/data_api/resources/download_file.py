# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from .._types import NOT_GIVEN, Body, Query, Headers, NoneType, NotGiven
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

__all__ = ["DownloadFileResource", "AsyncDownloadFileResource"]


class DownloadFileResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> DownloadFileResourceWithRawResponse:
        return DownloadFileResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> DownloadFileResourceWithStreamingResponse:
        return DownloadFileResourceWithStreamingResponse(self)

    def redirect(
        self,
        *,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> None:
        """Returns a 303 redirect to a pre-signed S3 URL."""
        extra_headers = {"Accept": "*/*", **(extra_headers or {})}
        return self._get(
            "/DownloadFile",
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=NoneType,
        )


class AsyncDownloadFileResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncDownloadFileResourceWithRawResponse:
        return AsyncDownloadFileResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncDownloadFileResourceWithStreamingResponse:
        return AsyncDownloadFileResourceWithStreamingResponse(self)

    async def redirect(
        self,
        *,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> None:
        """Returns a 303 redirect to a pre-signed S3 URL."""
        extra_headers = {"Accept": "*/*", **(extra_headers or {})}
        return await self._get(
            "/DownloadFile",
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=NoneType,
        )


class DownloadFileResourceWithRawResponse:
    def __init__(self, download_file: DownloadFileResource) -> None:
        self._download_file = download_file

        self.redirect = to_raw_response_wrapper(
            download_file.redirect,
        )


class AsyncDownloadFileResourceWithRawResponse:
    def __init__(self, download_file: AsyncDownloadFileResource) -> None:
        self._download_file = download_file

        self.redirect = async_to_raw_response_wrapper(
            download_file.redirect,
        )


class DownloadFileResourceWithStreamingResponse:
    def __init__(self, download_file: DownloadFileResource) -> None:
        self._download_file = download_file

        self.redirect = to_streamed_response_wrapper(
            download_file.redirect,
        )


class AsyncDownloadFileResourceWithStreamingResponse:
    def __init__(self, download_file: AsyncDownloadFileResource) -> None:
        self._download_file = download_file

        self.redirect = async_to_streamed_response_wrapper(
            download_file.redirect,
        )
