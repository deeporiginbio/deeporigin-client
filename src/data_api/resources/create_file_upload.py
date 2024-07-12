# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Optional

import httpx

from ..types import create_file_upload_create_params
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
from ..types.create_file_upload_create_response import CreateFileUploadCreateResponse

__all__ = ["CreateFileUploadResource", "AsyncCreateFileUploadResource"]


class CreateFileUploadResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> CreateFileUploadResourceWithRawResponse:
        return CreateFileUploadResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> CreateFileUploadResourceWithStreamingResponse:
        return CreateFileUploadResourceWithStreamingResponse(self)

    def create(
        self,
        *,
        content_length: str,
        content_type: Optional[str],
        name: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> CreateFileUploadCreateResponse:
        """Create a file upload URL.

        Typically this is creating a pre-signed S3 URL.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/CreateFileUpload",
            body=maybe_transform(
                {
                    "content_length": content_length,
                    "content_type": content_type,
                    "name": name,
                },
                create_file_upload_create_params.CreateFileUploadCreateParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=CreateFileUploadCreateResponse,
        )


class AsyncCreateFileUploadResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncCreateFileUploadResourceWithRawResponse:
        return AsyncCreateFileUploadResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncCreateFileUploadResourceWithStreamingResponse:
        return AsyncCreateFileUploadResourceWithStreamingResponse(self)

    async def create(
        self,
        *,
        content_length: str,
        content_type: Optional[str],
        name: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> CreateFileUploadCreateResponse:
        """Create a file upload URL.

        Typically this is creating a pre-signed S3 URL.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/CreateFileUpload",
            body=await async_maybe_transform(
                {
                    "content_length": content_length,
                    "content_type": content_type,
                    "name": name,
                },
                create_file_upload_create_params.CreateFileUploadCreateParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=CreateFileUploadCreateResponse,
        )


class CreateFileUploadResourceWithRawResponse:
    def __init__(self, create_file_upload: CreateFileUploadResource) -> None:
        self._create_file_upload = create_file_upload

        self.create = to_raw_response_wrapper(
            create_file_upload.create,
        )


class AsyncCreateFileUploadResourceWithRawResponse:
    def __init__(self, create_file_upload: AsyncCreateFileUploadResource) -> None:
        self._create_file_upload = create_file_upload

        self.create = async_to_raw_response_wrapper(
            create_file_upload.create,
        )


class CreateFileUploadResourceWithStreamingResponse:
    def __init__(self, create_file_upload: CreateFileUploadResource) -> None:
        self._create_file_upload = create_file_upload

        self.create = to_streamed_response_wrapper(
            create_file_upload.create,
        )


class AsyncCreateFileUploadResourceWithStreamingResponse:
    def __init__(self, create_file_upload: AsyncCreateFileUploadResource) -> None:
        self._create_file_upload = create_file_upload

        self.create = async_to_streamed_response_wrapper(
            create_file_upload.create,
        )
