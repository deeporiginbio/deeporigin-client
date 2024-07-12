# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import os
from typing import Any, Dict, Union, Mapping, cast
from typing_extensions import Self, Literal, override

import httpx

from . import resources, _exceptions
from ._qs import Querystring
from ._types import (
    NOT_GIVEN,
    Omit,
    Timeout,
    NotGiven,
    Transport,
    ProxiesTypes,
    RequestOptions,
)
from ._utils import (
    is_given,
    get_async_library,
)
from ._version import __version__
from ._streaming import Stream as Stream, AsyncStream as AsyncStream
from ._exceptions import APIStatusError
from ._base_client import (
    DEFAULT_MAX_RETRIES,
    SyncAPIClient,
    AsyncAPIClient,
)

__all__ = [
    "ENVIRONMENTS",
    "Timeout",
    "Transport",
    "ProxiesTypes",
    "RequestOptions",
    "resources",
    "DeeporiginData",
    "AsyncDeeporiginData",
    "Client",
    "AsyncClient",
]

ENVIRONMENTS: Dict[str, str] = {
    "edge": "https://os.edge.deeporigin.io/nucleus-api/api",
    "prod": "https://os.deeporigin.io/nucleus-api/api",
}


class DeeporiginData(SyncAPIClient):
    add_database_column: resources.AddDatabaseColumnResource
    archive_files: resources.ArchiveFilesResource
    configure_column_select_options: resources.ConfigureColumnSelectOptionsResource
    convert_id_format: resources.ConvertIDFormatResource
    databases: resources.DatabasesResource
    create_file_download_url: resources.CreateFileDownloadURLResource
    create_file_upload: resources.CreateFileUploadResource
    workspaces: resources.WorkspacesResource
    delete_database_column: resources.DeleteDatabaseColumnResource
    delete_rows: resources.DeleteRowsResource
    describe_code_execution: resources.DescribeCodeExecutionResource
    describe_database_stats: resources.DescribeDatabaseStatsResource
    describe_file: resources.DescribeFileResource
    describe_row: resources.DescribeRowResource
    download_file: resources.DownloadFileResource
    ensure_rows: resources.EnsureRowsResource
    execute_code: resources.ExecuteCodeResource
    execute_code_sync: resources.ExecuteCodeSyncResource
    get_code_execution_result: resources.GetCodeExecutionResultResource
    import_rows: resources.ImportRowsResource
    organization: resources.OrganizationResource
    chat_threads: resources.ChatThreadsResource
    files: resources.FilesResource
    mentions: resources.MentionsResource
    rows: resources.RowsResource
    sequences: resources.SequencesResource
    chat_messages: resources.ChatMessagesResource
    update_database: resources.UpdateDatabaseResource
    update_database_column: resources.UpdateDatabaseColumnResource
    update_workspace: resources.UpdateWorkspaceResource
    with_raw_response: DeeporiginDataWithRawResponse
    with_streaming_response: DeeporiginDataWithStreamedResponse

    # client options
    bearer_token: str
    org_id: str

    _environment: Literal["edge", "prod"] | NotGiven

    def __init__(
        self,
        *,
        bearer_token: str,
        org_id: str,
        environment: Literal["edge", "prod"] | NotGiven = NOT_GIVEN,
        base_url: str | httpx.URL | None | NotGiven = NOT_GIVEN,
        timeout: Union[float, Timeout, None, NotGiven] = NOT_GIVEN,
        max_retries: int = DEFAULT_MAX_RETRIES,
        default_headers: Mapping[str, str] | None = None,
        default_query: Mapping[str, object] | None = None,
        # Configure a custom httpx client.
        # We provide a `DefaultHttpxClient` class that you can pass to retain the default values we use for `limits`, `timeout` & `follow_redirects`.
        # See the [httpx documentation](https://www.python-httpx.org/api/#client) for more details.
        http_client: httpx.Client | None = None,
        # Enable or disable schema validation for data returned by the API.
        # When enabled an error APIResponseValidationError is raised
        # if the API responds with invalid data for the expected schema.
        #
        # This parameter may be removed or changed in the future.
        # If you rely on this feature, please open a GitHub issue
        # outlining your use-case to help us decide if it should be
        # part of our public interface in the future.
        _strict_response_validation: bool = False,
    ) -> None:
        """Construct a new synchronous deeporigin_data client instance."""
        self.bearer_token = bearer_token

        self.org_id = org_id

        self._environment = environment

        base_url_env = os.environ.get("DEEPORIGIN_DATA_BASE_URL")
        if is_given(base_url) and base_url is not None:
            # cast required because mypy doesn't understand the type narrowing
            base_url = cast("str | httpx.URL", base_url)  # pyright: ignore[reportUnnecessaryCast]
        elif is_given(environment):
            if base_url_env and base_url is not None:
                raise ValueError(
                    "Ambiguous URL; The `DEEPORIGIN_DATA_BASE_URL` env var and the `environment` argument are given. If you want to use the environment, you must pass base_url=None",
                )

            try:
                base_url = ENVIRONMENTS[environment]
            except KeyError as exc:
                raise ValueError(f"Unknown environment: {environment}") from exc
        elif base_url_env is not None:
            base_url = base_url_env
        else:
            self._environment = environment = "edge"

            try:
                base_url = ENVIRONMENTS[environment]
            except KeyError as exc:
                raise ValueError(f"Unknown environment: {environment}") from exc

        super().__init__(
            version=__version__,
            base_url=base_url,
            max_retries=max_retries,
            timeout=timeout,
            http_client=http_client,
            custom_headers=default_headers,
            custom_query=default_query,
            _strict_response_validation=_strict_response_validation,
        )

        self.add_database_column = resources.AddDatabaseColumnResource(self)
        self.archive_files = resources.ArchiveFilesResource(self)
        self.configure_column_select_options = resources.ConfigureColumnSelectOptionsResource(self)
        self.convert_id_format = resources.ConvertIDFormatResource(self)
        self.databases = resources.DatabasesResource(self)
        self.create_file_download_url = resources.CreateFileDownloadURLResource(self)
        self.create_file_upload = resources.CreateFileUploadResource(self)
        self.workspaces = resources.WorkspacesResource(self)
        self.delete_database_column = resources.DeleteDatabaseColumnResource(self)
        self.delete_rows = resources.DeleteRowsResource(self)
        self.describe_code_execution = resources.DescribeCodeExecutionResource(self)
        self.describe_database_stats = resources.DescribeDatabaseStatsResource(self)
        self.describe_file = resources.DescribeFileResource(self)
        self.describe_row = resources.DescribeRowResource(self)
        self.download_file = resources.DownloadFileResource(self)
        self.ensure_rows = resources.EnsureRowsResource(self)
        self.execute_code = resources.ExecuteCodeResource(self)
        self.execute_code_sync = resources.ExecuteCodeSyncResource(self)
        self.get_code_execution_result = resources.GetCodeExecutionResultResource(self)
        self.import_rows = resources.ImportRowsResource(self)
        self.organization = resources.OrganizationResource(self)
        self.chat_threads = resources.ChatThreadsResource(self)
        self.files = resources.FilesResource(self)
        self.mentions = resources.MentionsResource(self)
        self.rows = resources.RowsResource(self)
        self.sequences = resources.SequencesResource(self)
        self.chat_messages = resources.ChatMessagesResource(self)
        self.update_database = resources.UpdateDatabaseResource(self)
        self.update_database_column = resources.UpdateDatabaseColumnResource(self)
        self.update_workspace = resources.UpdateWorkspaceResource(self)
        self.with_raw_response = DeeporiginDataWithRawResponse(self)
        self.with_streaming_response = DeeporiginDataWithStreamedResponse(self)

    @property
    @override
    def qs(self) -> Querystring:
        return Querystring(array_format="comma")

    @property
    @override
    def auth_headers(self) -> dict[str, str]:
        bearer_token = self.bearer_token
        return {"Authorization": f"Bearer {bearer_token}"}

    @property
    @override
    def default_headers(self) -> dict[str, str | Omit]:
        return {
            **super().default_headers,
            "X-Stainless-Async": "false",
            "x-org-id": self.org_id,
            **self._custom_headers,
        }

    def copy(
        self,
        *,
        bearer_token: str | None = None,
        org_id: str | None = None,
        environment: Literal["edge", "prod"] | None = None,
        base_url: str | httpx.URL | None = None,
        timeout: float | Timeout | None | NotGiven = NOT_GIVEN,
        http_client: httpx.Client | None = None,
        max_retries: int | NotGiven = NOT_GIVEN,
        default_headers: Mapping[str, str] | None = None,
        set_default_headers: Mapping[str, str] | None = None,
        default_query: Mapping[str, object] | None = None,
        set_default_query: Mapping[str, object] | None = None,
        _extra_kwargs: Mapping[str, Any] = {},
    ) -> Self:
        """
        Create a new client instance re-using the same options given to the current client with optional overriding.
        """
        if default_headers is not None and set_default_headers is not None:
            raise ValueError("The `default_headers` and `set_default_headers` arguments are mutually exclusive")

        if default_query is not None and set_default_query is not None:
            raise ValueError("The `default_query` and `set_default_query` arguments are mutually exclusive")

        headers = self._custom_headers
        if default_headers is not None:
            headers = {**headers, **default_headers}
        elif set_default_headers is not None:
            headers = set_default_headers

        params = self._custom_query
        if default_query is not None:
            params = {**params, **default_query}
        elif set_default_query is not None:
            params = set_default_query

        http_client = http_client or self._client
        return self.__class__(
            bearer_token=bearer_token or self.bearer_token,
            org_id=org_id or self.org_id,
            base_url=base_url or self.base_url,
            environment=environment or self._environment,
            timeout=self.timeout if isinstance(timeout, NotGiven) else timeout,
            http_client=http_client,
            max_retries=max_retries if is_given(max_retries) else self.max_retries,
            default_headers=headers,
            default_query=params,
            **_extra_kwargs,
        )

    # Alias for `copy` for nicer inline usage, e.g.
    # client.with_options(timeout=10).foo.create(...)
    with_options = copy

    @override
    def _make_status_error(
        self,
        err_msg: str,
        *,
        body: object,
        response: httpx.Response,
    ) -> APIStatusError:
        if response.status_code == 400:
            return _exceptions.BadRequestError(err_msg, response=response, body=body)

        if response.status_code == 401:
            return _exceptions.AuthenticationError(err_msg, response=response, body=body)

        if response.status_code == 403:
            return _exceptions.PermissionDeniedError(err_msg, response=response, body=body)

        if response.status_code == 404:
            return _exceptions.NotFoundError(err_msg, response=response, body=body)

        if response.status_code == 409:
            return _exceptions.ConflictError(err_msg, response=response, body=body)

        if response.status_code == 422:
            return _exceptions.UnprocessableEntityError(err_msg, response=response, body=body)

        if response.status_code == 429:
            return _exceptions.RateLimitError(err_msg, response=response, body=body)

        if response.status_code >= 500:
            return _exceptions.InternalServerError(err_msg, response=response, body=body)
        return APIStatusError(err_msg, response=response, body=body)


class AsyncDeeporiginData(AsyncAPIClient):
    add_database_column: resources.AsyncAddDatabaseColumnResource
    archive_files: resources.AsyncArchiveFilesResource
    configure_column_select_options: resources.AsyncConfigureColumnSelectOptionsResource
    convert_id_format: resources.AsyncConvertIDFormatResource
    databases: resources.AsyncDatabasesResource
    create_file_download_url: resources.AsyncCreateFileDownloadURLResource
    create_file_upload: resources.AsyncCreateFileUploadResource
    workspaces: resources.AsyncWorkspacesResource
    delete_database_column: resources.AsyncDeleteDatabaseColumnResource
    delete_rows: resources.AsyncDeleteRowsResource
    describe_code_execution: resources.AsyncDescribeCodeExecutionResource
    describe_database_stats: resources.AsyncDescribeDatabaseStatsResource
    describe_file: resources.AsyncDescribeFileResource
    describe_row: resources.AsyncDescribeRowResource
    download_file: resources.AsyncDownloadFileResource
    ensure_rows: resources.AsyncEnsureRowsResource
    execute_code: resources.AsyncExecuteCodeResource
    execute_code_sync: resources.AsyncExecuteCodeSyncResource
    get_code_execution_result: resources.AsyncGetCodeExecutionResultResource
    import_rows: resources.AsyncImportRowsResource
    organization: resources.AsyncOrganizationResource
    chat_threads: resources.AsyncChatThreadsResource
    files: resources.AsyncFilesResource
    mentions: resources.AsyncMentionsResource
    rows: resources.AsyncRowsResource
    sequences: resources.AsyncSequencesResource
    chat_messages: resources.AsyncChatMessagesResource
    update_database: resources.AsyncUpdateDatabaseResource
    update_database_column: resources.AsyncUpdateDatabaseColumnResource
    update_workspace: resources.AsyncUpdateWorkspaceResource
    with_raw_response: AsyncDeeporiginDataWithRawResponse
    with_streaming_response: AsyncDeeporiginDataWithStreamedResponse

    # client options
    bearer_token: str
    org_id: str

    _environment: Literal["edge", "prod"] | NotGiven

    def __init__(
        self,
        *,
        bearer_token: str,
        org_id: str,
        environment: Literal["edge", "prod"] | NotGiven = NOT_GIVEN,
        base_url: str | httpx.URL | None | NotGiven = NOT_GIVEN,
        timeout: Union[float, Timeout, None, NotGiven] = NOT_GIVEN,
        max_retries: int = DEFAULT_MAX_RETRIES,
        default_headers: Mapping[str, str] | None = None,
        default_query: Mapping[str, object] | None = None,
        # Configure a custom httpx client.
        # We provide a `DefaultAsyncHttpxClient` class that you can pass to retain the default values we use for `limits`, `timeout` & `follow_redirects`.
        # See the [httpx documentation](https://www.python-httpx.org/api/#asyncclient) for more details.
        http_client: httpx.AsyncClient | None = None,
        # Enable or disable schema validation for data returned by the API.
        # When enabled an error APIResponseValidationError is raised
        # if the API responds with invalid data for the expected schema.
        #
        # This parameter may be removed or changed in the future.
        # If you rely on this feature, please open a GitHub issue
        # outlining your use-case to help us decide if it should be
        # part of our public interface in the future.
        _strict_response_validation: bool = False,
    ) -> None:
        """Construct a new async deeporigin_data client instance."""
        self.bearer_token = bearer_token

        self.org_id = org_id

        self._environment = environment

        base_url_env = os.environ.get("DEEPORIGIN_DATA_BASE_URL")
        if is_given(base_url) and base_url is not None:
            # cast required because mypy doesn't understand the type narrowing
            base_url = cast("str | httpx.URL", base_url)  # pyright: ignore[reportUnnecessaryCast]
        elif is_given(environment):
            if base_url_env and base_url is not None:
                raise ValueError(
                    "Ambiguous URL; The `DEEPORIGIN_DATA_BASE_URL` env var and the `environment` argument are given. If you want to use the environment, you must pass base_url=None",
                )

            try:
                base_url = ENVIRONMENTS[environment]
            except KeyError as exc:
                raise ValueError(f"Unknown environment: {environment}") from exc
        elif base_url_env is not None:
            base_url = base_url_env
        else:
            self._environment = environment = "edge"

            try:
                base_url = ENVIRONMENTS[environment]
            except KeyError as exc:
                raise ValueError(f"Unknown environment: {environment}") from exc

        super().__init__(
            version=__version__,
            base_url=base_url,
            max_retries=max_retries,
            timeout=timeout,
            http_client=http_client,
            custom_headers=default_headers,
            custom_query=default_query,
            _strict_response_validation=_strict_response_validation,
        )

        self.add_database_column = resources.AsyncAddDatabaseColumnResource(self)
        self.archive_files = resources.AsyncArchiveFilesResource(self)
        self.configure_column_select_options = resources.AsyncConfigureColumnSelectOptionsResource(self)
        self.convert_id_format = resources.AsyncConvertIDFormatResource(self)
        self.databases = resources.AsyncDatabasesResource(self)
        self.create_file_download_url = resources.AsyncCreateFileDownloadURLResource(self)
        self.create_file_upload = resources.AsyncCreateFileUploadResource(self)
        self.workspaces = resources.AsyncWorkspacesResource(self)
        self.delete_database_column = resources.AsyncDeleteDatabaseColumnResource(self)
        self.delete_rows = resources.AsyncDeleteRowsResource(self)
        self.describe_code_execution = resources.AsyncDescribeCodeExecutionResource(self)
        self.describe_database_stats = resources.AsyncDescribeDatabaseStatsResource(self)
        self.describe_file = resources.AsyncDescribeFileResource(self)
        self.describe_row = resources.AsyncDescribeRowResource(self)
        self.download_file = resources.AsyncDownloadFileResource(self)
        self.ensure_rows = resources.AsyncEnsureRowsResource(self)
        self.execute_code = resources.AsyncExecuteCodeResource(self)
        self.execute_code_sync = resources.AsyncExecuteCodeSyncResource(self)
        self.get_code_execution_result = resources.AsyncGetCodeExecutionResultResource(self)
        self.import_rows = resources.AsyncImportRowsResource(self)
        self.organization = resources.AsyncOrganizationResource(self)
        self.chat_threads = resources.AsyncChatThreadsResource(self)
        self.files = resources.AsyncFilesResource(self)
        self.mentions = resources.AsyncMentionsResource(self)
        self.rows = resources.AsyncRowsResource(self)
        self.sequences = resources.AsyncSequencesResource(self)
        self.chat_messages = resources.AsyncChatMessagesResource(self)
        self.update_database = resources.AsyncUpdateDatabaseResource(self)
        self.update_database_column = resources.AsyncUpdateDatabaseColumnResource(self)
        self.update_workspace = resources.AsyncUpdateWorkspaceResource(self)
        self.with_raw_response = AsyncDeeporiginDataWithRawResponse(self)
        self.with_streaming_response = AsyncDeeporiginDataWithStreamedResponse(self)

    @property
    @override
    def qs(self) -> Querystring:
        return Querystring(array_format="comma")

    @property
    @override
    def auth_headers(self) -> dict[str, str]:
        bearer_token = self.bearer_token
        return {"Authorization": f"Bearer {bearer_token}"}

    @property
    @override
    def default_headers(self) -> dict[str, str | Omit]:
        return {
            **super().default_headers,
            "X-Stainless-Async": f"async:{get_async_library()}",
            "x-org-id": self.org_id,
            **self._custom_headers,
        }

    def copy(
        self,
        *,
        bearer_token: str | None = None,
        org_id: str | None = None,
        environment: Literal["edge", "prod"] | None = None,
        base_url: str | httpx.URL | None = None,
        timeout: float | Timeout | None | NotGiven = NOT_GIVEN,
        http_client: httpx.AsyncClient | None = None,
        max_retries: int | NotGiven = NOT_GIVEN,
        default_headers: Mapping[str, str] | None = None,
        set_default_headers: Mapping[str, str] | None = None,
        default_query: Mapping[str, object] | None = None,
        set_default_query: Mapping[str, object] | None = None,
        _extra_kwargs: Mapping[str, Any] = {},
    ) -> Self:
        """
        Create a new client instance re-using the same options given to the current client with optional overriding.
        """
        if default_headers is not None and set_default_headers is not None:
            raise ValueError("The `default_headers` and `set_default_headers` arguments are mutually exclusive")

        if default_query is not None and set_default_query is not None:
            raise ValueError("The `default_query` and `set_default_query` arguments are mutually exclusive")

        headers = self._custom_headers
        if default_headers is not None:
            headers = {**headers, **default_headers}
        elif set_default_headers is not None:
            headers = set_default_headers

        params = self._custom_query
        if default_query is not None:
            params = {**params, **default_query}
        elif set_default_query is not None:
            params = set_default_query

        http_client = http_client or self._client
        return self.__class__(
            bearer_token=bearer_token or self.bearer_token,
            org_id=org_id or self.org_id,
            base_url=base_url or self.base_url,
            environment=environment or self._environment,
            timeout=self.timeout if isinstance(timeout, NotGiven) else timeout,
            http_client=http_client,
            max_retries=max_retries if is_given(max_retries) else self.max_retries,
            default_headers=headers,
            default_query=params,
            **_extra_kwargs,
        )

    # Alias for `copy` for nicer inline usage, e.g.
    # client.with_options(timeout=10).foo.create(...)
    with_options = copy

    @override
    def _make_status_error(
        self,
        err_msg: str,
        *,
        body: object,
        response: httpx.Response,
    ) -> APIStatusError:
        if response.status_code == 400:
            return _exceptions.BadRequestError(err_msg, response=response, body=body)

        if response.status_code == 401:
            return _exceptions.AuthenticationError(err_msg, response=response, body=body)

        if response.status_code == 403:
            return _exceptions.PermissionDeniedError(err_msg, response=response, body=body)

        if response.status_code == 404:
            return _exceptions.NotFoundError(err_msg, response=response, body=body)

        if response.status_code == 409:
            return _exceptions.ConflictError(err_msg, response=response, body=body)

        if response.status_code == 422:
            return _exceptions.UnprocessableEntityError(err_msg, response=response, body=body)

        if response.status_code == 429:
            return _exceptions.RateLimitError(err_msg, response=response, body=body)

        if response.status_code >= 500:
            return _exceptions.InternalServerError(err_msg, response=response, body=body)
        return APIStatusError(err_msg, response=response, body=body)


class DeeporiginDataWithRawResponse:
    def __init__(self, client: DeeporiginData) -> None:
        self.add_database_column = resources.AddDatabaseColumnResourceWithRawResponse(client.add_database_column)
        self.archive_files = resources.ArchiveFilesResourceWithRawResponse(client.archive_files)
        self.configure_column_select_options = resources.ConfigureColumnSelectOptionsResourceWithRawResponse(
            client.configure_column_select_options
        )
        self.convert_id_format = resources.ConvertIDFormatResourceWithRawResponse(client.convert_id_format)
        self.databases = resources.DatabasesResourceWithRawResponse(client.databases)
        self.create_file_download_url = resources.CreateFileDownloadURLResourceWithRawResponse(
            client.create_file_download_url
        )
        self.create_file_upload = resources.CreateFileUploadResourceWithRawResponse(client.create_file_upload)
        self.workspaces = resources.WorkspacesResourceWithRawResponse(client.workspaces)
        self.delete_database_column = resources.DeleteDatabaseColumnResourceWithRawResponse(
            client.delete_database_column
        )
        self.delete_rows = resources.DeleteRowsResourceWithRawResponse(client.delete_rows)
        self.describe_code_execution = resources.DescribeCodeExecutionResourceWithRawResponse(
            client.describe_code_execution
        )
        self.describe_database_stats = resources.DescribeDatabaseStatsResourceWithRawResponse(
            client.describe_database_stats
        )
        self.describe_file = resources.DescribeFileResourceWithRawResponse(client.describe_file)
        self.describe_row = resources.DescribeRowResourceWithRawResponse(client.describe_row)
        self.download_file = resources.DownloadFileResourceWithRawResponse(client.download_file)
        self.ensure_rows = resources.EnsureRowsResourceWithRawResponse(client.ensure_rows)
        self.execute_code = resources.ExecuteCodeResourceWithRawResponse(client.execute_code)
        self.execute_code_sync = resources.ExecuteCodeSyncResourceWithRawResponse(client.execute_code_sync)
        self.get_code_execution_result = resources.GetCodeExecutionResultResourceWithRawResponse(
            client.get_code_execution_result
        )
        self.import_rows = resources.ImportRowsResourceWithRawResponse(client.import_rows)
        self.organization = resources.OrganizationResourceWithRawResponse(client.organization)
        self.chat_threads = resources.ChatThreadsResourceWithRawResponse(client.chat_threads)
        self.files = resources.FilesResourceWithRawResponse(client.files)
        self.mentions = resources.MentionsResourceWithRawResponse(client.mentions)
        self.rows = resources.RowsResourceWithRawResponse(client.rows)
        self.sequences = resources.SequencesResourceWithRawResponse(client.sequences)
        self.chat_messages = resources.ChatMessagesResourceWithRawResponse(client.chat_messages)
        self.update_database = resources.UpdateDatabaseResourceWithRawResponse(client.update_database)
        self.update_database_column = resources.UpdateDatabaseColumnResourceWithRawResponse(
            client.update_database_column
        )
        self.update_workspace = resources.UpdateWorkspaceResourceWithRawResponse(client.update_workspace)


class AsyncDeeporiginDataWithRawResponse:
    def __init__(self, client: AsyncDeeporiginData) -> None:
        self.add_database_column = resources.AsyncAddDatabaseColumnResourceWithRawResponse(client.add_database_column)
        self.archive_files = resources.AsyncArchiveFilesResourceWithRawResponse(client.archive_files)
        self.configure_column_select_options = resources.AsyncConfigureColumnSelectOptionsResourceWithRawResponse(
            client.configure_column_select_options
        )
        self.convert_id_format = resources.AsyncConvertIDFormatResourceWithRawResponse(client.convert_id_format)
        self.databases = resources.AsyncDatabasesResourceWithRawResponse(client.databases)
        self.create_file_download_url = resources.AsyncCreateFileDownloadURLResourceWithRawResponse(
            client.create_file_download_url
        )
        self.create_file_upload = resources.AsyncCreateFileUploadResourceWithRawResponse(client.create_file_upload)
        self.workspaces = resources.AsyncWorkspacesResourceWithRawResponse(client.workspaces)
        self.delete_database_column = resources.AsyncDeleteDatabaseColumnResourceWithRawResponse(
            client.delete_database_column
        )
        self.delete_rows = resources.AsyncDeleteRowsResourceWithRawResponse(client.delete_rows)
        self.describe_code_execution = resources.AsyncDescribeCodeExecutionResourceWithRawResponse(
            client.describe_code_execution
        )
        self.describe_database_stats = resources.AsyncDescribeDatabaseStatsResourceWithRawResponse(
            client.describe_database_stats
        )
        self.describe_file = resources.AsyncDescribeFileResourceWithRawResponse(client.describe_file)
        self.describe_row = resources.AsyncDescribeRowResourceWithRawResponse(client.describe_row)
        self.download_file = resources.AsyncDownloadFileResourceWithRawResponse(client.download_file)
        self.ensure_rows = resources.AsyncEnsureRowsResourceWithRawResponse(client.ensure_rows)
        self.execute_code = resources.AsyncExecuteCodeResourceWithRawResponse(client.execute_code)
        self.execute_code_sync = resources.AsyncExecuteCodeSyncResourceWithRawResponse(client.execute_code_sync)
        self.get_code_execution_result = resources.AsyncGetCodeExecutionResultResourceWithRawResponse(
            client.get_code_execution_result
        )
        self.import_rows = resources.AsyncImportRowsResourceWithRawResponse(client.import_rows)
        self.organization = resources.AsyncOrganizationResourceWithRawResponse(client.organization)
        self.chat_threads = resources.AsyncChatThreadsResourceWithRawResponse(client.chat_threads)
        self.files = resources.AsyncFilesResourceWithRawResponse(client.files)
        self.mentions = resources.AsyncMentionsResourceWithRawResponse(client.mentions)
        self.rows = resources.AsyncRowsResourceWithRawResponse(client.rows)
        self.sequences = resources.AsyncSequencesResourceWithRawResponse(client.sequences)
        self.chat_messages = resources.AsyncChatMessagesResourceWithRawResponse(client.chat_messages)
        self.update_database = resources.AsyncUpdateDatabaseResourceWithRawResponse(client.update_database)
        self.update_database_column = resources.AsyncUpdateDatabaseColumnResourceWithRawResponse(
            client.update_database_column
        )
        self.update_workspace = resources.AsyncUpdateWorkspaceResourceWithRawResponse(client.update_workspace)


class DeeporiginDataWithStreamedResponse:
    def __init__(self, client: DeeporiginData) -> None:
        self.add_database_column = resources.AddDatabaseColumnResourceWithStreamingResponse(client.add_database_column)
        self.archive_files = resources.ArchiveFilesResourceWithStreamingResponse(client.archive_files)
        self.configure_column_select_options = resources.ConfigureColumnSelectOptionsResourceWithStreamingResponse(
            client.configure_column_select_options
        )
        self.convert_id_format = resources.ConvertIDFormatResourceWithStreamingResponse(client.convert_id_format)
        self.databases = resources.DatabasesResourceWithStreamingResponse(client.databases)
        self.create_file_download_url = resources.CreateFileDownloadURLResourceWithStreamingResponse(
            client.create_file_download_url
        )
        self.create_file_upload = resources.CreateFileUploadResourceWithStreamingResponse(client.create_file_upload)
        self.workspaces = resources.WorkspacesResourceWithStreamingResponse(client.workspaces)
        self.delete_database_column = resources.DeleteDatabaseColumnResourceWithStreamingResponse(
            client.delete_database_column
        )
        self.delete_rows = resources.DeleteRowsResourceWithStreamingResponse(client.delete_rows)
        self.describe_code_execution = resources.DescribeCodeExecutionResourceWithStreamingResponse(
            client.describe_code_execution
        )
        self.describe_database_stats = resources.DescribeDatabaseStatsResourceWithStreamingResponse(
            client.describe_database_stats
        )
        self.describe_file = resources.DescribeFileResourceWithStreamingResponse(client.describe_file)
        self.describe_row = resources.DescribeRowResourceWithStreamingResponse(client.describe_row)
        self.download_file = resources.DownloadFileResourceWithStreamingResponse(client.download_file)
        self.ensure_rows = resources.EnsureRowsResourceWithStreamingResponse(client.ensure_rows)
        self.execute_code = resources.ExecuteCodeResourceWithStreamingResponse(client.execute_code)
        self.execute_code_sync = resources.ExecuteCodeSyncResourceWithStreamingResponse(client.execute_code_sync)
        self.get_code_execution_result = resources.GetCodeExecutionResultResourceWithStreamingResponse(
            client.get_code_execution_result
        )
        self.import_rows = resources.ImportRowsResourceWithStreamingResponse(client.import_rows)
        self.organization = resources.OrganizationResourceWithStreamingResponse(client.organization)
        self.chat_threads = resources.ChatThreadsResourceWithStreamingResponse(client.chat_threads)
        self.files = resources.FilesResourceWithStreamingResponse(client.files)
        self.mentions = resources.MentionsResourceWithStreamingResponse(client.mentions)
        self.rows = resources.RowsResourceWithStreamingResponse(client.rows)
        self.sequences = resources.SequencesResourceWithStreamingResponse(client.sequences)
        self.chat_messages = resources.ChatMessagesResourceWithStreamingResponse(client.chat_messages)
        self.update_database = resources.UpdateDatabaseResourceWithStreamingResponse(client.update_database)
        self.update_database_column = resources.UpdateDatabaseColumnResourceWithStreamingResponse(
            client.update_database_column
        )
        self.update_workspace = resources.UpdateWorkspaceResourceWithStreamingResponse(client.update_workspace)


class AsyncDeeporiginDataWithStreamedResponse:
    def __init__(self, client: AsyncDeeporiginData) -> None:
        self.add_database_column = resources.AsyncAddDatabaseColumnResourceWithStreamingResponse(
            client.add_database_column
        )
        self.archive_files = resources.AsyncArchiveFilesResourceWithStreamingResponse(client.archive_files)
        self.configure_column_select_options = resources.AsyncConfigureColumnSelectOptionsResourceWithStreamingResponse(
            client.configure_column_select_options
        )
        self.convert_id_format = resources.AsyncConvertIDFormatResourceWithStreamingResponse(client.convert_id_format)
        self.databases = resources.AsyncDatabasesResourceWithStreamingResponse(client.databases)
        self.create_file_download_url = resources.AsyncCreateFileDownloadURLResourceWithStreamingResponse(
            client.create_file_download_url
        )
        self.create_file_upload = resources.AsyncCreateFileUploadResourceWithStreamingResponse(
            client.create_file_upload
        )
        self.workspaces = resources.AsyncWorkspacesResourceWithStreamingResponse(client.workspaces)
        self.delete_database_column = resources.AsyncDeleteDatabaseColumnResourceWithStreamingResponse(
            client.delete_database_column
        )
        self.delete_rows = resources.AsyncDeleteRowsResourceWithStreamingResponse(client.delete_rows)
        self.describe_code_execution = resources.AsyncDescribeCodeExecutionResourceWithStreamingResponse(
            client.describe_code_execution
        )
        self.describe_database_stats = resources.AsyncDescribeDatabaseStatsResourceWithStreamingResponse(
            client.describe_database_stats
        )
        self.describe_file = resources.AsyncDescribeFileResourceWithStreamingResponse(client.describe_file)
        self.describe_row = resources.AsyncDescribeRowResourceWithStreamingResponse(client.describe_row)
        self.download_file = resources.AsyncDownloadFileResourceWithStreamingResponse(client.download_file)
        self.ensure_rows = resources.AsyncEnsureRowsResourceWithStreamingResponse(client.ensure_rows)
        self.execute_code = resources.AsyncExecuteCodeResourceWithStreamingResponse(client.execute_code)
        self.execute_code_sync = resources.AsyncExecuteCodeSyncResourceWithStreamingResponse(client.execute_code_sync)
        self.get_code_execution_result = resources.AsyncGetCodeExecutionResultResourceWithStreamingResponse(
            client.get_code_execution_result
        )
        self.import_rows = resources.AsyncImportRowsResourceWithStreamingResponse(client.import_rows)
        self.organization = resources.AsyncOrganizationResourceWithStreamingResponse(client.organization)
        self.chat_threads = resources.AsyncChatThreadsResourceWithStreamingResponse(client.chat_threads)
        self.files = resources.AsyncFilesResourceWithStreamingResponse(client.files)
        self.mentions = resources.AsyncMentionsResourceWithStreamingResponse(client.mentions)
        self.rows = resources.AsyncRowsResourceWithStreamingResponse(client.rows)
        self.sequences = resources.AsyncSequencesResourceWithStreamingResponse(client.sequences)
        self.chat_messages = resources.AsyncChatMessagesResourceWithStreamingResponse(client.chat_messages)
        self.update_database = resources.AsyncUpdateDatabaseResourceWithStreamingResponse(client.update_database)
        self.update_database_column = resources.AsyncUpdateDatabaseColumnResourceWithStreamingResponse(
            client.update_database_column
        )
        self.update_workspace = resources.AsyncUpdateWorkspaceResourceWithStreamingResponse(client.update_workspace)


Client = DeeporiginData

AsyncClient = AsyncDeeporiginData
