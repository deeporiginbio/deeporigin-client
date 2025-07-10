# Generated API client - DO NOT EDIT MANUALLY
# This file was auto-generated from an OpenAPI schema

from .schema_api_httpx_client import FilePayload


def _get_head_shared_object_kwargs(
    filePath: str,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for headSharedObject"""
    _kwargs: dict[str, any] = {
        "method": "head",
        "url": f"/file-api/shared-data/{filePath}".format(filePath=filePath),
    }
    return _kwargs, "/file-api/shared-data/{filePath}"

def _get_get_do_object_kwargs(
    filePath: str,
    last_count: str | None = None,
    continuation_token: str | None = None,
    delimiter: str | None = None,
    recursive: str | None = None,
    max_keys: float | None = None,
    list_type: float | None = None,
    prefix: str | None = None,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for getDoObject"""
    _kwargs: dict[str, any] = {
        "method": "get",
        "url": f"/file-api/shared-data/{filePath}".format(filePath=filePath),
    }
    params: dict[str, any] = {}
    if last_count is not None:
        params["last-count"] = last_count
    if continuation_token is not None:
        params["continuation-token"] = continuation_token
    if delimiter is not None:
        params["delimiter"] = delimiter
    if recursive is not None:
        params["recursive"] = recursive
    if max_keys is not None:
        params["max-keys"] = max_keys
    if list_type is not None:
        params["list-type"] = list_type
    if prefix is not None:
        params["prefix"] = prefix
    if params:
        _kwargs["params"] = params
    return _kwargs, "/file-api/shared-data/{filePath}"

def _get_delete_shared_object_kwargs(
    filePath: str,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for deleteSharedObject"""
    _kwargs: dict[str, any] = {
        "method": "delete",
        "url": f"/file-api/shared-data/{filePath}".format(filePath=filePath),
    }
    return _kwargs, "/file-api/shared-data/{filePath}"

def _get_put_shared_object_kwargs(
    filePath: str,
    file: FilePayload | None = None,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for putSharedObject"""
    _kwargs: dict[str, any] = {
        "method": "put",
        "url": f"/file-api/shared-data/{filePath}".format(filePath=filePath),
    }
    if file is not None:
        _kwargs["files"] = {"file": file.to_tuple()[1]}
    return _kwargs, "/file-api/shared-data/{filePath}"

def _get_head_do_object_kwargs(
    filePath: str,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for headDoObject"""
    _kwargs: dict[str, any] = {
        "method": "head",
        "url": f"/file-api/do-data/{filePath}".format(filePath=filePath),
    }
    return _kwargs, "/file-api/do-data/{filePath}"

def _get_get_do_object_file_api_kwargs(
    filePath: str,
    last_count: str | None = None,
    continuation_token: str | None = None,
    delimiter: str | None = None,
    recursive: str | None = None,
    max_keys: float | None = None,
    list_type: float | None = None,
    prefix: str | None = None,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for getDoObject"""
    _kwargs: dict[str, any] = {
        "method": "get",
        "url": f"/file-api/do-data/{filePath}".format(filePath=filePath),
    }
    params: dict[str, any] = {}
    if last_count is not None:
        params["last-count"] = last_count
    if continuation_token is not None:
        params["continuation-token"] = continuation_token
    if delimiter is not None:
        params["delimiter"] = delimiter
    if recursive is not None:
        params["recursive"] = recursive
    if max_keys is not None:
        params["max-keys"] = max_keys
    if list_type is not None:
        params["list-type"] = list_type
    if prefix is not None:
        params["prefix"] = prefix
    if params:
        _kwargs["params"] = params
    return _kwargs, "/file-api/do-data/{filePath}"

def _get_delete_do_object_kwargs(
    filePath: str,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for deleteDoObject"""
    _kwargs: dict[str, any] = {
        "method": "delete",
        "url": f"/file-api/do-data/{filePath}".format(filePath=filePath),
    }
    return _kwargs, "/file-api/do-data/{filePath}"

def _get_put_do_object_kwargs(
    filePath: str,
    file: FilePayload | None = None,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for putDoObject"""
    _kwargs: dict[str, any] = {
        "method": "put",
        "url": f"/file-api/do-data/{filePath}".format(filePath=filePath),
    }
    if file is not None:
        _kwargs["files"] = {"file": file.to_tuple()[1]}
    return _kwargs, "/file-api/do-data/{filePath}"

def _get_get_signed_url_download_kwargs(
    orgFriendlyId: str,
    filePath: str,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for getSignedUrlDownload"""
    _kwargs: dict[str, any] = {
        "method": "get",
        "url": f"/file-api/signedUrl/download/{orgFriendlyId}/{filePath}".format(orgFriendlyId=orgFriendlyId, filePath=filePath),
    }
    return _kwargs, "/file-api/signedUrl/download/{orgFriendlyId}/{filePath}"

def _get_get_signed_url_upload_kwargs(
    orgFriendlyId: str,
    filePath: str,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for getSignedUrlUpload"""
    _kwargs: dict[str, any] = {
        "method": "get",
        "url": f"/file-api/signedUrl/upload/{orgFriendlyId}/{filePath}".format(orgFriendlyId=orgFriendlyId, filePath=filePath),
    }
    return _kwargs, "/file-api/signedUrl/upload/{orgFriendlyId}/{filePath}"

def _get_put_object_kwargs(
    orgFriendlyId: str,
    filePath: str,
    file: FilePayload | None = None,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for putObject"""
    _kwargs: dict[str, any] = {
        "method": "put",
        "url": f"/file-api/{orgFriendlyId}/{filePath}".format(orgFriendlyId=orgFriendlyId, filePath=filePath),
    }
    if file is not None:
        _kwargs["files"] = {"file": file.to_tuple()[1]}
    return _kwargs, "/file-api/{orgFriendlyId}/{filePath}"

def _get_head_object_kwargs(
    orgFriendlyId: str,
    filePath: str,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for headObject"""
    _kwargs: dict[str, any] = {
        "method": "head",
        "url": f"/file-api/{orgFriendlyId}/{filePath}".format(orgFriendlyId=orgFriendlyId, filePath=filePath),
    }
    return _kwargs, "/file-api/{orgFriendlyId}/{filePath}"

def _get_get_object_kwargs(
    orgFriendlyId: str,
    filePath: str,
    last_count: str | None = None,
    continuation_token: str | None = None,
    delimiter: str | None = None,
    recursive: str | None = None,
    max_keys: float | None = None,
    list_type: float | None = None,
    prefix: str | None = None,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for getObject"""
    _kwargs: dict[str, any] = {
        "method": "get",
        "url": f"/file-api/{orgFriendlyId}/{filePath}".format(orgFriendlyId=orgFriendlyId, filePath=filePath),
    }
    params: dict[str, any] = {}
    if last_count is not None:
        params["last-count"] = last_count
    if continuation_token is not None:
        params["continuation-token"] = continuation_token
    if delimiter is not None:
        params["delimiter"] = delimiter
    if recursive is not None:
        params["recursive"] = recursive
    if max_keys is not None:
        params["max-keys"] = max_keys
    if list_type is not None:
        params["list-type"] = list_type
    if prefix is not None:
        params["prefix"] = prefix
    if params:
        _kwargs["params"] = params
    return _kwargs, "/file-api/{orgFriendlyId}/{filePath}"

def _get_delete_object_kwargs(
    orgFriendlyId: str,
    filePath: str,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for deleteObject"""
    _kwargs: dict[str, any] = {
        "method": "delete",
        "url": f"/file-api/{orgFriendlyId}/{filePath}".format(orgFriendlyId=orgFriendlyId, filePath=filePath),
    }
    return _kwargs, "/file-api/{orgFriendlyId}/{filePath}"

# Manually added while Niels adds health endpoint
def _get_check_kwargs(
) -> tuple[dict[str, any], str]:
    """Get request kwargs for check"""
    _kwargs: dict[str, any] = {
        "method": "get",
        "url": "/file-api/health",
    }
    return _kwargs, "/file-api/health"
