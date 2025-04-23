# Generated API client - DO NOT EDIT MANUALLY
# This file was auto-generated from an OpenAPI schema

import httpx
from .schema_api_httpx_client import FilePayload


def _get_get_signed_url_download_kwargs(
    orgFriendlyId: str,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for getSignedUrlDownload"""
    _kwargs: dict[str, any] = {
        "method": "get",
        "url": f"/file-api/signedUrl/download/{orgFriendlyId}/*".format(orgFriendlyId=orgFriendlyId),
    }
    return _kwargs, "/file-api/signedUrl/download/{orgFriendlyId}/*"

def _get_get_signed_url_upload_kwargs(
    orgFriendlyId: str,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for getSignedUrlUpload"""
    _kwargs: dict[str, any] = {
        "method": "get",
        "url": f"/file-api/signedUrl/upload/{orgFriendlyId}/*".format(orgFriendlyId=orgFriendlyId),
    }
    return _kwargs, "/file-api/signedUrl/upload/{orgFriendlyId}/*"

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
    list_type: float | None = None,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for getObject"""
    _kwargs: dict[str, any] = {
        "method": "get",
        "url": f"/file-api/{orgFriendlyId}/{filePath}".format(orgFriendlyId=orgFriendlyId, filePath=filePath),
    }
    params: dict[str, any] = {}
    if list_type is not None:
        params["list-type"] = list_type
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

def _get_sync_objects_kwargs(
    orgFriendlyId: str,
    filePath: str,
    body: dict[str, any] | None = None,
) -> tuple[dict[str, any], str]:
    """Get request kwargs for syncObjects"""
    _kwargs: dict[str, any] = {
        "method": "post",
        "url": f"/file-api/{orgFriendlyId}/{filePath}".format(orgFriendlyId=orgFriendlyId, filePath=filePath),
    }
    if body is not None:
        _kwargs["json"] = body
        _kwargs["headers"] = {"Content-Type": "application/json"}
    return _kwargs, "/file-api/{orgFriendlyId}/{filePath}"

def _get_check_kwargs(
) -> tuple[dict[str, any], str]:
    """Get request kwargs for check"""
    _kwargs: dict[str, any] = {
        "method": "get",
        "url": "/file-api/health",
    }
    return _kwargs, "/file-api/health"
