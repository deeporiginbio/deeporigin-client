from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, TypeVar

from attrs import define as _attrs_define
from attrs import field as _attrs_field

from ..models.sync_file_schema_dto_provider import SyncFileSchemaDtoProvider

if TYPE_CHECKING:
    from ..models.sync_file_schema_dto_credentials import SyncFileSchemaDtoCredentials


T = TypeVar("T", bound="SyncFileSchemaDto")


@_attrs_define
class SyncFileSchemaDto:
    """
    Attributes:
        provider (SyncFileSchemaDtoProvider):
        credentials (SyncFileSchemaDtoCredentials):
        path (str):
    """

    provider: SyncFileSchemaDtoProvider
    credentials: "SyncFileSchemaDtoCredentials"
    path: str
    additional_properties: dict[str, Any] = _attrs_field(init=False, factory=dict)

    def to_dict(self) -> dict[str, Any]:
        provider = self.provider.value

        credentials = self.credentials.to_dict()

        path = self.path

        field_dict: dict[str, Any] = {}
        field_dict.update(self.additional_properties)
        field_dict.update(
            {
                "provider": provider,
                "credentials": credentials,
                "path": path,
            }
        )

        return field_dict

    @classmethod
    def from_dict(cls: type[T], src_dict: Mapping[str, Any]) -> T:
        from ..models.sync_file_schema_dto_credentials import SyncFileSchemaDtoCredentials

        d = dict(src_dict)
        provider = SyncFileSchemaDtoProvider(d.pop("provider"))

        credentials = SyncFileSchemaDtoCredentials.from_dict(d.pop("credentials"))

        path = d.pop("path")

        sync_file_schema_dto = cls(
            provider=provider,
            credentials=credentials,
            path=path,
        )

        sync_file_schema_dto.additional_properties = d
        return sync_file_schema_dto

    @property
    def additional_keys(self) -> list[str]:
        return list(self.additional_properties.keys())

    def __getitem__(self, key: str) -> Any:
        return self.additional_properties[key]

    def __setitem__(self, key: str, value: Any) -> None:
        self.additional_properties[key] = value

    def __delitem__(self, key: str) -> None:
        del self.additional_properties[key]

    def __contains__(self, key: str) -> bool:
        return key in self.additional_properties
