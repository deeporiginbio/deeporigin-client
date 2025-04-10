from collections.abc import Mapping
from typing import Any, TypeVar

from attrs import define as _attrs_define
from attrs import field as _attrs_field

T = TypeVar("T", bound="SyncFileSchemaDtoCredentials")


@_attrs_define
class SyncFileSchemaDtoCredentials:
    """
    Attributes:
        region (str):
        secret_key (str):
        access_key (str):
    """

    region: str
    secret_key: str
    access_key: str
    additional_properties: dict[str, Any] = _attrs_field(init=False, factory=dict)

    def to_dict(self) -> dict[str, Any]:
        region = self.region

        secret_key = self.secret_key

        access_key = self.access_key

        field_dict: dict[str, Any] = {}
        field_dict.update(self.additional_properties)
        field_dict.update(
            {
                "region": region,
                "secretKey": secret_key,
                "accessKey": access_key,
            }
        )

        return field_dict

    @classmethod
    def from_dict(cls: type[T], src_dict: Mapping[str, Any]) -> T:
        d = dict(src_dict)
        region = d.pop("region")

        secret_key = d.pop("secretKey")

        access_key = d.pop("accessKey")

        sync_file_schema_dto_credentials = cls(
            region=region,
            secret_key=secret_key,
            access_key=access_key,
        )

        sync_file_schema_dto_credentials.additional_properties = d
        return sync_file_schema_dto_credentials

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
