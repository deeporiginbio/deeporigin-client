from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, TypeVar, Union, cast

from attrs import define as _attrs_define
from attrs import field as _attrs_field

from ..types import UNSET, Unset

if TYPE_CHECKING:
    from ..models.check_response_503_details import CheckResponse503Details
    from ..models.check_response_503_error_type_0 import CheckResponse503ErrorType0
    from ..models.check_response_503_info_type_0 import CheckResponse503InfoType0


T = TypeVar("T", bound="CheckResponse503")


@_attrs_define
class CheckResponse503:
    """
    Attributes:
        status (Union[Unset, str]):  Example: error.
        info (Union['CheckResponse503InfoType0', None, Unset]):  Example: {'database': {'status': 'up'}}.
        error (Union['CheckResponse503ErrorType0', None, Unset]):  Example: {'redis': {'status': 'down', 'message':
            'Could not connect'}}.
        details (Union[Unset, CheckResponse503Details]):  Example: {'database': {'status': 'up'}, 'redis': {'status':
            'down', 'message': 'Could not connect'}}.
    """

    status: Union[Unset, str] = UNSET
    info: Union["CheckResponse503InfoType0", None, Unset] = UNSET
    error: Union["CheckResponse503ErrorType0", None, Unset] = UNSET
    details: Union[Unset, "CheckResponse503Details"] = UNSET
    additional_properties: dict[str, Any] = _attrs_field(init=False, factory=dict)

    def to_dict(self) -> dict[str, Any]:
        from ..models.check_response_503_error_type_0 import CheckResponse503ErrorType0
        from ..models.check_response_503_info_type_0 import CheckResponse503InfoType0

        status = self.status

        info: Union[None, Unset, dict[str, Any]]
        if isinstance(self.info, Unset):
            info = UNSET
        elif isinstance(self.info, CheckResponse503InfoType0):
            info = self.info.to_dict()
        else:
            info = self.info

        error: Union[None, Unset, dict[str, Any]]
        if isinstance(self.error, Unset):
            error = UNSET
        elif isinstance(self.error, CheckResponse503ErrorType0):
            error = self.error.to_dict()
        else:
            error = self.error

        details: Union[Unset, dict[str, Any]] = UNSET
        if not isinstance(self.details, Unset):
            details = self.details.to_dict()

        field_dict: dict[str, Any] = {}
        field_dict.update(self.additional_properties)
        field_dict.update({})
        if status is not UNSET:
            field_dict["status"] = status
        if info is not UNSET:
            field_dict["info"] = info
        if error is not UNSET:
            field_dict["error"] = error
        if details is not UNSET:
            field_dict["details"] = details

        return field_dict

    @classmethod
    def from_dict(cls: type[T], src_dict: Mapping[str, Any]) -> T:
        from ..models.check_response_503_details import CheckResponse503Details
        from ..models.check_response_503_error_type_0 import CheckResponse503ErrorType0
        from ..models.check_response_503_info_type_0 import CheckResponse503InfoType0

        d = dict(src_dict)
        status = d.pop("status", UNSET)

        def _parse_info(data: object) -> Union["CheckResponse503InfoType0", None, Unset]:
            if data is None:
                return data
            if isinstance(data, Unset):
                return data
            try:
                if not isinstance(data, dict):
                    raise TypeError()
                info_type_0 = CheckResponse503InfoType0.from_dict(data)

                return info_type_0
            except:  # noqa: E722
                pass
            return cast(Union["CheckResponse503InfoType0", None, Unset], data)

        info = _parse_info(d.pop("info", UNSET))

        def _parse_error(data: object) -> Union["CheckResponse503ErrorType0", None, Unset]:
            if data is None:
                return data
            if isinstance(data, Unset):
                return data
            try:
                if not isinstance(data, dict):
                    raise TypeError()
                error_type_0 = CheckResponse503ErrorType0.from_dict(data)

                return error_type_0
            except:  # noqa: E722
                pass
            return cast(Union["CheckResponse503ErrorType0", None, Unset], data)

        error = _parse_error(d.pop("error", UNSET))

        _details = d.pop("details", UNSET)
        details: Union[Unset, CheckResponse503Details]
        if isinstance(_details, Unset):
            details = UNSET
        else:
            details = CheckResponse503Details.from_dict(_details)

        check_response_503 = cls(
            status=status,
            info=info,
            error=error,
            details=details,
        )

        check_response_503.additional_properties = d
        return check_response_503

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
