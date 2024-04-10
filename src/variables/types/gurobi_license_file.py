from .secret_file_value import SecretFileValue


class GurobiLicenseFile(SecretFileValue):
    """Gurobi license file"""

    class Meta:
        platform_id = "Gurobi license file"

    @classmethod
    @property
    def FILENAME(cls) -> str:
        return "~/gurobi.lic"
