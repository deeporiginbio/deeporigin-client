from .secret_file_value import SecretFileValue


class MosekLicenseFile(SecretFileValue):
    """Mosek license file"""

    class Meta:
        platform_id = "Mosek license file"

    @classmethod
    @property
    def FILENAME(cls) -> str:
        return "~/mosek/mosek.lic"
