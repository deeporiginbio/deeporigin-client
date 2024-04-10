from .file import File


class SecretFile(File):
    """Secret file"""

    class Meta:
        platform_id = "Secret file"
