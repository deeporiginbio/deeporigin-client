import os
import pathlib
import subprocess
from pathlib import Path

__all__ = [
    "__version__",
]

SRC_DIR = pathlib.Path(__file__).parent

version_filename = os.path.join(SRC_DIR, "VERSION")

with open(version_filename, "r") as file:
    __version__ = file.read()

    if __version__ == "0.0.0.dev0":
        # in dev mode. we use git to get a "version number"
        # that will change with tags and commits
        process = subprocess.run(
            ["git", "describe", "--tags", "--dirty", "--long"],
            capture_output=True,
            universal_newlines=True,
            cwd=Path(__file__).parent.parent,
        )
        if process.returncode == 0:
            __version__ = process.stdout.strip()
