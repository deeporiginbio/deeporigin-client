import io
from contextlib import redirect_stderr, redirect_stdout

from beartype import beartype
from deeporigin import cli


@beartype
def _run_cli_command(argv: list[str], client) -> str:
    """helper function to run a CLI command, parse output and return"""
    stdout = io.StringIO()
    stderr = io.StringIO()

    with redirect_stdout(stdout), redirect_stderr(stderr):
        with cli.App(argv=argv) as app:
            app.client = client
            app.run()

    return stdout.getvalue().strip()
