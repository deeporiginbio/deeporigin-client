"""this module provides a mock client for testing purposes"""

from pathlib import Path

from deeporigin.platform import Client
from deeporigin.platform.recording import RequestRecorder


class MockClient(Client):
    """mock client to respond using recorded data stored in SQLite or JSON fixtures.

    If a SQLite database exists at tests/fixtures/requests.sqlite (default),
    responses will be served from there using sequencing (0,1,2,...) per
    request signature. Otherwise, falls back to legacy JSON fixtures in
    tests/fixtures/responses/<method>.json.
    """

    is_mock = True

    def __getattr__(self, name):
        """general purpose catch all for all methods"""

        # Don't intercept private attributes that should be handled normally
        if name.startswith("_"):
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no attribute '{name}'"
            )

        def method(*args, **kwargs):
            return self._return_response(name, *args, **kwargs)

        return method

    def _return_response(self, name, *args, **kwargs):
        """return stashed values, mimicking responses from the live instance exactly"""

        print(f"Returning response for {name} with kwargs: {kwargs}")

        # Prefer SQLite recordings when available
        sqlite_path = (
            Path(__file__).resolve().parent.parent
            / "tests"
            / "fixtures"
            / "requests.sqlite"
        )

        recorder = RequestRecorder(sqlite_path)
        # Maintain per-process state for sequencing
        if not hasattr(self, "_seq_state"):
            self._seq_state = {}
        try:
            data = recorder.fetch_next(
                method=name,
                kwargs=kwargs,
                state=self._seq_state,
            )
            return {"data": data}
        except KeyError:
            # If SQLite has no matching interaction, fall back to JSON fixtures for compatibility
            pass
