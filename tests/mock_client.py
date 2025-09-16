"""this module provides a mock client for testing purposes"""

import json
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

    def __getattr__(self, name):
        """general purpose catch all for all methods"""

        def method(*args, **kwargs):
            return self._return_response(name, *args, **kwargs)

        return method

    def _return_response(self, name, *args, **kwargs):
        """return stashed values, mimicking responses from the live instance exactly"""

        # Prefer SQLite recordings when available
        sqlite_path = (
            Path(__file__).resolve().parent.parent
            / "tests"
            / "fixtures"
            / "requests.sqlite"
        )

        if sqlite_path.exists():
            recorder = RequestRecorder(sqlite_path)
            # Maintain per-process state for sequencing
            if not hasattr(self, "_seq_state"):
                self._seq_state = {}
            try:
                data = recorder.fetch_next(
                    method=name, kwargs=kwargs, state=self._seq_state
                )
                return dict(data=data)
            except KeyError:
                # If SQLite has no matching interaction, fall back to JSON fixtures for compatibility
                pass

        stash_loc = (
            Path(__file__).resolve().parent.parent
            / "tests"
            / "fixtures"
            / "responses"
            / f"{name}.json"
        )

        if not stash_loc.exists():
            raise FileNotFoundError(
                f"Stash file for method: {name} not found. You may need to generate it by running tests against a live instance."
            )

        with open(stash_loc) as f:
            data = json.load(f)

        key = json.dumps(kwargs, sort_keys=True)
        if key in data.keys():
            return dict(data=data[key])
        else:
            print("Available keys are: ")
            for _key in data.keys():
                print(_key)
            raise KeyError(
                f"Could not find key called `{key}` for function `{name}`. You may need to generate it by running tests against a live instance using the same arguments."
            )
