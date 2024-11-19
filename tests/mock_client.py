"""this module provides a mock client for testing purposes"""

import json
from pathlib import Path


class MockClient:
    """mock client to respond with static data for testing
    purposes"""

    def __getattr__(self, name):
        """general purpose catch all for all methods"""

        def method(*args, **kwargs):
            return self._return_response(name, *args, **kwargs)

        return method

    def _return_response(self, name, *args, **kwargs):
        """return stashed values, mimicking responses from the live instance exactly"""

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
