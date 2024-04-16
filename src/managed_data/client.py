from dataclasses import dataclass
from typing import Union

import requests
from beartype import beartype
from deeporigin import cache_do_api_tokens, get_do_api_tokens
from deeporigin.config import get_value
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils import _nucleus_url


@dataclass
class DeepOriginClient:
    api_url = _nucleus_url()
    org_id = get_value()["organization_id"]

    api_access_token, api_refresh_token = get_do_api_tokens()
    cache_do_api_tokens(api_access_token, api_refresh_token)

    headers = {
        "accept": "application/json",
        "authorization": f"Bearer {api_access_token}",
        "content-type": "application/json",
        "x-org-id": org_id,
    }

    @beartype
    def invoke(
        self,
        endpoint: str,
        data: dict,
    ) -> Union[dict, list]:
        """core call to API endpoint"""

        response = requests.post(
            f"{self.api_url}{endpoint}",
            headers=self.headers,
            json=data,
        )

        return _check_response(response)


@beartype
def _check_response(response: requests.models.Response) -> Union[dict, list]:
    """utility function to check responses"""

    if response.status_code == 404:
        raise DeepOriginException("[Error 404] The requested resource was not found.")

    response.raise_for_status()
    response = response.json()

    if "error" in response:
        raise DeepOriginException(response["error"])

    if "data" in response:
        return response["data"]
    else:
        raise KeyError("`data` not in response")
