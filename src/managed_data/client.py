import json
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, Union
from urllib.parse import urljoin

import requests
from beartype import beartype
from deeporigin import auth
from deeporigin.exceptions import DeepOriginException


@dataclass
class Client(ABC):
    """client abstract base class. Actual clients inherit from this base class. Actual clients need to implement a `authenticate` and `invoke` method"""

    @abstractmethod
    def authenticate(self, refresh_tokens: bool = True):
        """authenticate method to be called before making any requests. This needs to be implemented in derived classes"""

        pass  # pragma: no cover

    @abstractmethod
    def invoke(self, endpoint: str, data: dict):
        """invoke method to be called to make requests. This needs to be implemented in derived classes"""

        pass  # pragma: no cover

    @abstractmethod
    def download(self, endpoint: str, data: dict):
        """abstract method to download a resource to disk"""

        pass  # pragma: no cover

    @abstractmethod
    def put(self, endpoint: str, data: dict):
        """abstract method to implement a PUT request"""

        pass  # pragma: no cover

    def __props__(self):
        """helper method to show all properties of the client. Do not use."""
        props = dict()
        for attribute in dir(self):
            if attribute.startswith("_"):
                continue

            try:
                if not callable(getattr(self, attribute)):
                    props[str(attribute)] = getattr(self, attribute)
            except Exception:
                pass

        return json.dumps(props, indent=4)

    def __repr__(self):
        """helper method to show all properties of the client. Do not use."""
        return self.__props__()

    def _repr_html_(self):
        """helper method to show all properties of the client. Do not use. This method is called when you invoke a variable in a Jupyter instance"""
        return self.__props__()

    def __str__(self):
        """helper method to show all properties of the client. Do not use."""
        return self.__props__()


@dataclass
class DeepOriginClient(Client):
    """client to interact with DeepOrigin API"""

    headers = dict()

    def __init__(
        self,
        *,
        api_url: Optional[str] = None,
        org_id: Optional[str] = None,
    ):
        """custom init method to create a DeepOriginClient without needing config files"""
        if not api_url:
            from deeporigin.utils import _nucleus_url

            api_url = _nucleus_url()

        if not org_id:
            from deeporigin.config import get_value

            org_id = get_value()["organization_id"]

        super().__init__()

    def authenticate(
        self,
        *,
        refresh_tokens: bool = True,
        access_token: Optional[str] = None,
    ):
        """authenticate to DeepOrigin API

        Authenticate to Deep Origin API. This needs to be called
        before making any requests. refresh_tokens=False allows
        for turning off refresh (saving one network request on
        startup, ), making the CLI faster

        Args:
            refresh_tokens (bool, optional): Whether to
            refresh tokens. Defaults to True.
            access_token: if provided, use this token instead of asking the authentication module for a token


        """

        if not access_token:
            tokens = auth.get_tokens(refresh=refresh_tokens)
            access_token = tokens["access"]

        self.headers = {
            "accept": "application/json",
            "authorization": f"Bearer {access_token}",
            "content-type": "application/json",
            "x-org-id": self.org_id,
        }

    def put(self, *args, **kwargs):
        """overridden put method to use requests to make PUT requests"""

        return requests.put(*args, **kwargs)

    def download(self, url: str, save_path: str) -> None:
        """concrete method to download a resource using GET and save to disk

        Args:
            url (str): url to download
            save_path (str): path to save file
        """

        with requests.get(url, stream=True) as response:
            if response.status_code != 200:
                raise DeepOriginException(f"Failed to download file from {url}")

            with open(save_path, "wb") as file:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:  # Filter out keep-alive new chunks
                        file.write(chunk)

    @beartype
    def invoke(
        self,
        endpoint: str,
        data: dict,
    ) -> Union[dict, list]:
        """core call to API endpoint"""

        url = urljoin(self.api_url, endpoint)

        response = requests.post(
            url,
            headers=self.headers,
            json=data,
        )

        try:
            response = _check_response(response)
        except Exception as e:
            print(f"URL: {url}")
            print(f"headers: {self.headers}")
            print(f"data: {data}")

            raise e

        return response


@beartype
def _check_response(response: requests.models.Response) -> Union[dict, list]:
    """utility function to check responses"""

    if response.status_code == 404:
        raise DeepOriginException(
            f"[Error 404] The requested resource was not found. The response was: {response.json()}"
        )
    elif response.status_code == 400:
        raise DeepOriginException(f"[Error 400] The response was: {response.json()}")

    response.raise_for_status()
    response = response.json()

    if "error" in response:
        raise DeepOriginException(response["error"])

    if "data" in response:
        return response["data"]
    else:
        raise KeyError("`data` not in response")
