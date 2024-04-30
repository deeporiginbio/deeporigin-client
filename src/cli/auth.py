import cement
from deeporigin import cache_do_api_tokens, get_do_api_tokens

__all__ = [
    "CONTROLLERS",
    "AuthController",
]


class AuthController(cement.Controller):
    class Meta:
        label = "authenticate"
        stacked_on = "base"
        stacked_type = "nested"
        help = "Authenticate to the Deep Origin platform"
        description = "Authenticate this client to the Deep Origin platform and save access tokens for later use."
        arguments = []

    @cement.ex(hide=True)
    def _default(self):
        access_token, refresh_token = get_do_api_tokens(verbose=True)

        cache_do_api_tokens(access_token, refresh_token)


CONTROLLERS = [
    AuthController,
]
