import os
import sys
import tempfile
import unittest
import unittest.mock

import yaml
from deeporigin import config
from deeporigin.exceptions import DeepOriginException


class TestCase(unittest.TestCase):
    def setUp(self):
        config.get_value.cache_clear()

    def tearDown(self):
        config.get_value.cache_clear()

    def test_env_vars(self):
        env = {
            "DEEP_ORIGIN_ORGANIZATION_ID": "a123",
            "DEEP_ORIGIN_BENCH_ID": "b123",
            "DEEP_ORIGIN_ENV": "env123",
            "DEEP_ORIGIN_API_ENDPOINT": "abc",
            "DEEP_ORIGIN_NUCLEUS_API_ROUTE": "nar",
            "DEEP_ORIGIN_AUTH_DOMAIN": "def",
            "DEEP_ORIGIN_AUTH_DEVICE_CODE_ENDPOINT": "def",
            "DEEP_ORIGIN_AUTH_TOKEN_ENDPOINT": "def",
            "DEEP_ORIGIN_AUTH_AUDIENCE": "ghi",
            "DEEP_ORIGIN_AUTH_GRANT_TYPE": "jklm",
            "DEEP_ORIGIN_AUTH_CLIENT_ID": "mno",
            "DEEP_ORIGIN_AUTH_CLIENT_SECRET": "mnop1",
            "DEEP_ORIGIN_LIST_BENCH_VARIABLES_QUERY_TEMPLATE": "mnop",
            "DEEP_ORIGIN_API_TOKENS_FILENAME": "mnopq",
            "DEEP_ORIGIN_VARIABLES_CACHE_FILENAME": "vwx",
            "DEEP_ORIGIN_AUTO_INSTALL_VARIABLES_FILENAME": "vwx",
            "DEEP_ORIGIN_GRAPHQL_API_ROUTE": "foo",
        }
        config.get_value.cache_clear()
        with unittest.mock.patch.dict("os.environ", env):
            value = config.get_value(user_config_filenames=tuple())
        expected_value = {
            "organization_id": "a123",
            "bench_id": "b123",
            "env": "env123",
            "api_endpoint": "abc",
            "nucleus_api_route": "nar",
            "auth_domain": "def",
            "auth_device_code_endpoint": "def",
            "auth_token_endpoint": "def",
            "auth_audience": "ghi",
            "graphql_api_route": "foo",
            "auth_grant_type": "jklm",
            "auth_client_id": "mno",
            "auth_client_secret": "mnop1",
            "list_bench_variables_query_template": "mnop",
            "api_tokens_filename": os.path.abspath("mnopq"),
            "variables_cache_filename": os.path.abspath("vwx"),
            "auto_install_variables_filename": os.path.abspath("vwx"),
            "feature_flags": None,
        }

        self.assertEqual(expected_value, value)

        env["DEEP_ORIGIN_FEATURE_FLAGS__VARIABLES"] = "false"
        config.get_value.cache_clear()
        with unittest.mock.patch.dict("os.environ", env):
            value = config.get_value(user_config_filenames=tuple())
        expected_value["feature_flags"] = {"variables": False}
        self.assertEqual(expected_value, value)

        env["DEEP_ORIGIN_FEATURE_FLAGS__VARIABLES"] = "true"
        config.get_value.cache_clear()
        with unittest.mock.patch.dict("os.environ", env):
            value = config.get_value(user_config_filenames=tuple())
        expected_value["feature_flags"]["variables"] = True
        self.assertEqual(expected_value, value)

    def test_file(self):
        env = {}

        user_config = {
            "organization_id": "a123",
            "bench_id": "b123",
            "env": "env123",
            "api_endpoint": "abc",
            "nucleus_api_route": "nar",
            "graphql_api_route": "foo",
            "auth_domain": "def",
            "auth_device_code_endpoint": "def",
            "auth_token_endpoint": "def",
            "auth_audience": "ghi",
            "auth_grant_type": "jklm",
            "auth_client_id": "mno",
            "auth_client_secret": "mnop1",
            "list_bench_variables_query_template": "mnop",
            "api_tokens_filename": os.path.join("~", "mnopq"),
            "variables_cache_filename": os.path.join("~", "vwx"),
            "auto_install_variables_filename": os.path.join("~", "vwx"),
        }
        _, user_config_filename = tempfile.mkstemp()
        with open(user_config_filename, "w") as file:
            yaml.dump(user_config, file, Dumper=yaml.CDumper)

        with unittest.mock.patch.dict("os.environ", env):
            value = config.get_value(user_config_filenames=(user_config_filename,))

        user_config["api_tokens_filename"] = os.path.expanduser(
            user_config["api_tokens_filename"]
        )
        user_config["variables_cache_filename"] = os.path.expanduser(
            user_config["variables_cache_filename"]
        )
        user_config["auto_install_variables_filename"] = os.path.expanduser(
            user_config["auto_install_variables_filename"]
        )
        user_config["feature_flags"] = None

        self.assertEqual(user_config, value)

        if not sys.platform.startswith("win"):
            os.remove(user_config_filename)

    def test_invalid(self):
        env = {}

        user_config = {
            "organization_id": "a123",
            "bench_id": "b123",
            "env": "env123",
            "api_endpoint": "abc",
            "auth_domain": "def",
            "auth_device_code_endpoint": "def",
            "auth_token_endpoint": "def",
            "auth_audience": "ghi",
            "auth_grant_type": "jklm",
            "auth_client_id": None,
            "auth_client_secret": "mnop1",
            "list_bench_variables_query_template": "mnop",
            "api_tokens_filename": os.path.join("~", "mnopq"),
            "variables_cache_filename": os.path.join("~", "vwx"),
            "auto_install_variables_filename": os.path.join("~", "vwx"),
        }
        _, user_config_filename = tempfile.mkstemp()
        with open(user_config_filename, "w") as file:
            yaml.dump(user_config, file, Dumper=yaml.CDumper)

        with unittest.mock.patch.dict("os.environ", env):
            with self.assertRaisesRegex(
                DeepOriginException, "configuration is not valid"
            ):
                config.get_value(user_config_filenames=(user_config_filename,))
