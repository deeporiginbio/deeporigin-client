import unittest
import unittest.mock

from deeporigin import config, context, feature_flags


class TestCase(unittest.TestCase):
    def setUp(self):
        config.get_value.cache_clear()
        context.get_value.cache_clear()
        feature_flags.get_value.cache_clear()

    def tearDown(self):
        config.get_value.cache_clear()
        context.get_value.cache_clear()
        feature_flags.get_value.cache_clear()

    def test_get_value(self):
        env = {
            "DEEP_ORIGIN_ORGANIZATION_ID": "a123",
            "DEEP_ORIGIN_BENCH_ID": "b123",
            "DEEP_ORIGIN_ENV": "env123",
            "DEEP_ORIGIN_API_ENDPOINT": "c123",
            "DEEP_ORIGIN_USER_ID": "auth0|xxxxxxxxxxxxxxxxxxxxxxxx",
            "DEEP_ORIGIN_AUTH_DOMAIN": "xxx",
            "DEEP_ORIGIN_AUTH_DEVICE_CODE_ENDPOINT": "xxx",
            "DEEP_ORIGIN_AUTH_TOKEN_ENDPOINT": "xxx",
            "DEEP_ORIGIN_AUTH_AUDIENCE": "xxx",
            "DEEP_ORIGIN_AUTH_GRANT_TYPE": "xxx",
            "DEEP_ORIGIN_AUTH_CLIENT_ID": "xxx",
            "DEEP_ORIGIN_AUTH_CLIENT_SECRET": "xxx",
            "DEEP_ORIGIN_LIST_BENCH_VARIABLES_QUERY_TEMPLATE": "xxx",
        }

        feature_flags.get_value.cache_clear()
        with unittest.mock.patch.dict("os.environ", env):
            value = feature_flags.get_value()
        expected_value = feature_flags.FeatureFlags(
            variables=True,
        )
        self.assertEqual(expected_value, value)

    def test_get_value_env_var_override(self):
        env = {
            "DEEP_ORIGIN_ORGANIZATION_ID": "a123",
            "DEEP_ORIGIN_BENCH_ID": "b123",
            "DEEP_ORIGIN_ENV": "env123",
            "DEEP_ORIGIN_API_ENDPOINT": "c123",
            "DEEP_ORIGIN_USER_ID": "auth0|xxxxxxxxxxxxxxxxxxxxxxxx",
            "DEEP_ORIGIN_AUTH_DOMAIN": "xxx",
            "DEEP_ORIGIN_AUTH_DEVICE_CODE_ENDPOINT": "xxx",
            "DEEP_ORIGIN_AUTH_TOKEN_ENDPOINT": "xxx",
            "DEEP_ORIGIN_AUTH_AUDIENCE": "xxx",
            "DEEP_ORIGIN_AUTH_GRANT_TYPE": "xxx",
            "DEEP_ORIGIN_AUTH_CLIENT_ID": "xxx",
            "DEEP_ORIGIN_AUTH_CLIENT_SECRET": "xxx",
            "DEEP_ORIGIN_LIST_BENCH_VARIABLES_QUERY_TEMPLATE": "xxx",
            "DEEP_ORIGIN_FEATURE_FLAGS__VARIABLES": "false",
        }

        config.get_value.cache_clear()
        feature_flags.get_value.cache_clear()
        with unittest.mock.patch.dict("os.environ", env):
            value = feature_flags.get_value()
        expected_value = feature_flags.FeatureFlags(
            variables=False,
        )
        self.assertEqual(expected_value, value)

        env["DEEP_ORIGIN_FEATURE_FLAGS__VARIABLES"] = "true"
        config.get_value.cache_clear()
        feature_flags.get_value.cache_clear()
        with unittest.mock.patch.dict("os.environ", env):
            value = feature_flags.get_value()
        expected_value.variables = True
        self.assertEqual(expected_value, value)
