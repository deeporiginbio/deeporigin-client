import copy
import io
import os
import pathlib
import shutil
import sys
import tempfile
import unittest
import unittest.mock
from contextlib import redirect_stderr, redirect_stdout

import crontab
import pydantic
from deeporigin import auth, cli, config, feature_flags, utils, variables
from deeporigin.exceptions import DeepOriginException
from deeporigin.variables import core as variables_core
from deeporigin.variables import types as variables_types
from deeporigin.warnings import DeepOriginWarning


@unittest.skipIf(sys.platform.startswith("win"), "Test skipped on Windows")
class TestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(
            os.path.join(os.path.dirname(__file__), "fixtures", "ssh_key"), "r"
        ) as file:
            cls.TEST_SSH_PRIVATE_KEY = file.read()

        with open(
            os.path.join(os.path.dirname(__file__), "fixtures", "ssh_key.pub"), "r"
        ) as file:
            cls.TEST_SSH_PUBLIC_KEY = file.read()

        with open(
            os.path.join(os.path.dirname(__file__), "fixtures", "private-key.gpg"), "r"
        ) as file:
            cls.TEST_PRIVATE_GPG_KEY = file.read()

    def setUp(self):
        config.get_value.cache_clear()
        feature_flags.get_value.cache_clear()
        auth.get_tokens.cache_clear()

        _, self.api_tokens_filename = tempfile.mkstemp(suffix=".yml")
        os.remove(self.api_tokens_filename)

        _, self.variables_cache_filename = tempfile.mkstemp(suffix=".yml")
        os.remove(self.variables_cache_filename)

        _, self.auto_install_variables_filename = tempfile.mkstemp(suffix=".sh")
        os.remove(self.auto_install_variables_filename)

        env = {
            "DEEP_ORIGIN_ORGANIZATION_ID": "a123",
            "DEEP_ORIGIN_BENCH_ID": "b123",
            "DEEP_ORIGIN_ENV": "env123",
            "DEEP_ORIGIN_API_ENDPOINT": "c123",
            "DEEP_ORIGIN_API_TOKENS_FILENAME": self.api_tokens_filename,
            "DEEP_ORIGIN_VARIABLES_CACHE_FILENAME": self.variables_cache_filename,
            "DEEP_ORIGIN_AUTO_INSTALL_VARIABLES_FILENAME": self.auto_install_variables_filename,
            "DEEP_ORIGIN_USER_ID": "auth0|xxxxxxxxxxxxxxxxxxxxxxxx",
            "DEEP_ORIGIN_AUTH_DOMAIN": "xxx",
            "DEEP_ORIGIN_AUTH_DEVICE_CODE_ENDPOINT": "xxx",
            "DEEP_ORIGIN_AUTH_TOKEN_ENDPOINT": "xxx",
            "DEEP_ORIGIN_AUTH_AUDIENCE": "xxx",
            "DEEP_ORIGIN_LIST_BENCH_VARIABLES_QUERY_TEMPLATE": "xxx",
            "DEEP_ORIGIN_AUTH_GRANT_TYPE": "xxx",
            "DEEP_ORIGIN_AUTH_CLIENT_ID": "xxx",
            "DEEP_ORIGIN_AUTH_CLIENT_SECRET": "xxx",
        }
        with unittest.mock.patch.dict("os.environ", env):
            config.get_value()
            feature_flags.get_value()

    def tearDown(self):
        config.get_value.cache_clear()
        feature_flags.get_value.cache_clear()
        auth.get_tokens.cache_clear()
        if os.path.isfile(self.api_tokens_filename):
            os.remove(self.api_tokens_filename)
        if os.path.isfile(self.variables_cache_filename):
            os.remove(self.variables_cache_filename)
        if os.path.isfile(self.auto_install_variables_filename):
            os.remove(self.auto_install_variables_filename)

    def test_get_variable_types_by_values(self):
        types = variables.get_variable_types_by_values(
            ["EnvironmentVariable", "SecretEnvironmentVariable"]
        )
        self.assertEqual(
            [
                variables.VariableType.EnvironmentVariable,
                variables.VariableType.SecretEnvironmentVariable,
            ],
            types,
        )

    def test_get_variable_types_by_values_error(self):
        with self.assertRaisesRegex(
            DeepOriginException,
            "The following types of variables and secrets are not valid",
        ):
            variables.get_variable_types_by_values(["NotDefined"])

    def test_deserialize_variable(self):
        self.assertEqual(
            variables_types.EnvironmentVariable(
                drn="drn:0",
                name="a",
                key="b",
                value="c",
            ),
            variables_core.deserialize_variable(
                {
                    "drn": "drn:0",
                    "type": "Non-secret environment variables",
                    "name": "a",
                    "key": "b",
                    "value": "c",
                }
            ),
        )

        self.assertEqual(
            variables_types.File(
                drn="drn:1",
                name="a",
                filename="/b",
                value="c",
            ),
            variables_core.deserialize_variable(
                {
                    "drn": "drn:1",
                    "type": "Configuration file",
                    "name": "a",
                    "filename": "/b",
                    "value": "c",
                }
            ),
        )

        self.assertEqual(
            variables_types.SecretEnvironmentVariable(
                drn="drn:2",
                name="a",
                key="b",
                value="c",
            ),
            variables_core.deserialize_variable(
                {
                    "drn": "drn:2",
                    "type": "Secret environment variables",
                    "name": "a",
                    "key": "b",
                    "value": "c",
                }
            ),
        )

        self.assertEqual(
            variables_types.SecretFile(
                drn="drn:3",
                name="a",
                filename="/b",
                value="c",
            ),
            variables_core.deserialize_variable(
                {
                    "drn": "drn:3",
                    "type": "Secret file",
                    "name": "a",
                    "filename": "/b",
                    "value": "c",
                }
            ),
        )

        self.assertEqual(
            variables_types.AwsProfile(
                drn="drn:4",
                name="a",
                profile="b",
                access_key_id="C" * 20,
                secret_access_key="D" * 40,
                region="us-west-1",
            ),
            variables_core.deserialize_variable(
                {
                    "drn": "drn:4",
                    "type": "AWS Credentials",
                    "name": "a",
                    "profileName": "b",
                    "accessKey": "C" * 20,
                    "secretKey": "D" * 40,
                    "defaultRegion": "us-west-1",
                }
            ),
        )

        self.assertEqual(
            variables_types.GitHttpCredentials(
                drn="drn:5",
                name="a",
                repository="github.com/org/repo-b",
                username="c",
                password="d",
            ),
            variables_core.deserialize_variable(
                {
                    "drn": "drn:5",
                    "type": "Git HTTP credentials",
                    "name": "a",
                    "repository": "github.com/org/repo-b",
                    "username": "c",
                    "password": "d",
                }
            ),
        )

        self.assertEqual(
            variables_types.PrivateGpgKey(
                drn="drn:6",
                name="a",
                value=self.TEST_PRIVATE_GPG_KEY,
            ),
            variables_core.deserialize_variable(
                {
                    "drn": "drn:6",
                    "type": "GPG signing keys",
                    "name": "a",
                    "value": self.TEST_PRIVATE_GPG_KEY,
                }
            ),
        )

        self.assertEqual(
            variables_types.PrivateSshKey(
                drn="drn:7",
                name="a",
                filename="/b",
                value=self.TEST_SSH_PRIVATE_KEY,
            ),
            variables_core.deserialize_variable(
                {
                    "drn": "drn:7",
                    "type": "SSH private keys",
                    "name": "a",
                    "filename": "/b",
                    "value": self.TEST_SSH_PRIVATE_KEY,
                }
            ),
        )

        with self.assertRaisesRegex(DeepOriginException, "not a valid type"):
            (
                variables_core.deserialize_variable(
                    {
                        "drn": "drn:8",
                        "type": "Unknown",
                        "name": "",
                        "filename": "/b",
                        "value": "",
                    }
                ),
            )

    def test_get_variables_from_do_platform(self):
        config.get_value.cache_clear()
        feature_flags.get_value.cache_clear()
        auth.get_tokens.cache_clear()

        env = {
            "DEEP_ORIGIN_ORGANIZATION_ID": "a123",
            "DEEP_ORIGIN_BENCH_ID": "b123",
            "DEEP_ORIGIN_ENV": "env123",
            "DEEP_ORIGIN_API_ENDPOINT": "c123",
            "DEEP_ORIGIN_API_TOKENS_FILENAME": self.api_tokens_filename,
            "DEEP_ORIGIN_VARIABLES_CACHE_FILENAME": self.variables_cache_filename,
            "DEEP_ORIGIN_AUTO_INSTALL_VARIABLES_FILENAME": self.auto_install_variables_filename,
            "DEEP_ORIGIN_USER_ID": "auth0|xxxxxxxxxxxxxxxxxxxxxxxx",
            "DEEP_ORIGIN_AUTH_DOMAIN": "xxx",
            "DEEP_ORIGIN_AUTH_DEVICE_CODE_ENDPOINT": "xxx",
            "DEEP_ORIGIN_AUTH_TOKEN_ENDPOINT": "xxx",
            "DEEP_ORIGIN_AUTH_AUDIENCE": "xxx",
            "DEEP_ORIGIN_LIST_BENCH_VARIABLES_QUERY_TEMPLATE": "xxx",
            "DEEP_ORIGIN_AUTH_GRANT_TYPE": "xxx",
            "DEEP_ORIGIN_AUTH_CLIENT_ID": "xxx",
            "DEEP_ORIGIN_AUTH_CLIENT_SECRET": "xxx",
        }
        with unittest.mock.patch.dict("os.environ", env):
            config.get_value()
            feature_flags.get_value()

        responses = [
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "device_code": "xxx",
                    "user_code": "xxx",
                    "verification_uri_complete": "xxx",
                    "interval": 5,
                },
            ),
            unittest.mock.MagicMock(
                status_code=200,
                json=lambda: {
                    "access_token": "xxx",
                    "refresh_token": "xxx",
                },
            ),
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "data": {
                        "listBenchSecrets": [
                            {
                                "drn": "drn:1",
                                "type": "Non-secret environment variables",
                                "name": "abc",
                                "key": "def",
                                "value": "ghi",
                            },
                            {
                                "drn": "drn:2",
                                "type": "Secret file",
                                "name": "jkl",
                                "filename": "~/mno",
                                "value": "pqr",
                            },
                        ],
                    },
                },
            ),
        ]
        with unittest.mock.patch("requests.post", side_effect=responses):
            vars = variables_core.get_variables_from_do_platform()
            expected_vars = [
                variables_types.EnvironmentVariable(
                    drn="drn:1",
                    name="abc",
                    key="def",
                    value="ghi",
                ),
                variables_types.SecretFile(
                    drn="drn:2",
                    name="jkl",
                    filename="~/mno",
                    value="pqr",
                ),
            ]
            self.assertEqual(expected_vars, vars)

    def test_get_export_variables_from_local(self):
        if os.path.isfile(self.api_tokens_filename):
            os.remove(self.api_tokens_filename)

        if os.path.isfile(self.variables_cache_filename):
            os.remove(self.variables_cache_filename)

        if os.path.isfile(self.auto_install_variables_filename):
            os.remove(self.auto_install_variables_filename)

        self.assertEqual({}, variables_core.get_variables_from_local())

        vars = {
            "EnvironmentVariable": {
                "abc": {
                    "name": "var-1",
                    "key": "abc",
                    "value": "123",
                },
                "def": {
                    "name": "var-2",
                    "key": "def",
                    "value": "456",
                },
            }
        }

        variables_core.export_variables_to_local(vars)

        self.assertEqual(vars, variables_core.get_variables_from_local())

    def test_get_variable_local_key(self):
        self.assertEqual(
            "drn:1",
            variables_core.get_variable_local_key(
                variables_types.EnvironmentVariable(
                    drn="drn:1", name="var", key="abc", value=""
                )
            ),
        )
        self.assertEqual(
            "drn:2",
            variables_core.get_variable_local_key(
                variables_types.SecretEnvironmentVariable(
                    drn="drn:2", name="var", key="abc", value=""
                )
            ),
        )
        self.assertEqual(
            "drn:3",
            variables_core.get_variable_local_key(
                variables_types.File(drn="drn:3", name="var", filename="/abc", value="")
            ),
        )
        self.assertEqual(
            "drn:4",
            variables_core.get_variable_local_key(
                variables_types.SecretFile(
                    drn="drn:4", name="var", filename="/abc", value=""
                )
            ),
        )
        self.assertEqual(
            "drn:5",
            variables_core.get_variable_local_key(
                variables_types.PrivateSshKey(
                    drn="drn:5",
                    name="var",
                    filename="/abc",
                    value=self.TEST_SSH_PRIVATE_KEY,
                )
            ),
        )
        self.assertEqual(
            "drn:6",
            variables_core.get_variable_local_key(
                variables_types.PrivateGpgKey(
                    drn="drn:6", name="var", value=self.TEST_PRIVATE_GPG_KEY
                )
            ),
        )
        self.assertEqual(
            "drn:7",
            variables_core.get_variable_local_key(
                variables_types.AwsProfile(
                    drn="drn:7",
                    name="var",
                    profile="abc",
                    access_key_id="X" * 20,
                    secret_access_key="Y" * 40,
                    region="us-west-1",
                )
            ),
        )
        self.assertEqual(
            "drn:8",
            variables_core.get_variable_local_key(
                variables_types.GitHttpCredentials(
                    drn="drn:8",
                    name="var",
                    repository="github.com/org/repo-abc",
                    username="def",
                    password="",
                )
            ),
        )

    def test_is_variable_modified(self):
        if os.path.isfile(self.api_tokens_filename):
            os.remove(self.api_tokens_filename)

        if os.path.isfile(self.variables_cache_filename):
            os.remove(self.variables_cache_filename)

        if os.path.isfile(self.auto_install_variables_filename):
            os.remove(self.auto_install_variables_filename)

        var = variables_types.EnvironmentVariable(
            drn="drn:1", name="var-1", key="abc", value="123"
        )
        self.assertTrue(variables_core.is_variable_modified(var))
        self.assertFalse(variables_core.is_variable_modified(var))
        var = variables_types.EnvironmentVariable(
            drn="drn:1", name="var-1", key="abc", value="456"
        )
        self.assertTrue(variables_core.is_variable_modified(var))
        var = variables_types.EnvironmentVariable(
            drn="drn:1", name="var-1", key="def", value="456"
        )
        self.assertTrue(variables_core.is_variable_modified(var))
        self.assertFalse(variables_core.is_variable_modified(var))
        var = variables_types.SecretEnvironmentVariable(
            drn="drn:1", name="var-1", key="def", value="456"
        )
        self.assertTrue(variables_core.is_variable_modified(var))

        var = variables_types.SecretFile(
            drn="drn:2", name="var-2", filename="/abc", value="123"
        )
        self.assertTrue(variables_core.is_variable_modified(var))
        self.assertFalse(variables_core.is_variable_modified(var))
        var = variables_types.SecretFile(
            drn="drn:2", name="var-2", filename="/abc", value="456"
        )
        self.assertTrue(variables_core.is_variable_modified(var))
        var = variables_types.SecretFile(
            drn="drn:2", name="var-2", filename="/def", value="456"
        )
        self.assertTrue(variables_core.is_variable_modified(var))
        self.assertFalse(variables_core.is_variable_modified(var))
        var = variables_types.File(
            drn="drn:2", name="var-2", filename="/def", value="456"
        )
        self.assertTrue(variables_core.is_variable_modified(var))

    def test_install_environment_variable_variable(self):
        user_home_dirname = tempfile.mkdtemp()
        env_filename = variables_types.EnvironmentVariable.get_env_filename(
            user_home_dirname=user_home_dirname
        )

        self.assertEqual(
            {},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            ),
        )

        var = variables_types.SecretEnvironmentVariable(
            drn="drn:1", name="var-1", key="var_1", value="val_1"
        )
        var.install(None, user_home_dirname)
        self.assertEqual(
            {"var_1": "val_1"},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            ),
        )

        var.install(None, user_home_dirname)
        self.assertEqual(
            {"var_1": "val_1"},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            ),
        )

        var.value = 'val\n"2'
        var.install(None, user_home_dirname)
        self.assertEqual(
            {"var_1": 'val\n"2'},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            ),
        )

        bashrc_filename = utils.expand_user(
            os.path.join("~", ".bashrc"), user_home_dirname
        )
        with open(bashrc_filename, "r") as file:
            lines = file.read().split("\n")
        self.assertEqual(10, len(lines))
        self.assertEqual("######## BEGIN DEEP ORIGIN CLI ########", lines[-8])
        self.assertEqual("######## END DEEP ORIGIN CLI ########", lines[-2])

        var.uninstall(user_home_dirname)
        self.assertEqual(
            {},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            ),
        )

        var.install(None, user_home_dirname)
        variables_core.export_variables_to_local(
            {
                var.drn: dict(type=var.__class__.__name__, **var.__dict__),
            }
        )
        self.assertEqual(
            set([var.drn]), set(variables_core.get_variables_from_local().keys())
        )
        self.assertEqual(
            {"var_1": 'val\n"2'},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            ),
        )

        variables.uninstall_variables(
            types=[variables.VariableType.AwsProfile],
            user_home_dirname=user_home_dirname,
        )
        self.assertEqual(
            set([var.drn]), set(variables_core.get_variables_from_local().keys())
        )
        self.assertEqual(
            {"var_1": 'val\n"2'},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            ),
        )

        variables.uninstall_variables(user_home_dirname=user_home_dirname)
        self.assertEqual({}, variables_core.get_variables_from_local())
        self.assertEqual(
            {},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            ),
        )

        shutil.rmtree(user_home_dirname)

    def test_install_environment_variable_variable_direct_modification(self):
        user_home_dirname = tempfile.mkdtemp()
        env_filename = variables_types.EnvironmentVariable.get_env_filename(
            user_home_dirname=user_home_dirname
        )

        self.assertEqual(
            {},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            ),
        )

        var = variables_types.SecretEnvironmentVariable(
            drn="drn:1", name="var-1", key="var_1", value="val_1"
        )
        var.install(None, user_home_dirname)
        self.assertEqual(
            {"var_1": "val_1"},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            ),
        )

        prev_var = copy.deepcopy(var)
        prev_var.value = "val_2"
        var_2 = copy.deepcopy(var)
        var_2.value = "val_3"
        with self.assertWarnsRegex(DeepOriginWarning, "modified directly"):
            var_2.install(prev_var, user_home_dirname)
        self.assertEqual(
            {"var_1": "val_1"},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            ),
        )

        prev_var = copy.deepcopy(var)
        var_2 = copy.deepcopy(var)
        var_2.value = "val_3"
        var_2.install(prev_var, user_home_dirname)
        self.assertEqual(
            {"var_1": "val_3"},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            ),
        )

        shutil.rmtree(user_home_dirname)

    def test_install_file_variable(self):
        user_home_dirname = tempfile.mkdtemp()

        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.PrivateSshKey.read_ssh_config(
            user_home_dirname=user_home_dirname
        )
        self.assertEqual([], other_config)
        self.assertEqual([], do_config)
        self.assertFalse(in_do_config)

        var = variables_types.PrivateSshKey(
            drn="drn:1",
            name="var-1",
            filename=os.path.join("~", ".ssh", "var_1"),
            value=self.TEST_SSH_PRIVATE_KEY,
        )
        var.install(None, user_home_dirname)
        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.PrivateSshKey.read_ssh_config(
            variable=var, user_home_dirname=user_home_dirname
        )
        self.assertEqual(["\n", "\n"], other_config)
        self.assertEqual(1, len(do_config))
        self.assertTrue(do_config[0].startswith("IdentityFile "))
        self.assertTrue(in_do_config)

        var.install(None, user_home_dirname)
        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.PrivateSshKey.read_ssh_config(
            variable=var, user_home_dirname=user_home_dirname
        )
        self.assertEqual(["\n", "\n"], other_config)
        self.assertEqual(1, len(do_config))
        self.assertTrue(do_config[0].startswith("IdentityFile "))
        self.assertTrue(in_do_config)

        var.value = "val_2"
        var.install(None, user_home_dirname)
        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.PrivateSshKey.read_ssh_config(
            variable=var, user_home_dirname=user_home_dirname
        )
        self.assertEqual(["\n", "\n"], other_config)
        self.assertEqual(1, len(do_config))
        self.assertTrue(do_config[0].startswith("IdentityFile "))
        self.assertTrue(in_do_config)

        var.filename = os.path.join("~", ".ssh", "var_2")
        var.install(None, user_home_dirname)
        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.PrivateSshKey.read_ssh_config(
            variable=var, user_home_dirname=user_home_dirname
        )
        self.assertEqual(["\n", "\n"], other_config)
        self.assertEqual(2, len(do_config))
        self.assertTrue(do_config[0].startswith("IdentityFile "))
        self.assertTrue(in_do_config)

        self.assertTrue(
            os.path.isfile(os.path.join(user_home_dirname, ".ssh", "var_2"))
        )
        var.uninstall(user_home_dirname)
        self.assertFalse(
            os.path.isfile(os.path.join(user_home_dirname, ".ssh", "var_2"))
        )
        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.PrivateSshKey.read_ssh_config(
            variable=var, user_home_dirname=user_home_dirname
        )
        self.assertEqual(["\n", "\n"], other_config)
        self.assertEqual(1, len(do_config))
        self.assertTrue(do_config[0].startswith("IdentityFile "))
        self.assertFalse(in_do_config)

        shutil.rmtree(user_home_dirname)

    def test_install_file_variable_direct_modification(self):
        user_home_dirname = tempfile.mkdtemp()

        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.PrivateSshKey.read_ssh_config(
            user_home_dirname=user_home_dirname
        )
        self.assertEqual([], other_config)
        self.assertEqual([], do_config)
        self.assertFalse(in_do_config)

        var = variables_types.PrivateSshKey(
            drn="drn:1",
            name="var-1",
            filename=os.path.join("~", ".ssh", "var_1"),
            value=self.TEST_SSH_PRIVATE_KEY,
        )
        var.install(None, user_home_dirname)
        with open(os.path.join(user_home_dirname, ".ssh", "var_1"), "r") as file:
            self.assertEqual(var.value, file.read())

        prev_var = copy.deepcopy(var)
        prev_var.value = "val_2"
        var_2 = copy.deepcopy(var)
        var_2.value = "val_3"
        with self.assertWarnsRegex(DeepOriginWarning, "modified directly"):
            var_2.install(prev_var, user_home_dirname)
        with open(os.path.join(user_home_dirname, ".ssh", "var_1"), "r") as file:
            self.assertEqual(var.value, file.read())

        prev_var = copy.deepcopy(var)
        var_2 = copy.deepcopy(var)
        var_2.value = "val_3"
        var_2.install(prev_var, user_home_dirname)
        with open(os.path.join(user_home_dirname, ".ssh", "var_1"), "r") as file:
            self.assertEqual(var_2.value, file.read())

        shutil.rmtree(user_home_dirname)

    def test_install_gpg_private_key_variable(self):
        user_home_dirname = tempfile.mkdtemp()

        var = variables_types.PrivateGpgKey(
            drn="drn:1", name="var-1", value=self.TEST_PRIVATE_GPG_KEY
        )
        with unittest.mock.patch(
            "subprocess.run", return_value=unittest.mock.MagicMock(returncode=0)
        ):
            var.install(None, user_home_dirname)

        with unittest.mock.patch(
            "subprocess.run", return_value=unittest.mock.MagicMock(returncode=1)
        ):
            with self.assertRaisesRegex(DeepOriginException, "could not be imported"):
                var.install(None, user_home_dirname)

        with unittest.mock.patch(
            "subprocess.run", return_value=unittest.mock.MagicMock(returncode=0)
        ):
            var.uninstall(user_home_dirname)

        with unittest.mock.patch(
            "subprocess.run", return_value=unittest.mock.MagicMock(returncode=1)
        ):
            with self.assertRaisesRegex(DeepOriginException, "could not be inspected"):
                var.uninstall(user_home_dirname)

        with unittest.mock.patch(
            "subprocess.run",
            side_effect=[
                unittest.mock.MagicMock(returncode=0),
                unittest.mock.MagicMock(returncode=1),
            ],
        ):
            with self.assertRaisesRegex(
                DeepOriginException, "could not be uninstalled"
            ):
                var.uninstall(user_home_dirname)

        shutil.rmtree(user_home_dirname)

    def test_install_aws_credentials_variable(self):
        user_home_dirname = tempfile.mkdtemp()

        var = variables_types.AwsProfile(
            drn="drn:1",
            name="var-1",
            profile="profile_1",
            access_key_id="X" * 20,
            secret_access_key="Y" * 40,
            region="us-west-2",
        )
        with unittest.mock.patch(
            "subprocess.run", return_value=unittest.mock.MagicMock(returncode=0)
        ):
            var.install(None, user_home_dirname)

        with unittest.mock.patch(
            "subprocess.run", return_value=unittest.mock.MagicMock(returncode=1)
        ):
            with self.assertRaisesRegex(DeepOriginException, "could not be imported"):
                var.install(None, user_home_dirname)

        with unittest.mock.patch(
            "subprocess.run", return_value=unittest.mock.MagicMock(returncode=0)
        ):
            var.uninstall(user_home_dirname)

        with unittest.mock.patch(
            "subprocess.run",
            side_effect=[
                unittest.mock.MagicMock(
                    returncode=0,
                    stdout=unittest.mock.MagicMock(decode=lambda: var.access_key_id),
                ),
                unittest.mock.MagicMock(
                    returncode=0,
                    stdout=unittest.mock.MagicMock(
                        decode=lambda: var.secret_access_key
                    ),
                ),
                unittest.mock.MagicMock(
                    returncode=0,
                    stdout=unittest.mock.MagicMock(decode=lambda: var.region),
                ),
                unittest.mock.MagicMock(returncode=1),
            ],
        ):
            with self.assertRaisesRegex(
                DeepOriginException, "could not be uninstalled"
            ):
                var.uninstall(user_home_dirname)

        shutil.rmtree(user_home_dirname)

    def test_install_git_http_credentials_variable(self):
        user_home_dirname = tempfile.mkdtemp()

        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.GitHttpCredentials.read_git_credentials(
            user_home_dirname=user_home_dirname
        )
        self.assertEqual([], other_config)
        self.assertEqual([], do_config)
        self.assertFalse(in_do_config)

        var = variables_types.GitHttpCredentials(
            drn="drn:1",
            name="var-1",
            repository="github.com/org/repo",
            username="john-doe",
            password="YYY",
        )
        var.install(None, user_home_dirname)
        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.GitHttpCredentials.read_git_credentials(
            variable=var, user_home_dirname=user_home_dirname
        )
        self.assertEqual(["\n", "\n"], other_config)
        self.assertEqual(1, len(do_config))
        self.assertTrue(do_config[0].startswith("https://"))
        self.assertTrue(in_do_config)

        var.install(None, user_home_dirname)
        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.GitHttpCredentials.read_git_credentials(
            variable=var, user_home_dirname=user_home_dirname
        )
        self.assertEqual(["\n", "\n"], other_config)
        self.assertEqual(1, len(do_config))
        self.assertTrue(do_config[0].startswith("https://"))
        self.assertTrue(in_do_config)

        var.username = "jane-doe"
        var.install(None, user_home_dirname)
        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.GitHttpCredentials.read_git_credentials(
            variable=var, user_home_dirname=user_home_dirname
        )
        self.assertEqual(["\n", "\n"], other_config)
        self.assertEqual(2, len(do_config))
        self.assertTrue(do_config[0].startswith("https://"))
        self.assertTrue(in_do_config)

        var.uninstall(user_home_dirname)
        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.GitHttpCredentials.read_git_credentials(
            variable=var, user_home_dirname=user_home_dirname
        )
        self.assertEqual(["\n", "\n"], other_config)
        self.assertEqual(1, len(do_config))
        self.assertTrue(do_config[0].startswith("https://"))
        self.assertFalse(in_do_config)

        shutil.rmtree(user_home_dirname)

    def test_install_git_http_credentials_variable_direct_modification(self):
        user_home_dirname = tempfile.mkdtemp()

        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.GitHttpCredentials.read_git_credentials(
            user_home_dirname=user_home_dirname
        )
        self.assertEqual([], other_config)
        self.assertEqual([], do_config)
        self.assertFalse(in_do_config)

        var = variables_types.GitHttpCredentials(
            drn="drn:1",
            name="var-1",
            repository="github.com/org/repo",
            username="john-doe",
            password="YYY",
        )
        var.install(None, user_home_dirname)
        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.GitHttpCredentials.read_git_credentials(
            variable=var, user_home_dirname=user_home_dirname
        )
        self.assertTrue(in_do_config)

        prev_var = copy.deepcopy(var)
        prev_var.username = "jane-doe"
        var_2 = copy.deepcopy(var)
        var_2.password = "XXX"
        with self.assertWarnsRegex(DeepOriginWarning, "modified directly"):
            var_2.install(prev_var, user_home_dirname)
        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.GitHttpCredentials.read_git_credentials(
            variable=var_2, user_home_dirname=user_home_dirname
        )
        self.assertFalse(in_do_config)

        prev_var = copy.deepcopy(var)
        var_2 = copy.deepcopy(var)
        var_2.password = "XXX"
        var_2.install(prev_var, user_home_dirname)
        (
            other_config,
            do_config,
            in_do_config,
        ) = variables_types.GitHttpCredentials.read_git_credentials(
            variable=var_2, user_home_dirname=user_home_dirname
        )
        self.assertTrue(in_do_config)

        shutil.rmtree(user_home_dirname)

    def test_expand_user(self):
        user_home_dirname = "/home/user"
        self.assertEqual(user_home_dirname, utils.expand_user("~", user_home_dirname))
        self.assertEqual(
            os.path.join(user_home_dirname, "subdir"),
            utils.expand_user(os.path.join("~", "subdir"), user_home_dirname),
        )
        self.assertEqual("subdir", utils.expand_user("subdir", user_home_dirname))
        self.assertEqual("subdir~", utils.expand_user("subdir~", user_home_dirname))
        self.assertEqual("~subdir", utils.expand_user("~subdir", user_home_dirname))

        user_home_dirname = str(pathlib.Path().home())
        self.assertEqual(
            os.path.expanduser("~"), utils.expand_user("~", user_home_dirname)
        )
        self.assertEqual(
            os.path.expanduser(os.path.join("~", "subdir")),
            utils.expand_user(os.path.join("~", "subdir"), user_home_dirname),
        )
        self.assertEqual(
            os.path.expanduser("subdir"),
            utils.expand_user("subdir", user_home_dirname),
        )
        self.assertEqual(
            os.path.expanduser("subdir~"),
            utils.expand_user("subdir~", user_home_dirname),
        )
        self.assertEqual(
            os.path.expanduser("~subdir"),
            utils.expand_user("~subdir", user_home_dirname),
        )

    def test_install_variables(self):
        user_home_dirname = tempfile.mkdtemp()

        auth.get_tokens.cache_clear()
        # os.remove(self.api_tokens_filename)
        api_responses = [
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "device_code": "xxx",
                    "user_code": "xxx",
                    "verification_uri_complete": "xxx",
                    "interval": 5,
                },
            ),
            unittest.mock.MagicMock(
                status_code=200,
                json=lambda: {
                    "access_token": "xxx",
                    "refresh_token": "xxx",
                },
            ),
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "data": {
                        "listBenchSecrets": [
                            {
                                "drn": "drn:1",
                                "type": "Non-secret environment variables",
                                "name": "abc",
                                "key": "def",
                                "value": "ghi",
                            },
                            {
                                "drn": "drn:2",
                                "type": "Secret file",
                                "name": "jkl",
                                "filename": "~/mno",
                                "value": "pqr",
                            },
                        ],
                    },
                },
            ),
        ]

        with unittest.mock.patch("requests.post", side_effect=api_responses):
            modification = variables.install_variables(
                user_home_dirname=user_home_dirname
            )
        self.assertEqual(
            set(["EnvironmentVariable", "SecretFile"]),
            set(
                [
                    m["type"]
                    for m in modification.values()
                    if m["status"] == variables.VariableStatus.added
                ]
            ),
        )

        auth.get_tokens.cache_clear()
        if os.path.isfile(self.api_tokens_filename):
            os.remove(self.api_tokens_filename)
        with unittest.mock.patch("requests.post", side_effect=api_responses):
            modification = variables.install_variables(
                user_home_dirname=user_home_dirname
            )
        self.assertEqual(
            set([variables.VariableStatus.unmodified]),
            set(
                [
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "EnvironmentVariable"
                ]
            ),
        )
        self.assertEqual(
            set([variables.VariableStatus.unmodified]),
            set(
                [
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "SecretFile"
                ]
            ),
        )

        auth.get_tokens.cache_clear()
        if os.path.isfile(self.api_tokens_filename):
            os.remove(self.api_tokens_filename)
        api_responses[-1] = unittest.mock.MagicMock(
            raise_for_status=lambda: None,
            json=lambda: {
                "data": {
                    "listBenchSecrets": [
                        {
                            "drn": "drn:1",
                            "type": "Non-secret environment variables",
                            "name": "abc",
                            "key": "def",
                            "value": "xyz",
                        },
                        {
                            "drn": "drn:2",
                            "type": "Secret file",
                            "name": "jkl",
                            "filename": "~/mno",
                            "value": "uvw",
                        },
                    ],
                },
            },
        )
        with unittest.mock.patch("requests.post", side_effect=api_responses):
            modification = variables.install_variables(
                user_home_dirname=user_home_dirname
            )
        self.assertEqual(
            set([variables.VariableStatus.modified]),
            set(
                [
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "EnvironmentVariable"
                ]
            ),
        )
        self.assertEqual(
            set([variables.VariableStatus.modified]),
            set(
                [
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "SecretFile"
                ]
            ),
        )

        auth.get_tokens.cache_clear()
        if os.path.isfile(self.api_tokens_filename):
            os.remove(self.api_tokens_filename)

        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
            with cli.App(argv=["variables", "install"]) as app:
                with unittest.mock.patch("requests.post", side_effect=api_responses):
                    app.run()
                stdout = stdout_capture.getvalue().strip()

        self.assertIn("No variables were modified", stdout)

        auth.get_tokens.cache_clear()
        if os.path.isfile(self.api_tokens_filename):
            os.remove(self.api_tokens_filename)
        api_responses[-1] = unittest.mock.MagicMock(
            raise_for_status=lambda: None,
            json=lambda: {
                "data": {
                    "listBenchSecrets": [
                        {
                            "drn": "drn:1",
                            "type": "Non-secret environment variables",
                            "name": "abc",
                            "key": "def",
                            "value": "xyz2",
                        },
                        {
                            "drn": "drn:2",
                            "type": "Secret file",
                            "name": "jkl",
                            "filename": "~/mno",
                            "value": "uvw2",
                        },
                    ],
                },
            },
        )

        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
            with cli.App(argv=["variables", "install"]) as app:
                with unittest.mock.patch("requests.post", side_effect=api_responses):
                    app.run()
                stdout = stdout_capture.getvalue().strip()

        self.assertIn("2 variables were modified", stdout)

        auth.get_tokens.cache_clear()
        if os.path.isfile(self.api_tokens_filename):
            os.remove(self.api_tokens_filename)
        api_responses[-1] = unittest.mock.MagicMock(
            raise_for_status=lambda: None,
            json=lambda: {
                "data": {
                    "listBenchSecrets": [
                        {
                            "drn": "drn:1",
                            "type": "Non-secret environment variables",
                            "name": "abc",
                            "key": "def",
                            "value": "xyz2-3",
                        },
                        {
                            "drn": "drn:2",
                            "type": "Secret file",
                            "name": "jkl",
                            "filename": "~/mno",
                            "value": "uvw2",
                        },
                    ],
                },
            },
        )

        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
            with cli.App(argv=["variables", "install"]) as app:
                with unittest.mock.patch("requests.post", side_effect=api_responses):
                    app.run()
                stdout = stdout_capture.getvalue().strip()

        self.assertIn("1 variables were modified", stdout)
        self.assertIn("1 variables were unmodified", stdout)

        auth.get_tokens.cache_clear()
        if os.path.isfile(self.api_tokens_filename):
            os.remove(self.api_tokens_filename)
        api_responses[-1] = unittest.mock.MagicMock(
            raise_for_status=lambda: None,
            json=lambda: {
                "data": {
                    "listBenchSecrets": [
                        {
                            "drn": "drn:2",
                            "type": "Secret file",
                            "name": "jkl",
                            "filename": "~/mno",
                            "value": "uvw2",
                        },
                        {
                            "drn": "drn:3",
                            "type": "Non-secret environment variables",
                            "name": "abc",
                            "key": "def",
                            "value": "xyz2-3",
                        },
                    ],
                },
            },
        )

        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
            with cli.App(argv=["variables", "install"]) as app:
                with unittest.mock.patch("requests.post", side_effect=api_responses):
                    app.run()
                stdout = stdout_capture.getvalue().strip()

        self.assertIn("No variables were modified", stdout)
        self.assertIn("1 variables were added", stdout)
        self.assertIn("1 variables were deleted", stdout)
        self.assertIn("1 variables were unmodified", stdout)

        with self.assertRaisesRegex(DeepOriginException, "are not valid"):
            with cli.App(argv=["variables", "install", "--type", "not-defined"]) as app:
                app.run()

        feature_flags.get_value().variables = False
        with cli.App(argv=["variables", "install"]) as app:
            with self.assertWarnsRegex(
                feature_flags.FeatureNotAvailableWarning, "not yet available"
            ):
                app.run()

        feature_flags.get_value().variables = True
        with cli.App(argv=["variables", "uninstall"]) as app:
            app.run()

        feature_flags.get_value().variables = False
        with cli.App(argv=["variables", "uninstall"]) as app:
            with self.assertWarnsRegex(
                feature_flags.FeatureNotAvailableWarning, "not yet available"
            ):
                app.run()

    def test_install_variables_2(self):
        user_home_dirname = tempfile.mkdtemp()
        env_filename = variables_types.EnvironmentVariable.get_env_filename(
            user_home_dirname=user_home_dirname
        )

        auth.get_tokens.cache_clear()
        # os.remove(self.api_tokens_filename)
        api_responses = [
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "device_code": "xxx",
                    "user_code": "xxx",
                    "verification_uri_complete": "xxx",
                    "interval": 5,
                },
            ),
            unittest.mock.MagicMock(
                status_code=200,
                json=lambda: {
                    "access_token": "xxx",
                    "refresh_token": "xxx",
                },
            ),
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "data": {
                        "listBenchSecrets": [
                            {
                                "drn": "drn:1",
                                "type": "Non-secret environment variables",
                                "name": "abc",
                                "key": "def",
                                "value": "ghi",
                            },
                            {
                                "drn": "drn:2",
                                "type": "Secret environment variables",
                                "name": "abc2",
                                "key": "def2",
                                "value": "ghi2",
                            },
                            {
                                "drn": "drn:3",
                                "type": "Secret file",
                                "name": "jkl",
                                "filename": "~/mno",
                                "value": "pqr",
                            },
                        ],
                    },
                },
            ),
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "data": {
                        "listBenchSecrets": [
                            {
                                "drn": "drn:2",
                                "type": "Secret environment variables",
                                "name": "abc2",
                                "key": "def2",
                                "value": "ghi3",
                            },
                        ],
                    },
                },
            ),
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "data": {
                        "listBenchSecrets": [
                            {
                                "drn": "drn:2",
                                "type": "Secret environment variables",
                                "name": "abc2",
                                "key": "def3",
                                "value": "ghi4",
                            },
                        ],
                    },
                },
            ),
        ]

        with unittest.mock.patch("requests.post", side_effect=api_responses):
            # 1
            modification = variables.install_variables(
                user_home_dirname=user_home_dirname
            )
            self.assertEqual(
                set([variables.VariableStatus.added]),
                set(
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "EnvironmentVariable"
                ),
            )
            self.assertEqual(
                set([variables.VariableStatus.added]),
                set(
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "SecretEnvironmentVariable"
                ),
            )
            self.assertEqual(
                set([variables.VariableStatus.added]),
                set(
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "SecretFile"
                ),
            )

            env = variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            )
            self.assertEqual({"def": "ghi", "def2": "ghi2"}, env)

            self.assertTrue(os.path.isfile(os.path.join(user_home_dirname, "mno")))

            # 2
            modification = variables.install_variables(
                user_home_dirname=user_home_dirname
            )
            self.assertEqual(
                set([variables.VariableStatus.deleted]),
                set(
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "EnvironmentVariable"
                ),
            )
            self.assertEqual(
                set([variables.VariableStatus.modified]),
                set(
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "SecretEnvironmentVariable"
                ),
            )
            self.assertEqual(
                set([variables.VariableStatus.deleted]),
                set(
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "SecretFile"
                ),
            )

            env = variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            )
            self.assertEqual({"def2": "ghi3"}, env)

            self.assertFalse(os.path.isfile(os.path.join(user_home_dirname, "mno")))

            # 3
            modification = variables.install_variables(
                user_home_dirname=user_home_dirname
            )
            self.assertEqual(
                set([]),
                set(
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "EnvironmentVariable"
                ),
            )
            self.assertEqual(
                set([variables.VariableStatus.modified]),
                set(
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "SecretEnvironmentVariable"
                ),
            )
            self.assertEqual(
                set([]),
                set(
                    m["status"]
                    for m in modification.values()
                    if m["type"] == "SecretFile"
                ),
            )

            env = variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                env_filename
            )
            self.assertEqual({"def3": "ghi4"}, env)

            self.assertFalse(os.path.isfile(os.path.join(user_home_dirname, "mno")))

    def test_install_invalid_variables(self):
        auth.get_tokens.cache_clear()
        # os.remove(self.api_tokens_filename)
        api_responses = [
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "device_code": "xxx",
                    "user_code": "xxx",
                    "verification_uri_complete": "xxx",
                    "interval": 5,
                },
            ),
            unittest.mock.MagicMock(
                status_code=200,
                json=lambda: {
                    "access_token": "xxx",
                    "refresh_token": "xxx",
                },
            ),
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "data": {
                        "listBenchSecrets": [
                            {
                                "drn": "drn:1",
                                "type": "Non-secret environment variables",
                                "name": "abc",
                                "key": "",
                                "value": "ghi",
                            },
                            {
                                "drn": "drn:2",
                                "type": "Secret file",
                                "name": "jkl",
                                "filename": "mno",
                                "value": "pqr",
                            },
                        ],
                    },
                },
            ),
        ]

        with self.assertRaisesRegex(DeepOriginException, "2 variables are not valid."):
            with cli.App(argv=["variables", "install"]) as app:
                with unittest.mock.patch("requests.post", side_effect=api_responses):
                    app.run()

    def test_enable_disable_variable_auto_updating(self):
        cron_job_id = variables_core.get_auto_install_variables_cronjob_id()

        responses = [
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "device_code": "xxx",
                    "user_code": "xxx",
                    "verification_uri_complete": "xxx",
                    "interval": 0.001,
                },
            ),
            unittest.mock.MagicMock(
                status_code=403,
                json=lambda: {
                    "error": "authorization_pending",
                },
            ),
            unittest.mock.MagicMock(
                status_code=200,
                json=lambda: {
                    "access_token": "xxx",
                    "refresh_token": "xxx",
                },
            ),
        ]
        with unittest.mock.patch("requests.post", side_effect=responses):
            variables.enable_variable_auto_updating(
                user=False, org=False, cli="deeporigin"
            )

        with crontab.CronTab(user=True) as cron_tab:
            self.assertIn(
                True,
                set(job.is_enabled() for job in cron_tab.find_comment(cron_job_id)),
            )

        with open(self.api_tokens_filename, "w"):
            pass
        variables.disable_variable_auto_updating()
        self.assertFalse(os.path.isfile(self.api_tokens_filename))

        with crontab.CronTab(user=True) as cron_tab:
            self.assertNotIn(
                True,
                set(job.is_enabled() for job in cron_tab.find_comment(cron_job_id)),
            )

        with cli.App(argv=["variables", "auto-install"]) as app:
            app.run()

        with cli.App(argv=["variables", "auto-install", "--disable"]) as app:
            app.run()

        with self.assertRaisesRegex(DeepOriginException, "are not valid"):
            with cli.App(
                argv=["variables", "auto-install", "--type", "not-defined"]
            ) as app:
                app.run()

        feature_flags.get_value().variables = False
        with cli.App(argv=["variables", "auto-install"]) as app:
            with self.assertWarnsRegex(
                feature_flags.FeatureNotAvailableWarning, "not yet available"
            ):
                app.run()

    def test_get_auto_install_variables_cronjob_id(self):
        self.assertTrue(
            variables_core.get_auto_install_variables_cronjob_id().startswith(
                "deep-origin-install-variables-"
            )
        )

    def test_get_remove_do_api_tokens(self):
        config.get_value.cache_clear()
        feature_flags.get_value.cache_clear()
        auth.get_tokens.cache_clear()

        _, self.api_tokens_filename = tempfile.mkstemp(suffix=".yml")
        os.remove(self.api_tokens_filename)

        env = {
            "DEEP_ORIGIN_ORGANIZATION_ID": "a123",
            "DEEP_ORIGIN_BENCH_ID": "b123",
            "DEEP_ORIGIN_ENV": "env123",
            "DEEP_ORIGIN_API_ENDPOINT": "c123",
            "DEEP_ORIGIN_API_TOKENS_FILENAME": self.api_tokens_filename,
            "DEEP_ORIGIN_VARIABLES_CACHE_FILENAME": self.variables_cache_filename,
            "DEEP_ORIGIN_AUTO_INSTALL_VARIABLES_FILENAME": self.auto_install_variables_filename,
            "DEEP_ORIGIN_USER_ID": "auth0|xxxxxxxxxxxxxxxxxxxxxxxx",
            "DEEP_ORIGIN_AUTH_DOMAIN": "xxx",
            "DEEP_ORIGIN_AUTH_DEVICE_CODE_ENDPOINT": "xxx",
            "DEEP_ORIGIN_AUTH_TOKEN_ENDPOINT": "xxx",
            "DEEP_ORIGIN_AUTH_AUDIENCE": "xxx",
            "DEEP_ORIGIN_LIST_BENCH_VARIABLES_QUERY_TEMPLATE": "xxx",
            "DEEP_ORIGIN_AUTH_GRANT_TYPE": "xxx",
            "DEEP_ORIGIN_AUTH_CLIENT_ID": "xxx",
            "DEEP_ORIGIN_AUTH_CLIENT_SECRET": "xxx",
        }
        with unittest.mock.patch.dict("os.environ", env):
            config.get_value()
            feature_flags.get_value()

        self.assertFalse(os.path.isfile(self.api_tokens_filename))

        responses = [
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "device_code": "xxx",
                    "user_code": "xxx",
                    "verification_uri_complete": "xxx",
                    "interval": 0.001,
                },
            ),
            unittest.mock.MagicMock(
                status_code=403,
                json=lambda: {
                    "error": "authorization_pending",
                },
            ),
            unittest.mock.MagicMock(
                status_code=200,
                json=lambda: {
                    "access_token": "xxx",
                    "refresh_token": "xxx",
                },
            ),
        ]
        with unittest.mock.patch("requests.post", side_effect=responses):
            auth.get_tokens()

        self.assertTrue(os.path.isfile(self.api_tokens_filename))

        responses = [
            unittest.mock.MagicMock(
                raise_for_status=lambda: None,
                json=lambda: {
                    "device_code": "xxx",
                    "user_code": "xxx",
                    "verification_uri_complete": "xxx",
                    "interval": 0.001,
                },
            ),
            unittest.mock.MagicMock(
                status_code=400,
                json=lambda: {
                    "error": "authorization_pending",
                },
            ),
        ]
        with unittest.mock.patch("requests.post", side_effect=responses):
            with self.assertRaisesRegex(
                DeepOriginException, "Sign in to the Deep Origin platform failed"
            ):
                auth.authenticate()

        auth.get_tokens.cache_clear()
        response = unittest.mock.MagicMock(
            raise_for_status=lambda: None,
            json=lambda: {
                "access_token": "xxx",
                "refresh_token": "xxx",
            },
        )
        with unittest.mock.patch("requests.post", return_value=response):
            auth.get_tokens()

        auth.remove_cached_tokens()
        self.assertFalse(os.path.isfile(self.api_tokens_filename))

    def test_validate_drn(self):
        variables_types.EnvironmentVariable(
            drn="drn:123", name="123", key="a", value="z"
        )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.EnvironmentVariable(
                drn="123", name="123", key="a", value="z"
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.EnvironmentVariable(drn=123, name="123", key="a", value="z")

    def test_validate_name(self):
        variables_types.EnvironmentVariable(
            drn="drn:123", name="123", key="a", value="z"
        )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.EnvironmentVariable(
                drn="drn:123", name="", key="a", value="z"
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.EnvironmentVariable(
                drn="drn:123", name=123, key="a", value="z"
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.EnvironmentVariable(
                drn="drn:123", name=None, key="a", value="z"
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.EnvironmentVariable(drn="drn:123", key="a", value="z")

    def test_validate_environment_variable(self):
        variables_types.EnvironmentVariable(
            drn="drn:123", name="123", key="a123_", value="z"
        )
        variables_types.EnvironmentVariable(
            drn="drn:123", name="123", key="a123_", value=""
        )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.EnvironmentVariable(
                drn="drn:123", name="123", key="123a_", value="z"
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.EnvironmentVariable(
                drn="drn:123", name="123", key="a123_", value=None
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.EnvironmentVariable(drn="drn:123", name="123", key="a123_")

        with self.assertRaises(pydantic.ValidationError):
            variables_types.EnvironmentVariable(
                drn="drn:123", name="123", key="a123_", value=123
            )

    def test_validate_file(self):
        variables_types.File(drn="drn:123", name="123", filename="/abc", value="z")
        variables_types.File(drn="drn:123", name="123", filename="~/abc", value="z")
        variables_types.File(drn="drn:123", name="123", filename="/abc", value="")

        with self.assertRaises(pydantic.ValidationError):
            variables_types.File(drn="drn:123", name="123", filename="abc", value="z")

        with self.assertRaises(pydantic.ValidationError):
            variables_types.File(drn="drn:123", name="123", filename=123, value="z")

        with self.assertRaises(pydantic.ValidationError):
            variables_types.File(drn="drn:123", name="123", filename=None, value="z")

        with self.assertRaises(pydantic.ValidationError):
            variables_types.File(drn="drn:123", name="123", filename="/abc", value=123)

        with self.assertRaises(pydantic.ValidationError):
            variables_types.File(drn="drn:123", name="123", filename="/abc", value=None)

        with self.assertRaises(pydantic.ValidationError):
            variables_types.File(drn="drn:123", name="123", filename="/abc")

    def test_validate_private_ssh_key(self):
        variables_types.PrivateSshKey(
            drn="drn:123",
            name="123",
            filename="/.ssh/id_rsa",
            value=self.TEST_SSH_PRIVATE_KEY,
        )
        variables_types.PrivateSshKey(
            drn="drn:123",
            name="123",
            filename="~/.ssh/id_rsa",
            value=self.TEST_SSH_PRIVATE_KEY,
        )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.PrivateSshKey(
                drn="drn:123",
                name="123",
                filename="/.ssh/id_rsa.pub",
                value=self.TEST_SSH_PUBLIC_KEY,
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.PrivateSshKey(
                drn="drn:123", name="123", filename="/.ssh/id_rsa", value=""
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.PrivateSshKey(
                drn="drn:123", name="123", filename="/.ssh/id_rsa", value=123
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.PrivateSshKey(
                drn="drn:123", name="123", filename="/.ssh/id_rsa", value=None
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.PrivateSshKey(
                drn="drn:123", name="123", filename="/.ssh/id_rsa"
            )

    def test_validate_private_gpg_key(self):
        variables_types.PrivateGpgKey(
            drn="drn:123",
            name="123",
            value=self.TEST_PRIVATE_GPG_KEY,
        )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.PrivateGpgKey(drn="drn:123", name="123", value="123")

        with self.assertRaises(pydantic.ValidationError):
            variables_types.PrivateGpgKey(drn="drn:123", name="123", value=123)

        with self.assertRaises(pydantic.ValidationError):
            variables_types.PrivateGpgKey(drn="drn:123", name="123", value=None)

        with self.assertRaises(pydantic.ValidationError):
            variables_types.PrivateGpgKey(drn="drn:123", name="123")

    def test_validate_aws_profile(self):
        variables_types.AwsProfile(
            drn="drn:123",
            name="123",
            profile="user",
            access_key_id="X" * 20,
            secret_access_key="Y" * 40,
            region="us-west-1",
        )
        variables_types.AwsProfile(
            drn="drn:123",
            name="123",
            profile="user",
            access_key_id="X" * 20,
            secret_access_key="y" * 40,
            region="us-west-1",
        )
        variables_types.AwsProfile(
            drn="drn:123",
            name="123",
            profile="user",
            access_key_id="X" * 20,
            secret_access_key="Y" * 40,
            region=None,
        )
        variables_types.AwsProfile(
            drn="drn:123",
            name="123",
            profile="user",
            access_key_id="X" * 20,
            secret_access_key="Y" * 40,
        )

        # profile
        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="",
                access_key_id="X" * 20,
                secret_access_key="Y" * 40,
                region="us-west-1",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile=None,
                access_key_id="X" * 20,
                secret_access_key="Y" * 40,
                region="us-west-1",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                access_key_id="X" * 20,
                secret_access_key="Y" * 40,
                region="us-west-1",
            )

        # access key ID
        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                access_key_id="x" * 20,
                secret_access_key="Y" * 40,
                region="us-west-1",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                access_key_id="X" * 19,
                secret_access_key="Y" * 40,
                region="us-west-1",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                access_key_id="",
                secret_access_key="Y" * 40,
                region="us-west-1",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                access_key_id=None,
                secret_access_key="Y" * 40,
                region="us-west-1",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                secret_access_key="Y" * 40,
                region="us-west-1",
            )

        # secret access key
        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                access_key_id="X" * 20,
                secret_access_key="Y" * 39,
                region="us-west-1",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                access_key_id="X" * 20,
                secret_access_key="",
                region="us-west-1",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                access_key_id="X" * 20,
                secret_access_key=None,
                region="us-west-1",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                access_key_id="X" * 20,
                region="us-west-1",
            )

        # region
        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                access_key_id="X" * 20,
                secret_access_key="Y" * 40,
                region="US-WEST-1",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                access_key_id="X" * 20,
                secret_access_key="Y" * 40,
                region="us",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                access_key_id="X" * 20,
                secret_access_key="Y" * 40,
                region=1,
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.AwsProfile(
                drn="drn:123",
                name="123",
                profile="user",
                access_key_id="X" * 20,
                secret_access_key="Y" * 40,
                region="",
            )

    def test_validate_git_http_credentials(self):
        variables_types.GitHttpCredentials(
            drn="drn:123",
            name="123",
            repository="github.com/org/repo",
            username="user",
            password="pass",
        )
        variables_types.GitHttpCredentials(
            drn="drn:123",
            name="123",
            repository="https://github.com/org/repo",
            username="user",
            password="pass",
        )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GitHttpCredentials(
                drn="drn:123",
                name="123",
                repository="https://https://github.com/org/repo",
                username="user",
                password="pass",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GitHttpCredentials(
                drn="drn:123",
                name="123",
                repository="",
                username="user",
                password="pass",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GitHttpCredentials(
                drn="drn:123",
                name="123",
                repository=123,
                username="user",
                password="pass",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GitHttpCredentials(
                drn="drn:123",
                name="123",
                repository=None,
                username="user",
                password="pass",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GitHttpCredentials(
                drn="drn:123", name="123", username="user", password="pass"
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GitHttpCredentials(
                drn="drn:123",
                name="123",
                repository="github.com/org/repo",
                username="",
                password="pass",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GitHttpCredentials(
                drn="drn:123",
                name="123",
                repository="github.com/org/repo",
                username=123,
                password="pass",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GitHttpCredentials(
                drn="drn:123",
                name="123",
                repository="github.com/org/repo",
                username=None,
                password="pass",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GitHttpCredentials(
                drn="drn:123",
                name="123",
                repository="github.com/org/repo",
                password="pass",
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GitHttpCredentials(
                drn="drn:123",
                name="123",
                repository="github.com/org/repo",
                username="user",
                password=123,
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GitHttpCredentials(
                drn="drn:123",
                name="123",
                repository="github.com/org/repo",
                username="user",
                password=None,
            )

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GitHttpCredentials(
                drn="drn:123",
                name="123",
                repository="github.com/org/repo",
                username="user",
            )

    def test_validate_open_ai_api_key(self):
        variables_types.OpenAiApiKey(drn="drn:123", name="123", value="z")

        with self.assertRaises(pydantic.ValidationError):
            variables_types.OpenAiApiKey(drn="drn:123", name="123", value="")

        with self.assertRaises(pydantic.ValidationError):
            variables_types.OpenAiApiKey(drn="drn:123", name="123", value=123)

        with self.assertRaises(pydantic.ValidationError):
            variables_types.OpenAiApiKey(drn="drn:123", name="123", value=None)

        with self.assertRaises(pydantic.ValidationError):
            variables_types.OpenAiApiKey(drn="drn:123", name="123")

    def test_validate_anthropic_api_key(self):
        variables_types.OpenAiApiKey(drn="drn:123", name="123", value="z")

        with self.assertRaises(pydantic.ValidationError):
            variables_types.OpenAiApiKey(drn="drn:123", name="123", value="")

        with self.assertRaises(pydantic.ValidationError):
            variables_types.OpenAiApiKey(drn="drn:123", name="123", value=123)

        with self.assertRaises(pydantic.ValidationError):
            variables_types.OpenAiApiKey(drn="drn:123", name="123", value=None)

        with self.assertRaises(pydantic.ValidationError):
            variables_types.OpenAiApiKey(drn="drn:123", name="123")

    def test_validate_file_value(self):
        variables_types.GurobiLicenseFile(drn="drn:123", name="123", value="z")
        variables_types.GurobiLicenseFile(drn="drn:123", name="123", value="")

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GurobiLicenseFile(drn="drn:123", name="123", value=123)

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GurobiLicenseFile(drn="drn:123", name="123", value=None)

        with self.assertRaises(pydantic.ValidationError):
            variables_types.GurobiLicenseFile(drn="drn:123", name="123")

    def test_validate_express_license_file(self):
        variables_types.XpressLicenseFile(drn="drn:123", name="123", value="z")
        variables_types.XpressLicenseFile(drn="drn:123", name="123", value="")

        with self.assertRaises(pydantic.ValidationError):
            variables_types.XpressLicenseFile(drn="drn:123", name="123", value=123)

        with self.assertRaises(pydantic.ValidationError):
            variables_types.XpressLicenseFile(drn="drn:123", name="123", value=None)

        with self.assertRaises(pydantic.ValidationError):
            variables_types.XpressLicenseFile(drn="drn:123", name="123")

    def test_secret_env_var_value_key(self):
        var = variables_types.OpenAiApiKey(drn="drn:123", name="123", value="z")
        self.assertEqual(var.KEY, "OPENAI_API_KEY")

        var = variables_types.AnthropicApiKey(drn="drn:123", name="123", value="z")
        self.assertEqual(var.KEY, "ANTHROPIC_API_KEY")

    def test_secret_file_value_filename(self):
        var = variables_types.GurobiLicenseFile(drn="drn:123", name="123", value="z")
        self.assertEqual(var.FILENAME, "~/gurobi.lic")

        var = variables_types.MosekLicenseFile(drn="drn:123", name="123", value="z")
        self.assertEqual(var.FILENAME, "~/mosek/mosek.lic")

        var = variables_types.GurobiLicenseFile(drn="drn:123", name="123", value="z")
        self.assertEqual(var.FILENAME, "~/gurobi.lic")

    def test_xpress_license_file_filename_key(self):
        var = variables_types.XpressLicenseFile(drn="drn:123", name="123", value="z")
        self.assertEqual(var.FILENAME, "~/xpauth.xpr")

        var = variables_types.XpressLicenseFile(drn="drn:123", name="123", value="z")
        self.assertEqual(var.KEY, "XPAUTH_PATH")

    def test_install_secret_env_var_value(self):
        user_home_dirname = tempfile.mkdtemp()
        user_env_filename = variables_types.EnvironmentVariable.get_env_filename(
            is_user_variable=True,
            user_home_dirname=user_home_dirname,
        )
        system_env_filename = variables_types.EnvironmentVariable.get_env_filename(
            is_user_variable=False,
            user_home_dirname=user_home_dirname,
        )

        expected_system_env = {}
        self.assertEqual(
            {},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                user_env_filename
            ),
        )
        self.assertEqual(
            expected_system_env,
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                system_env_filename
            ),
        )

        var_open_ai = variables_types.OpenAiApiKey(drn="drn:123", name="123", value="x")
        var_open_ai.install(None, user_home_dirname)
        expected_system_env[var_open_ai.KEY] = var_open_ai.value
        self.assertEqual(
            {},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                user_env_filename
            ),
        )
        self.assertEqual(
            expected_system_env,
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                system_env_filename
            ),
        )

        var_anthropic = variables_types.AnthropicApiKey(
            drn="drn:123", name="123", value="y"
        )
        var_anthropic.install(None, user_home_dirname)
        expected_system_env[var_anthropic.KEY] = var_anthropic.value
        self.assertEqual(
            {},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                user_env_filename
            ),
        )
        self.assertEqual(
            expected_system_env,
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                system_env_filename
            ),
        )

        var_xpress = variables_types.XpressLicenseFile(
            drn="drn:123", name="123", value="z"
        )
        var_xpress.install(None, user_home_dirname)
        expected_system_env[var_xpress.KEY] = utils.expand_user(
            var_xpress.FILENAME, user_home_dirname=user_home_dirname
        )
        self.assertEqual(
            {},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                user_env_filename
            ),
        )
        self.assertEqual(
            expected_system_env,
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                system_env_filename
            ),
        )

        var_open_ai.uninstall(user_home_dirname=user_home_dirname)
        expected_system_env.pop(var_open_ai.KEY)
        self.assertEqual(
            {},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                user_env_filename
            ),
        )
        self.assertEqual(
            expected_system_env,
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                system_env_filename
            ),
        )

        var_anthropic.uninstall(user_home_dirname=user_home_dirname)
        expected_system_env.pop(var_anthropic.KEY)
        self.assertEqual(
            {},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                user_env_filename
            ),
        )
        self.assertEqual(
            expected_system_env,
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                system_env_filename
            ),
        )

        var_xpress.uninstall(user_home_dirname=user_home_dirname)
        expected_system_env.pop(var_xpress.KEY)
        self.assertEqual(
            {},
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                user_env_filename
            ),
        )
        self.assertEqual(
            expected_system_env,
            variables_types.EnvironmentVariable.read_deep_origin_environment_variables(
                system_env_filename
            ),
        )

    def test_install_secret_file_value(self):
        user_home_dirname = tempfile.mkdtemp()

        expected_filenames = []
        self.assertEqual(
            sorted(expected_filenames), sorted(os.listdir(user_home_dirname))
        )

        var_gurobi = variables_types.GurobiLicenseFile(
            drn="drn:1",
            name="var-1",
            value="x",
        )
        var_gurobi.install(None, user_home_dirname)
        expected_filenames.append(var_gurobi.FILENAME.replace("~/", ""))
        self.assertEqual(
            sorted(expected_filenames), sorted(os.listdir(user_home_dirname))
        )
        with open(
            utils.expand_user(var_gurobi.FILENAME, user_home_dirname), "r"
        ) as file:
            self.assertEqual(var_gurobi.value, file.read())

        var_mosek = variables_types.MosekLicenseFile(
            drn="drn:1",
            name="var-1",
            value="y",
        )
        var_mosek.install(None, user_home_dirname)
        expected_filenames.append(os.path.dirname(var_mosek.FILENAME.replace("~/", "")))
        self.assertEqual(
            sorted(expected_filenames), sorted(os.listdir(user_home_dirname))
        )
        with open(
            utils.expand_user(var_mosek.FILENAME, user_home_dirname), "r"
        ) as file:
            self.assertEqual(var_mosek.value, file.read())

        var_xpress = variables_types.XpressLicenseFile(
            drn="drn:1",
            name="var-1",
            value="z",
        )
        var_xpress.install(None, user_home_dirname)
        expected_filenames.append(var_xpress.FILENAME.replace("~/", ""))
        expected_filenames.append(".bashrc")
        expected_filenames.append(".deeporigin")
        self.assertEqual(
            sorted(expected_filenames), sorted(os.listdir(user_home_dirname))
        )
        with open(
            utils.expand_user(var_xpress.FILENAME, user_home_dirname), "r"
        ) as file:
            self.assertEqual(var_xpress.value, file.read())

        var_gurobi.uninstall(user_home_dirname=user_home_dirname)
        expected_filenames.remove(var_gurobi.FILENAME.replace("~/", ""))
        self.assertEqual(
            sorted(expected_filenames), sorted(os.listdir(user_home_dirname))
        )
        self.assertFalse(
            os.path.isfile(utils.expand_user(var_gurobi.FILENAME, user_home_dirname))
        )

        var_mosek.uninstall(user_home_dirname=user_home_dirname)
        self.assertEqual(
            sorted(expected_filenames), sorted(os.listdir(user_home_dirname))
        )
        self.assertFalse(
            os.path.isfile(utils.expand_user(var_mosek.FILENAME, user_home_dirname))
        )

        var_xpress.uninstall(user_home_dirname=user_home_dirname)
        expected_filenames.remove(var_xpress.FILENAME.replace("~/", ""))
        self.assertEqual(
            sorted(expected_filenames), sorted(os.listdir(user_home_dirname))
        )
        self.assertFalse(
            os.path.isfile(utils.expand_user(var_xpress.FILENAME, user_home_dirname))
        )
