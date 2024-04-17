import io
import json
import unittest
import unittest.mock
from contextlib import redirect_stderr, redirect_stdout

from deeporigin import cli, context


class TestCase(unittest.TestCase):
    def setUp(self):
        context.get_value.cache_clear()

    def tearDown(self):
        context.get_value.cache_clear()

    def test_get_project_paths(self):
        context.get_value.cache_clear()
        env = {
            "DEEP_ORIGIN_BENCH_ID": "cute-bear-123",
            "DEEP_ORIGIN_USER_ID": "john-doe",
            "DEEP_ORIGIN_ORGANIZATION_ID": "acme-bio",
            "DEEP_ORIGIN_HOST": json.dumps(
                {
                    "id": "aws",
                    "properties": {
                        "region": "us-west-2",
                    },
                }
            ),
            "DEEP_ORIGIN_HARDWARE": json.dumps(
                {
                    "cpu": {
                        "quantityVcpu": 8,
                        "architecture": "amd64",
                        "memoryGB": 32,
                    },
                    "gpu": {
                        "quantityGpu": 1,
                        "architecture": {
                            "vendor": "nvidia",
                            "microarchitecture": "tesla",
                            "model": "t4",
                        },
                        "memoryGB": 4,
                    },
                }
            ),
            "DEEP_ORIGIN_ENV": "env123",
            "DEEP_ORIGIN_DEBUG": "true",
        }
        with unittest.mock.patch.dict("os.environ", env):
            value = context.get_value()
        expected_value = context.Context(
            bench_id="cute-bear-123",
            user_id="john-doe",
            org_id="acme-bio",
            host=context.HostContext(
                id="aws",
                properties={"region": "us-west-2"},
            ),
            hardware=context.HardwareContext(
                cpu=context.CpuContext(
                    quantity_vcpu=8,
                    architecture="amd64",
                    memory_gb=32,
                ),
                gpu=context.GpuContext(
                    quantity_gpu=1,
                    architecture=context.GpuArchitectureContext(
                        vendor="nvidia",
                        microarchitecture="tesla",
                        model="t4",
                    ),
                    memory_gb=4,
                ),
            ),
            env="env123",
            debug=True,
        )
        self.assertEqual(expected_value, value)

        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
            with cli.App(argv=["context"]) as app:
                app.run()
            stdout = stdout_capture.getvalue().strip()
        expected_stdout = "\n".join(
            [
                "Bench ID: cute-bear-123",
                "User ID: john-doe",
                "Organization ID: acme-bio",
                "Host:",
                "  ID: aws",
                "  Properties:",
                "    region: us-west-2",
                "Hardware:",
                "  Processing:",
                "    Quantity: 8.0 vCPU",
                "    Architecture: amd64",
                "    Memory: 32.0 GB",
                "  Accelerated processing:",
                "    Quantity: 1.0 GPU",
                "    Architecture:",
                "      Vendor: nvidia",
                "      Microarchitecture: tesla",
                "      Model: t4",
                "    Memory: 4.0 GB",
                "Environment: env123",
                "Debug: True",
            ]
        )
        self.assertEqual(expected_stdout, stdout)

    def test_get_project_paths_nested_null(self):
        context.get_value.cache_clear()
        env = {
            "DEEP_ORIGIN_HOST": json.dumps(
                {
                    "id": None,
                    "properties": None,
                }
            ),
            "DEEP_ORIGIN_HARDWARE": json.dumps(
                {
                    "cpu": {
                        "quantityVcpu": None,
                        "architecture": None,
                        "memoryGB": None,
                    },
                    "gpu": {
                        "quantityGpu": None,
                        "architecture": None,
                        "memoryGB": None,
                    },
                }
            ),
        }
        with unittest.mock.patch.dict("os.environ", env):
            value = context.get_value()
        expected_value = context.Context(
            bench_id=None,
            user_id=None,
            org_id=None,
            host=context.HostContext(
                id=None,
                properties=None,
            ),
            hardware=context.HardwareContext(
                cpu=context.CpuContext(
                    quantity_vcpu=None,
                    architecture=None,
                    memory_gb=None,
                ),
                gpu=context.GpuContext(
                    quantity_gpu=None,
                    architecture=None,
                    memory_gb=None,
                ),
            ),
            env=None,
            debug=False,
        )
        self.assertEqual(expected_value, value)

        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
            with cli.App(argv=["context"]) as app:
                app.run()

            stdout = stdout_capture.getvalue().strip()

        expected_stdout = "\n".join(
            [
                "Bench ID: None",
                "User ID: None",
                "Organization ID: None",
                "Host:",
                "  ID: None",
                "  Properties: None",
                "Hardware:",
                "  Processing:",
                "    Quantity: N/A",
                "    Architecture: N/A",
                "    Memory: N/A",
                "  Accelerated processing:",
                "    Quantity: None",
                "    Architecture: None",
                "    Memory: None",
                "Environment: None",
                "Debug: False",
            ]
        )
        self.assertEqual(expected_stdout, stdout)

    def test_get_project_paths_null_1(self):
        context.get_value.cache_clear()
        env = {
            "DEEP_ORIGIN_HARDWARE": json.dumps(
                {
                    "cpu": None,
                    "gpu": None,
                }
            ),
        }
        with unittest.mock.patch.dict("os.environ", env):
            value = context.get_value()
        expected_value = context.Context(
            bench_id=None,
            user_id=None,
            org_id=None,
            host=None,
            hardware=context.HardwareContext(
                cpu=None,
                gpu=None,
            ),
            env=None,
            debug=False,
        )
        self.assertEqual(expected_value, value)

        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
            with cli.App(argv=["context"]) as app:
                app.run()

            stdout = stdout_capture.getvalue().strip()
        expected_stdout = "\n".join(
            [
                "Bench ID: None",
                "User ID: None",
                "Organization ID: None",
                "Host: None",
                "Hardware:",
                "  Processing: N/A",
                "  Accelerated processing: None",
                "Environment: None",
                "Debug: False",
            ]
        )
        self.assertEqual(expected_stdout, stdout)

    def test_get_project_paths_null_2(self):
        context.get_value.cache_clear()
        env = {}
        with unittest.mock.patch.dict("os.environ", env):
            value = context.get_value()
        expected_value = context.Context(
            bench_id=None,
            user_id=None,
            org_id=None,
            host=None,
            hardware=None,
            env=None,
            debug=False,
        )
        self.assertEqual(expected_value, value)

        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
            with cli.App(argv=["context"]) as app:
                app.run()

            stdout = stdout_capture.getvalue().strip()

        expected_stdout = "\n".join(
            [
                "Bench ID: None",
                "User ID: None",
                "Organization ID: None",
                "Host: None",
                "Hardware: None",
                "Environment: None",
                "Debug: False",
            ]
        )
        self.assertEqual(expected_stdout, stdout)
