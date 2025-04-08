from contextlib import redirect_stderr, redirect_stdout
import io
import json
import unittest
import unittest.mock

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
            "DEEP_ORIGIN_COMPUTE_CLUSTER": json.dumps(
                {
                    "id": "us-west-2.aws.bench.deeporigin.io",
                    "properties": {
                        "provider": "aws",
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
            workstation_id="cute-bear-123",
            user_id="john-doe",
            org_id="acme-bio",
            compute_cluster=context.ComputeClusterContext(
                id="us-west-2.aws.bench.deeporigin.io",
                properties={"provider": "aws", "region": "us-west-2"},
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

        expected_str = "\n".join(
            [
                "Workstation ID: cute-bear-123",
                "Compute cluster:",
                "  ID: us-west-2.aws.bench.deeporigin.io",
                "  Properties:",
                "    provider: aws",
                "    region: us-west-2",
                "User ID: john-doe",
                "Organization ID: acme-bio",
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
        self.assertEqual(
            expected_str, value.__str__(compute_cluster=True, hardware=True)
        )

    def test_get_project_paths_nested_null(self):
        context.get_value.cache_clear()
        env = {
            "DEEP_ORIGIN_COMPUTE_CLUSTER": json.dumps(
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
            workstation_id=None,
            user_id=None,
            org_id=None,
            compute_cluster=context.ComputeClusterContext(
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

        expected_str = "\n".join(
            [
                "Workstation ID: None",
                "Compute cluster:",
                "  ID: None",
                "  Properties: None",
                "User ID: None",
                "Organization ID: None",
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
        self.assertEqual(
            expected_str, value.__str__(compute_cluster=True, hardware=True)
        )

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
            workstation_id=None,
            user_id=None,
            org_id=None,
            compute_cluster=None,
            hardware=context.HardwareContext(
                cpu=None,
                gpu=None,
            ),
            env=None,
            debug=False,
        )
        self.assertEqual(expected_value, value)

        expected_str = "\n".join(
            [
                "Workstation ID: None",
                "Compute cluster: None",
                "User ID: None",
                "Organization ID: None",
                "Hardware:",
                "  Processing: N/A",
                "  Accelerated processing: None",
                "Environment: None",
                "Debug: False",
            ]
        )
        self.assertEqual(
            expected_str, value.__str__(compute_cluster=True, hardware=True)
        )

    def test_get_project_paths_null_2(self):
        context.get_value.cache_clear()
        env = {}
        with unittest.mock.patch.dict("os.environ", env):
            value = context.get_value()
        expected_value = context.Context(
            workstation_id=None,
            user_id=None,
            org_id=None,
            compute_cluster=None,
            hardware=None,
            env=None,
            debug=False,
        )
        self.assertEqual(expected_value, value)

        expected_str_with_cluster_hardware = "\n".join(
            [
                "Workstation ID: None",
                "Compute cluster: None",
                "User ID: None",
                "Organization ID: None",
                "Hardware: None",
                "Environment: None",
                "Debug: False",
            ]
        )
        self.assertEqual(
            expected_str_with_cluster_hardware,
            value.__str__(compute_cluster=True, hardware=True),
        )

        expected_str_without_cluster_hardware = "\n".join(
            [
                "Workstation ID: None",
                "User ID: None",
                "Organization ID: None",
                "Environment: None",
                "Debug: False",
            ]
        )
        self.assertEqual(
            expected_str_without_cluster_hardware,
            value.__str__(compute_cluster=False, hardware=False),
        )

        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()

        with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
            with cli.App(argv=["context"]) as app:
                app.run()

            stdout = stdout_capture.getvalue().strip()

        expected_stdout = expected_str_without_cluster_hardware
        self.assertEqual(expected_stdout, stdout)
