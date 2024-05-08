import dataclasses
import functools
import json
import os
import typing

__all__ = [
    "get_value",
    "Context",
    "HostContext",
    "HardwareContext",
    "CpuContext",
    "GpuContext",
    "GpuArchitectureContext",
]


@dataclasses.dataclass
class HostContext:
    """Host for a Deep Origin ComputeBench"""

    id: str
    properties: dict[str, typing.Any]


@dataclasses.dataclass
class CpuContext:
    """CPU for a Deep Origin ComputeBench"""

    quantity_vcpu: float
    architecture: str
    memory_gb: float


@dataclasses.dataclass
class GpuArchitectureContext:
    """GPU architecture for a Deep Origin ComputeBench"""

    vendor: str
    microarchitecture: str
    model: str


@dataclasses.dataclass
class GpuContext:
    """GPU for a Deep Origin ComputeBench"""

    quantity_gpu: int
    architecture: GpuArchitectureContext
    memory_gb: float


@dataclasses.dataclass
class HardwareContext:
    """Hardware for a Deep Origin ComputeBench"""

    cpu: CpuContext
    gpu: typing.Optional[GpuContext]


@dataclasses.dataclass
class Context:
    """Context for a Deep Origin ComputeBench"""

    bench_id: str
    user_id: str
    org_id: str
    host: HostContext
    hardware: HardwareContext
    env: str
    debug: bool


@functools.cache
def get_value():
    """returns a context from environment variables"""
    host_str = os.getenv("DEEP_ORIGIN_HOST", None)
    if host_str:
        host_dict = json.loads(host_str)
        host = HostContext(
            id=host_dict.get("id", None), properties=host_dict.get("properties", None)
        )
    else:
        host = None

    hardware_str = os.getenv("DEEP_ORIGIN_HARDWARE", None)
    if hardware_str:
        hardware_dict = json.loads(hardware_str)

        cpu_dict = hardware_dict.get("cpu", None)
        if cpu_dict:
            cpu = CpuContext(
                quantity_vcpu=cpu_dict.get("quantityVcpu", None),
                architecture=cpu_dict.get("architecture", None),
                memory_gb=cpu_dict.get("memoryGB", None),
            )
        else:
            cpu = None

        gpu_dict = hardware_dict.get("gpu", None)
        if gpu_dict:
            gpu_architecture_dict = gpu_dict.get("architecture", None)
            if gpu_architecture_dict:
                architecture = GpuArchitectureContext(
                    vendor=gpu_architecture_dict.get("vendor", None),
                    microarchitecture=gpu_architecture_dict.get(
                        "microarchitecture", None
                    ),
                    model=gpu_architecture_dict.get("model", None),
                )
            else:
                architecture = None

            gpu = GpuContext(
                quantity_gpu=gpu_dict.get("quantityGpu", None),
                architecture=architecture,
                memory_gb=gpu_dict.get("memoryGB", None),
            )
        else:
            gpu = None

        hardware = HardwareContext(
            cpu=cpu,
            gpu=gpu,
        )
    else:
        hardware = None

    return Context(
        bench_id=os.getenv("DEEP_ORIGIN_BENCH_ID", None),
        user_id=os.getenv("DEEP_ORIGIN_USER_ID", None),
        org_id=os.getenv("DEEP_ORIGIN_ORGANIZATION_ID", None),
        host=host,
        hardware=hardware,
        env=os.getenv("DEEP_ORIGIN_ENV", None),
        debug=os.getenv("DEEP_ORIGIN_DEBUG", "").lower() == "true",
    )
