"""this module provides functions used by the CLI to provide context within a Deep Origin workstation"""

import dataclasses
import functools
import json
import os
import typing

__all__ = [
    "get_value",
    "Context",
    "ComputeClusterContext",
    "HardwareContext",
    "CpuContext",
    "GpuContext",
    "GpuArchitectureContext",
]


@dataclasses.dataclass
class ComputeClusterContext:
    """Host compute cluster for a Deep Origin workstation"""

    id: str
    properties: dict[str, typing.Any]

    def __str__(self) -> str:
        str_val = []

        str_val.append(f"ID: {self.id}")
        if self.properties:
            str_val.append("Properties:")
            for key, value in self.properties.items():
                str_val.append(f"  {key}: {value}")
        else:
            str_val.append(f"Properties: {None}")

        return "\n".join(str_val)


@dataclasses.dataclass
class CpuContext:
    """CPU for a Deep Origin workstation"""

    quantity_vcpu: float
    architecture: str
    memory_gb: float

    def __str__(self) -> str:
        str_val = []

        if self.quantity_vcpu is not None:
            str_val.append(f"Quantity: {self.quantity_vcpu:.1f} vCPU")
        else:
            str_val.append("Quantity: N/A")
        str_val.append(f"Architecture: {self.architecture or 'N/A'}")
        if self.memory_gb is not None:
            str_val.append(f"Memory: {self.memory_gb:.1f} GB")
        else:
            str_val.append("Memory: N/A")

        return "\n".join(str_val)


@dataclasses.dataclass
class GpuArchitectureContext:
    """GPU architecture for a Deep Origin workstation"""

    vendor: str
    microarchitecture: str
    model: str

    def __str__(self) -> str:
        str_val = []

        str_val.append(f"Vendor: {self.vendor}")
        str_val.append(f"Microarchitecture: {self.microarchitecture}")
        str_val.append(f"Model: {self.model}")

        return "\n".join(str_val)


@dataclasses.dataclass
class GpuContext:
    """GPU for a Deep Origin"""

    quantity_gpu: int
    architecture: GpuArchitectureContext
    memory_gb: float

    def __str__(self) -> str:
        str_val = []

        if self.quantity_gpu is not None:
            str_val.append(f"Quantity: {self.quantity_gpu:.1f} GPU")
        else:
            str_val.append("Quantity: None")
        if self.architecture:
            str_val.append("Architecture:")
            str_val.append("  " + str(self.architecture).replace("\n", "\n  "))
        else:
            str_val.append(f"Architecture: {None}")
        if self.memory_gb is not None:
            str_val.append(f"Memory: {self.memory_gb:.1f} GB")
        else:
            str_val.append("Memory: None")

        return "\n".join(str_val)


@dataclasses.dataclass
class HardwareContext:
    """Hardware for a Deep Origin workstation"""

    cpu: CpuContext
    gpu: typing.Optional[GpuContext]

    def __str__(self) -> str:
        str_val = []

        # CPU
        if self.cpu:
            str_val.append("Processing:")
            str_val.append("  " + str(self.cpu).replace("\n", "\n  "))
        else:
            str_val.append("Processing: N/A")

        # GPU
        if self.gpu:
            str_val.append("Accelerated processing:")
            str_val.append("  " + str(self.gpu).replace("\n", "\n  "))

        else:
            str_val.append(f"Accelerated processing: {None}")

        return "\n".join(str_val)


@dataclasses.dataclass
class Context:
    """Context for a Deep Origin workstation"""

    workstation_id: str
    user_id: str
    org_id: str
    compute_cluster: ComputeClusterContext
    hardware: HardwareContext
    env: str
    debug: bool

    def __str__(
        self,
        workstation: bool = True,
        compute_cluster: bool = False,
        user: bool = True,
        org: bool = True,
        hardware: bool = False,
        env: bool = True,
        debug: bool = True,
    ) -> str:
        str_val = []

        # workstation id
        if workstation:
            str_val.append(f"Workstation ID: {self.workstation_id}")

        # host compute cluster
        if compute_cluster:
            if self.compute_cluster:
                str_val.append("Compute cluster:")
                str_val.append("  " + str(self.compute_cluster).replace("\n", "\n  "))
            else:
                str_val.append(f"Compute cluster: {None}")

        # user and organization
        if user:
            str_val.append(f"User ID: {self.user_id}")

        if org:
            str_val.append(f"Organization ID: {self.org_id}")

        # hardware
        if hardware:
            if self.hardware:
                str_val.append("Hardware:")
                str_val.append("  " + str(self.hardware).replace("\n", "\n  "))
            else:
                str_val.append(f"Hardware: {None}")

        # environment
        if env:
            str_val.append(f"Environment: {self.env}")

        # debug
        if debug:
            str_val.append(f"Debug: {self.debug}")

        return "\n".join(str_val)


@functools.cache
def get_value():
    """returns a context from environment variables"""
    compute_cluster_str = os.getenv("DEEP_ORIGIN_COMPUTE_CLUSTER", None)
    if compute_cluster_str:
        compute_cluster_dict = json.loads(compute_cluster_str)
        compute_cluster = ComputeClusterContext(
            id=compute_cluster_dict.get("id", None),
            properties=compute_cluster_dict.get("properties", None),
        )
    else:
        compute_cluster = None

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
        workstation_id=os.getenv("DEEP_ORIGIN_BENCH_ID", None),
        user_id=os.getenv("DEEP_ORIGIN_USER_ID", None),
        org_id=os.getenv("DEEP_ORIGIN_ORGANIZATION_ID", None),
        compute_cluster=compute_cluster,
        hardware=hardware,
        env=os.getenv("DEEP_ORIGIN_ENV", None),
        debug=os.getenv("DEEP_ORIGIN_DEBUG", "").lower() == "true",
    )
