import cement

from ..context import (
    get_value as get_context,
)

__all__ = [
    "CONTROLLERS",
    "ContextController",
]


class ContextController(cement.Controller):
    """Context controller"""

    class Meta:
        label = "context"
        stacked_on = "base"
        stacked_type = "nested"
        help = "Get the context for the Deep Origin ComputeBench"
        description = "Get the context for the Deep Origin ComputeBench, such as the ID of the bench, user and organization, host, and hardware blueprint."
        arguments = []

    @cement.ex(hide=True)
    def _default(self):
        """Default action. returns the context"""
        context = get_context()

        # bench id
        print(f"Bench ID: {context.bench_id}")

        # user and organization
        print(f"User ID: {context.user_id}")
        print(f"Organization ID: {context.org_id}")

        # host
        if context.host:
            print("Host:")
            print(f"  ID: {context.host.id}")
            if context.host.properties:
                print("  Properties:")
                for key, value in context.host.properties.items():
                    print(f"    {key}: {value}")
            else:
                print(f"  Properties: {None}")
        else:
            print(f"Host: {None}")

        # hardware
        if context.hardware:
            print("Hardware:")

            # CPU
            if context.hardware.cpu:
                print("  Processing:")
                if context.hardware.cpu.quantity_vcpu is not None:
                    print(
                        f"    Quantity: {context.hardware.cpu.quantity_vcpu:.1f} vCPU"
                    )
                else:
                    print("    Quantity: N/A")
                print(f"    Architecture: {context.hardware.cpu.architecture or 'N/A'}")
                if context.hardware.cpu.memory_gb is not None:
                    print(f"    Memory: {context.hardware.cpu.memory_gb:.1f} GB")
                else:
                    print("    Memory: N/A")
            else:
                print("  Processing: N/A")

            # GPU
            if context.hardware.gpu:
                print("  Accelerated processing:")
                if context.hardware.gpu.quantity_gpu is not None:
                    print(f"    Quantity: {context.hardware.gpu.quantity_gpu:.1f} GPU")
                else:
                    print("    Quantity: None")
                if context.hardware.gpu.architecture:
                    print("    Architecture:")
                    print(f"      Vendor: {context.hardware.gpu.architecture.vendor}")
                    print(
                        f"      Microarchitecture: {context.hardware.gpu.architecture.microarchitecture}"
                    )
                    print(f"      Model: {context.hardware.gpu.architecture.model}")
                else:
                    print(f"    Architecture: {None}")
                if context.hardware.gpu.memory_gb is not None:
                    print(f"    Memory: {context.hardware.gpu.memory_gb:.1f} GB")
                else:
                    print("    Memory: None")
            else:
                print(f"  Accelerated processing: {None}")
        else:
            print(f"Hardware: {None}")

        # environment
        print(f"Environment: {context.env}")

        # debug
        print(f"Debug: {context.debug}")


CONTROLLERS = [
    ContextController,
]
