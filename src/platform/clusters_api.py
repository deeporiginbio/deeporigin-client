"""bridge module to interact with the platform tools api"""

import functools
import sys
from typing import Optional

from beartype import beartype

from deeporigin.platform.utils import _add_functions_to_module

__all__ = _add_functions_to_module(
    module=sys.modules[__name__],
    api_name="ClustersApi",
)


@functools.cache
@beartype
def _get_cluster_id(
    *,
    client=None,
    org_friendly_id: Optional[str] = None,
) -> str:
    """gets a valid cluster ID to run tools on

    this defaults to pulling us-west-2"""

    available_clusters = list_clusters(  # noqa: F821
        client=client,
        org_friendly_id=org_friendly_id,
    )

    cluster = [
        cluster
        for cluster in available_clusters
        if "us-west-2" in cluster.attributes.name
    ]
    cluster = cluster[0]
    return cluster.id
