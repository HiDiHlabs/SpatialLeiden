from importlib.metadata import PackageNotFoundError, version

from ._multiplex_leiden import multiplex_leiden, spatialleiden
from ._resolution_search import (
    search_resolution,
    search_resolution_latent,
    search_resolution_spatial,
)
from ._utils import distance2connectivity

try:
    __version__ = version("spatialleiden")
except PackageNotFoundError:
    __version__ = "unknown version"

del PackageNotFoundError, version


__all__ = [
    "multiplex_leiden",
    "spatialleiden",
    "search_resolution",
    "search_resolution_latent",
    "search_resolution_spatial",
    "distance2connectivity",
]
