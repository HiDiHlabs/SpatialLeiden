from typing import TypeAlias

import leidenalg as la
import numpy as np
from anndata import AnnData
from igraph import Graph
from numpy.typing import NDArray
from scipy.sparse import find, sparray, spmatrix

_GraphArray: TypeAlias = sparray | spmatrix | np.ndarray


def _build_igraph(adjacency: _GraphArray, *, directed: bool = True) -> Graph:
    # adapted from scanpy https://github.com/scverse/scanpy
    sources, targets, weights = find(adjacency)
    g = Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])
    g.add_edges(list(zip(sources, targets)))
    g.es["weight"] = weights
    if g.vcount() != adjacency.shape[0]:
        raise RuntimeError(
            f"The constructed graph has only {g.vcount()} nodes. "
            "Your adjacency matrix contained redundant nodes."
        )
    return g


def leiden_multiplex(
    latent_neighbors: _GraphArray,
    spatial_neighbors: _GraphArray,
    *,
    directed: tuple[bool, bool] = (True, True),
    use_weights: bool = True,
    n_iterations: int = -1,
    partition_type=la.RBConfigurationVertexPartition,
    layer_weights: tuple[int, int] = (1, 1),
    latent_partition_kwargs: dict | None = None,
    spatial_partition_kwargs: dict | None = None,
    seed: int = 42,
) -> NDArray[np.integer]:

    adjacency_latent = _build_igraph(latent_neighbors, directed=directed[0])
    adjacency_spatial = _build_igraph(spatial_neighbors, directed=directed[1])

    # parameterise the partitions
    if spatial_partition_kwargs is None:
        spatial_partition_kwargs = dict()
    if latent_partition_kwargs is None:
        latent_partition_kwargs = dict()

    if use_weights:
        spatial_partition_kwargs["weights"] = "weight"
        latent_partition_kwargs["weights"] = "weight"

    latent_part = partition_type(adjacency_latent, **latent_partition_kwargs)
    spatial_part = partition_type(adjacency_spatial, **spatial_partition_kwargs)

    optimiser = la.Optimiser()
    optimiser.set_rng_seed(seed)

    _ = optimiser.optimise_partition_multiplex(
        [latent_part, spatial_part],
        layer_weights=list(layer_weights),
        n_iterations=n_iterations,
    )

    return np.array(latent_part.membership)


def spatialleiden(
    adata: AnnData,
    *,
    resolution: tuple[int, int] = (1, 1),
    latent_neighbors: _GraphArray | None = None,
    spatial_neighbors: _GraphArray | None = None,
    key_added: str = "spatialleiden",
    directed: tuple[bool, bool] = (True, True),
    use_weights: bool = True,
    n_iterations: int = -1,
    partition_type=la.RBConfigurationVertexPartition,
    layer_weights: tuple[int, int] = (1, 1),
    latent_distance_key: str = "connectivities",
    spatial_distance_key: str = "spatial_connectivities",
    latent_partition_kwargs: dict | None = None,
    spatial_partition_kwargs: dict | None = None,
    seed: int = 42,
):

    if latent_neighbors is None:
        latent_distances = adata.obsp[latent_distance_key]
    if spatial_neighbors is None:
        spatial_distances = adata.obsp[spatial_distance_key]

    if latent_partition_kwargs is None:
        latent_partition_kwargs = dict()
    if spatial_partition_kwargs is None:
        spatial_partition_kwargs = dict()

    spatial_partition_kwargs["resolution"], latent_partition_kwargs["resolution"] = (
        resolution
    )

    cluster = leiden_multiplex(
        latent_distances,
        spatial_distances,
        directed=directed,
        use_weights=use_weights,
        n_iterations=n_iterations,
        partition_type=partition_type,
        layer_weights=layer_weights,
        spatial_partition_kwargs=spatial_partition_kwargs,
        latent_partition_kwargs=latent_partition_kwargs,
        seed=seed,
    )

    adata.obs[key_added] = cluster
    adata.obs[key_added] = adata.obs[key_added].astype("category")
