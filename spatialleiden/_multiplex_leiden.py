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


def multiplex_leiden(
    latent_neighbors: _GraphArray,
    spatial_neighbors: _GraphArray,
    *,
    directed: tuple[bool, bool] = (True, True),
    use_weights: bool = True,
    n_iterations: int = -1,
    partition_type=la.RBConfigurationVertexPartition,
    layer_ratio: float = 1,
    latent_partition_kwargs: dict | None = None,
    spatial_partition_kwargs: dict | None = None,
    seed: int = 42,
) -> NDArray[np.integer]:
    """
    Partition the nodes using multiplex leiden clustering.

    Parameters
    ----------
    latent_neighbors : scipy.sparse.sparray | scipy.sparse.spmatrix | numpy.ndarray
        Matrix of row-wise neighbor definitions in the latent space
        i.e. c\ :sub:`ij` is the connectivity of i :math:`\\to` j.
    spatial_neighbors : scipy.sparse.sparray | scipy.sparse.spmatrix | numpy.ndarray
        Matrix of row-wise neighbor definitions in the topological space
        i.e. c\ :sub:`ij` is the connectivity of i :math:`\\to` j.
    directed : tuple[bool, bool], optional
        Whether to use a directed graph for latent and topological neighbors,
        respectively.
    use_weights : bool, optional
        Whether to use weights for the edges.
    n_iterations : int, optional
        Number of iterations to run the Leiden algorithm. If the number is negative it
        is run until convergence.
    partition_type : optional
        A :py:class:`leidenalg.VertexPartition.MutableVertexPartition` to be used.
    layer_ratio : tuple[int, int], optional
        The ratio of the weighting of the layers in latent and topological space.
        A higher ratio will increase relevance of the topological neighbors and lead to
        more spatially homogeneous clusters.
    latent_partition_kwargs : dict | None, optional
        Keyword arguments for the latent space partition.
    spatial_partition_kwargs : dict | None, optional
        Keyword arguments for the topological space partition.
    seed : int, optional
        Randoem seed.

    Returns
    -------
    numpy.ndarray[numpy.integer]
        Cluster assignment.
    """

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
        layer_weights=[1, 1 * layer_ratio],
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
    layer_ratio: float = 1,
    latent_distance_key: str = "connectivities",
    spatial_distance_key: str = "spatial_connectivities",
    latent_partition_kwargs: dict | None = None,
    spatial_partition_kwargs: dict | None = None,
    seed: int = 42,
):
    """
    Perform SpatialLeiden clustering.

    This is a wrapper around :py:func:`spatialleiden.multiplex_leiden` that can directly
    use :py:class:`anndata.AnnData` as input and to save results.

    Parameters
    ----------
    adata : anndata.AnnData
    resolution : tuple[int, int], optional
        resolution for the latent space and topological space layer, respectively.
    latent_neighbors : scipy.sparse.sparray | scipy.sparse.spmatrix | numpy.ndarray
        Matrix of row-wise neighbor definitions in the latent space
        i.e. c\ :sub:`ij` is the connectivity of i :math:`\\to` j.
    spatial_neighbors : scipy.sparse.sparray | scipy.sparse.spmatrix | numpy.ndarray
        Matrix of row-wise neighbor definitions in the topological space
        i.e. c\ :sub:`ij` is the connectivity of i :math:`\\to` j.
    key_added : str, optional
        Key to store the clustering results in :py:attr:`anndata.AnnData.obs`
    directed : tuple[bool, bool], optional
        Whether to use a directed graph for latent and topological neighbors,
        respectively.
    use_weights : bool, optional
        Whether to use weights for the edges.
    n_iterations : int, optional
        Number of iterations to run the Leiden algorithm. If the number is negative it
        runs until convergence.
    partition_type : optional
        A :py:class:`leidenalg.VertexPartition.MutableVertexPartition` to be used.
    layer_ratio : tuple[int, int], optional
        The ratio of the weighting of the layers in latent and topological space.
        A higher ratio will increase relevance of the topological neighbors and lead to
        more spatially homogeneous clusters.
    latent_distance_key : str, optional
        Key to use for the latent neighbor connectivities in
        :py:attr:`anndata.AnnData.obsp`. Only used if `latent_neighbors` is `None`.
    spatial_distance_key : str, optional
        Key to use for the spatial neighbor connectivities in
        :py:attr:`anndata.AnnData.obsp`. Only used if `spatial_neighbors` is `None`.
    latent_partition_kwargs : dict | None, optional
        Keyword arguments for the latent space partition.
    spatial_partition_kwargs : dict | None, optional
        Keyword arguments for the topological space partition.
    seed : int, optional
        Random seed.
    """

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

    cluster = multiplex_leiden(
        latent_distances,
        spatial_distances,
        directed=directed,
        use_weights=use_weights,
        n_iterations=n_iterations,
        partition_type=partition_type,
        layer_ratio=layer_ratio,
        spatial_partition_kwargs=spatial_partition_kwargs,
        latent_partition_kwargs=latent_partition_kwargs,
        seed=seed,
    )

    adata.obs[key_added] = cluster
    adata.obs[key_added] = adata.obs[key_added].astype("category")
