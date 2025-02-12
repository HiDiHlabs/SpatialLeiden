from collections.abc import Collection, Iterable, Mapping
from types import NoneType
from typing import Any, Type, TypeAlias, TypeVar

import leidenalg as la
import numpy as np
from anndata import AnnData
from igraph import Graph
from leidenalg.VertexPartition import MutableVertexPartition
from mudata import MuData
from numpy.typing import NDArray
from scipy.sparse import find, sparray, spmatrix

_GraphArray: TypeAlias = sparray | spmatrix | np.ndarray


def _build_igraph(adjacency: _GraphArray, *, directed: bool = True) -> Graph:
    sources, targets, weights = find(adjacency)
    g = Graph(
        n=adjacency.shape[0],
        edges=zip(sources, targets),
        directed=directed,
        edge_attrs={"weight": weights},
    )
    return g


def multiplex_leiden(
    *neighbors: _GraphArray,
    directed: bool | Collection[bool] = True,
    use_weights: bool | Collection[bool] = True,
    n_iterations: int = -1,
    partition_type: Type[MutableVertexPartition] = la.RBConfigurationVertexPartition,
    layer_weights: float | Collection[float] = 1,
    partition_kwargs: dict | None | Collection[dict | None] = None,
    seed: int = 42,
) -> NDArray[np.integer]:
    """
    Partition the nodes using multiplex leiden clustering.

    For more information on multiplex leiden clustering read the
    `leidenalg documentation <https://leidenalg.readthedocs.io/en/stable/multiplex.html>`_.

    Parameters
    ----------
    neighbors : scipy.sparse.sparray | scipy.sparse.spmatrix | numpy.ndarray
        Matrices of row-wise neighbor definitions for the different layers
        i.e. c\\ :sub:`ij` is the connectivity of i :math:`\\to` j.
    directed : bool | collections.abc.Collection[bool], optional
        Whether to use a directed graph for each layer, respectively.
    use_weights : bool | collections.abc.Collection[bool], optional
        Whether to use weights for the edges of each layer, respectively.
    n_iterations : int, optional
        Number of iterations to run the Leiden algorithm. If the number is negative it
        is run until convergence.
    partition_type : optional
        A :py:class:`leidenalg.VertexPartition.MutableVertexPartition` to be used.
    layer_weights : float | collections.abc.Collection[float], optional
        The weights of each layer, respectively.
    partition_kwargs : dict | None | collections.abc.Collection[dict | None], optional
        Keyword arguments for the partition of each layer.
    seed : int, optional
        Random seed.

    Returns
    -------
    numpy.ndarray[numpy.integer]
        Cluster assignment.
    """

    def check_length(x, type, n) -> Collection | list:
        if isinstance(x, type):
            x = [x] * n
        elif len(x) != n:
            raise ValueError("")
        return x

    n_layers = len(neighbors)

    directed = check_length(directed, bool, n_layers)
    use_weights = check_length(use_weights, bool, n_layers)
    layer_weights = check_length(layer_weights, float, n_layers)
    partition_kwargs = check_length(partition_kwargs, (dict, NoneType), n_layers)

    layers = [
        _build_igraph(n, directed=d) for n, d in zip(neighbors, directed, strict=True)
    ]

    # parameterise the partitions
    partition_kwargs_ls = list()
    for p_kwargs, with_weight in zip(partition_kwargs, use_weights, strict=True):
        p_kwargs = p_kwargs if p_kwargs is not None else dict()
        if with_weight:
            p_kwargs["weights"] = "weight"
        partition_kwargs_ls.append(p_kwargs)

    partitions = [
        partition_type(layer, **kwargs)
        for layer, kwargs in zip(layers, partition_kwargs_ls, strict=True)
    ]

    optimiser = la.Optimiser()
    optimiser.set_rng_seed(seed)

    _ = optimiser.optimise_partition_multiplex(
        partitions,
        layer_weights=list(layer_weights),
        n_iterations=n_iterations,
    )

    return np.array(partitions[0].membership)


def spatialleiden(
    adata: AnnData,
    *,
    resolution: tuple[float, float] = (1, 1),
    latent_neighbors: _GraphArray | None = None,
    spatial_neighbors: _GraphArray | None = None,
    key_added: str = "spatialleiden",
    directed: tuple[bool, bool] = (True, True),
    use_weights: tuple[bool, bool] = (True, True),
    n_iterations: int = -1,
    partition_type: Type[MutableVertexPartition] = la.RBConfigurationVertexPartition,
    layer_ratio: float = 1,
    latent_neighbors_key: str = "connectivities",
    spatial_neighbors_key: str = "spatial_connectivities",
    latent_partition_kwargs: dict[str, Any] | None = None,
    spatial_partition_kwargs: dict[str, Any] | None = None,
    seed: int = 42,
):
    """
    Perform SpatialLeiden clustering.

    This is a wrapper around :py:func:`spatialleiden.multiplex_leiden` that uses
    :py:class:`anndata.AnnData` as input and works with two layers; one latent space
    and one spatial layer.

    Parameters
    ----------
    adata : anndata.AnnData
    resolution : tuple[float, float], optional
        Resolution for the latent space and spatial layer, respectively.
    latent_neighbors : scipy.sparse.sparray | scipy.sparse.spmatrix | numpy.ndarray
        Matrix of row-wise neighbor definitions in the latent space layer
        i.e. c\\ :sub:`ij` is the connectivity of i :math:`\\to` j.
    spatial_neighbors : scipy.sparse.sparray | scipy.sparse.spmatrix | numpy.ndarray
        Matrix of row-wise neighbor definitions in the spatial layer
        i.e. c\\ :sub:`ij` is the connectivity of i :math:`\\to` j.
    key_added : str, optional
        Key to store the clustering results in :py:attr:`anndata.AnnData.obs`
    directed : tuple[bool, bool], optional
        Whether to use a directed graph for latent space and spatial neighbors,
        respectively.
    use_weights : tuple[bool, bool], optional
        Whether to use weights for the edges for latent space and spatial neighbors,
        respectively.
    n_iterations : int, optional
        Number of iterations to run the Leiden algorithm. If the number is negative it
        runs until convergence.
    partition_type : typing.Type[leidenalg.VertexPartition.MutableVertexPartition], optional
        A :py:class:`leidenalg.VertexPartition.MutableVertexPartition` to be used.
    layer_ratio : float, optional
        The ratio of the weighting of the layers; latent space vs spatial.
        A higher ratio will increase relevance of the spatial neighbors and lead to
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
        Keyword arguments for the spatial partition.
    seed : int, optional
        Random seed.
    """

    if latent_neighbors is None:
        latent_connectivities = adata.obsp[latent_neighbors_key]
    if spatial_neighbors is None:
        spatial_connectivities = adata.obsp[spatial_neighbors_key]

    if latent_partition_kwargs is None:
        latent_partition_kwargs = dict()
    if spatial_partition_kwargs is None:
        spatial_partition_kwargs = dict()

    latent_partition_kwargs["resolution_parameter"] = resolution[0]
    spatial_partition_kwargs["resolution_parameter"] = resolution[1]

    cluster = multiplex_leiden(
        latent_connectivities,
        spatial_connectivities,
        directed=directed,
        use_weights=use_weights,
        n_iterations=n_iterations,
        partition_type=partition_type,
        layer_weights=[1, layer_ratio],
        partition_kwargs=[latent_partition_kwargs, spatial_partition_kwargs],
        seed=seed,
    )

    adata.obs[key_added] = cluster
    adata.obs[key_added] = adata.obs[key_added].astype("category")


def spatialleiden_multiomics(
    mdata: MuData,
    *,
    resolution: float | Mapping[str, float] = 1,
    key_added: str = "spatialleiden",
    directed: bool | Mapping[str, bool] = True,
    use_weights: bool | Mapping[str, bool] = True,
    n_iterations: int = -1,
    partition_type: Type[MutableVertexPartition] = la.RBConfigurationVertexPartition,
    layer_weights: float | Mapping[str, float] = 1,
    neighbors_key: str | Mapping[str, str] = "connectivities",
    spatial_neighbors_key: str = "spatial_connectivities",
    partition_kwargs: None | Mapping[str, dict[str, Any]] = None,
    seed: int = 42,
):
    """
    Perform Multi-Omics SpatialLeiden clustering.

    This is a wrapper around :py:func:`spatialleiden.multiplex_leiden` that uses
    :py:class:`mudata.MuData` as input and works with multiple layers; one for each
    modality and one for the spatial layer.

    Parameters
    ----------
    mdata : mudata.MuData
    resolution : float, collections.abc.Mapping[str, float], optional
        Resolution for the neighbor graphs of the different modalities and the spactial
        layer.
    key_added : str, optional
        Key to store the clustering results in :py:attr:`mudata.MuData.obs`.
    directed: bool | collections.abc.Mapping[str, bool], optional
        Whether to use a directed graph for the neighbor graphs of the modalities and
        the spatial layer.
    use_weights: bool | collections.abc.Mapping[str, bool], optional
        Whether to use a weighted edges for the neighbor graphs of the modalities and
        the spatial layer.
    n_iterations : int, optional
        Number of iterations to run the Leiden algorithm. If the number is negative it
        runs until convergence.
    partition_type: typing.Type[leidenalg.VertexPartition.MutableVertexPartition], optional
        A :py:class:`leidenalg.VertexPartition.MutableVertexPartition` to be used.
    layer_weights: float | collections.abc.Mapping[str, float], optional
        The weighting of the different layers.
    neighbors_key: str | collections.abc.Mapping[str, str], optional
        Key(s) used to lookup the neighbor graphs for the different modalities in the
        corresponding :py:attr:`anndata.Anndata.obsp`.
    spatial_neighbors_key: str, optional
        Key used to lookup the spatial neighbors graph in :py:attr:`mudata.Mudata.obsp`.
    partition_kwargs: None | collections.abc.Mapping[str, dict[str, typing.Any]], optional
        Keyword arguments for the modality and spatial partitions.
    seed : int, optional
        Random seed.
    """
    T = TypeVar("T")

    def value_or_mapping_to_list(
        x: T | Mapping[str, T], /, modalities: Iterable[str]
    ) -> list[T]:
        if isinstance(x, Mapping):
            return [x[m] for m in modalities]
        else:
            return [x for _ in modalities]

    spatial_connectivities = mdata.obsp[spatial_neighbors_key]
    if isinstance(neighbors_key, str):
        modal_connectivities = [
            mdata.mod[m].obsp[neighbors_key] for m in mdata.mod_names
        ]
    else:
        modal_connectivities = [
            mdata.mod[m].obsp[neighbors_key[m]] for m in mdata.mod_names
        ]

    modalities = mdata.mod_names + ["spatial"]

    partition_kwargs_ls: list[dict[str, Any]]
    if partition_kwargs is None:
        partition_kwargs_ls = [dict() for _ in modalities]
    else:
        partition_kwargs_ls = [partition_kwargs.get(m, dict()) for m in modalities]

    if isinstance(resolution, (int, float)):
        for p_kwargs in partition_kwargs_ls:
            p_kwargs["resolution_parameter"] = resolution
    else:
        for mod, r in resolution.items():
            i = modalities.index(mod)
            partition_kwargs_ls[i]["resolution_parameter"] = r

    cluster = multiplex_leiden(
        *modal_connectivities,
        spatial_connectivities,
        directed=value_or_mapping_to_list(directed, modalities),
        use_weights=value_or_mapping_to_list(use_weights, modalities),
        n_iterations=n_iterations,
        partition_type=partition_type,
        layer_weights=value_or_mapping_to_list(layer_weights, modalities),
        partition_kwargs=partition_kwargs_ls,
        seed=seed,
    )

    mdata.obs[key_added] = cluster
    mdata.obs[key_added] = mdata.obs[key_added].astype("category")
