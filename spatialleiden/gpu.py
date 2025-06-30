from __future__ import annotations

from collections.abc import Collection, Iterable, Mapping
from typing import TYPE_CHECKING, TypeAlias, TypeVar

import numpy as np
from anndata import AnnData
from mudata import MuData
from numpy.typing import NDArray
from scipy.sparse import find, sparray, spmatrix

if TYPE_CHECKING:
    from cugraph import Graph

DEFAULT_ITERATIONS: int = 100

_GraphArray: TypeAlias = sparray | spmatrix | np.ndarray


def _build_cugraph(*adjacencies, use_weights: bool = True) -> Graph:
    import cudf
    from cugraph import Graph

    def _get_edgelist(adjacency, layer) -> cudf.DataFrame:
        sources, targets, weights = find(adjacency)
        return cudf.DataFrame(
            {
                "source": sources,
                "destination": targets,
                "weights": weights,
                "layer": layer,
            }
        )

    g = Graph(directed=False)

    edge_list = cudf.concat(
        [_get_edgelist(adjacency, i) for i, adjacency in enumerate(adjacencies)]
    )
    # with warnings.catch_warnings():
    #     warnings.simplefilter("ignore")
    if use_weights:
        g.from_cudf_edgelist(
            edge_list,
            source="source",
            destination="destination",
            weight="weights",
            edge_type="layer",
        )
    else:
        g.from_cudf_edgelist(
            edge_list, source="source", destination="destination", edge_type="layer"
        )
    return g


def multiplex_leiden(
    *neighbors: _GraphArray,
    resolutions: float | Collection[float] = 1,
    use_weights: bool = True,
    n_iterations: int = DEFAULT_ITERATIONS,
    layer_weights: float | Collection[float] = 1,
    random_state: int = 42,
    **kwargs,
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
    use_weights : bool | collections.abc.Collection[bool], optional
        Whether to use weights for the edges of each layer, respectively.
    n_iterations : int, optional
        Number of iterations to run the Leiden algorithm. If the number is negative it
        is run until convergence.
    layer_weights : float | collections.abc.Collection[float], optional
        The weights of each layer, respectively.
    random_state : int, optional
        Random seed.

    Returns
    -------
    numpy.ndarray[numpy.integer]
        Cluster assignment.
    """
    import cugraph

    def check_length(x, type, n: int) -> Collection:
        if isinstance(x, type):
            x = [x] * n
        elif len(x) != n:
            raise ValueError("Incorrect number of parameters")
        return x

    n_layers = len(neighbors)

    resolutions = check_length(resolutions, float, n_layers)
    layer_weights = check_length(layer_weights, float, n_layers)

    graph = _build_cugraph(*neighbors, use_weights=use_weights)

    partitioning, modularity = cugraph.leiden(
        graph,
        max_iter=n_iterations,
        resolutions=resolutions,
        layer_weights=layer_weights,
        random_state=random_state,
        **kwargs,
    )

    return partitioning["partition"].values_host


def leiden(
    adata: AnnData,
    *,
    resolution: float = 1,
    neighbors: _GraphArray | None = None,
    key_added: str = "leiden",
    use_weights: bool = True,
    n_iterations: int = DEFAULT_ITERATIONS,
    neighbors_key: str = "connectivities",
    random_state: int = 42,
    **kwargs,
):
    """
    Perform Leiden clustering.

    Parameters
    ----------
    adata : anndata.AnnData
    resolution : float, optional
        Resolution for the partition. Controls the coarseness of the clustering.
    neighbors : scipy.sparse.sparray | scipy.sparse.spmatrix | numpy.ndarray
        Matrix of row-wise neighbor definitions
        i.e. c\\ :sub:`ij` is the connectivity of i :math:`\\to` j.
    key_added : str, optional
        Key to store the clustering results in :py:attr:`anndata.AnnData.obs`
    directed : bool, optional
        Whether to use a directed graph.
    use_weights : bool, optional
        Whether to use weights for the edges.
    n_iterations : int, optional
        Number of iterations to run the Leiden algorithm. If the number is negative it
        runs until convergence.
    partition_type : typing.Type[leidenalg.VertexPartition.MutableVertexPartition], optional
        A :py:class:`leidenalg.VertexPartition.MutableVertexPartition` to be used.
    neighbors_key : str, optional
        Key to use for the neighbor connectivities in
        :py:attr:`anndata.AnnData.obsp`. Only used if `neighbors` is `None`.
    random_state : int, optional
        Random seed.
    kwargs
        Keyword arguments for :py:func:`cugraph.leiden`.
    """
    import cugraph

    if neighbors is None:
        neighbors = adata.obsp[neighbors_key]

    graph = _build_cugraph(neighbors, use_weights=use_weights)

    partitioning, modularity = cugraph.leiden(
        graph,
        max_iter=n_iterations,
        resolution=resolution,
        random_state=random_state,
        **kwargs,
    )

    adata.obs[key_added] = partitioning["partition"].values_host
    adata.obs[key_added] = adata.obs[key_added].astype("category")


def spatialleiden(
    adata: AnnData,
    *,
    resolution: float | tuple[float, float] = 1,
    latent_neighbors: _GraphArray | None = None,
    spatial_neighbors: _GraphArray | None = None,
    key_added: str = "spatialleiden",
    use_weights: bool = True,
    n_iterations: int = DEFAULT_ITERATIONS,
    layer_ratio: float = 1,
    latent_neighbors_key: str = "connectivities",
    spatial_neighbors_key: str = "spatial_connectivities",
    random_state: int = 42,
    **kwargs,
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
    use_weights : tuple[bool, bool], optional
        Whether to use weights for the edges for latent space and spatial neighbors,
        respectively.
    n_iterations : int, optional
        Number of iterations to run the Leiden algorithm. If the number is negative it
        runs until convergence.
    layer_ratio : float, optional
        The ratio of the weighting of the layers; latent space vs spatial.
        A higher ratio will increase relevance of the spatial neighbors and lead to
        more spatially homogeneous clusters.
    latent_neighbors_key : str, optional
        Key to use for the latent neighbor connectivities in
        :py:attr:`anndata.AnnData.obsp`. Only used if `latent_neighbors` is `None`.
    spatial_neighbors_key : str, optional
        Key to use for the spatial neighbor connectivities in
        :py:attr:`anndata.AnnData.obsp`. Only used if `spatial_neighbors` is `None`.
    random_state : int, optional
        Random seed.
    kwargs
    """

    if latent_neighbors is None:
        latent_neighbors = adata.obsp[latent_neighbors_key]
    if spatial_neighbors is None:
        spatial_neighbors = adata.obsp[spatial_neighbors_key]

    adata.obs[key_added] = multiplex_leiden(
        latent_neighbors,
        spatial_neighbors,
        resolutions=resolution,
        use_weights=use_weights,
        n_iterations=n_iterations,
        layer_weights=[1, layer_ratio],
        random_state=random_state,
        **kwargs,
    )

    adata.obs[key_added] = adata.obs[key_added].astype("category")


def spatialleiden_multimodal(
    mdata: MuData,
    *,
    resolution: float | Mapping[str, float] = 1,
    key_added: str = "spatialleiden",
    use_weights: bool = True,
    n_iterations: int = DEFAULT_ITERATIONS,
    layer_weights: float | Mapping[str, float] = 1,
    neighbors_key: str | Mapping[str, str] = "connectivities",
    spatial_neighbors_key: str = "spatial_connectivities",
    random_state: int = 42,
    **kwargs,
):
    """
    Perform multimodal SpatialLeiden clustering.

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
    use_weights: bool | collections.abc.Mapping[str, bool], optional
        Whether to use a weighted edges for the neighbor graphs of the modalities and
        the spatial layer.
    n_iterations : int, optional
        Number of iterations to run the Leiden algorithm. If the number is negative it
        runs until convergence.
    layer_weights: float | collections.abc.Mapping[str, float], optional
        The weighting of the different layers.
    neighbors_key: str | collections.abc.Mapping[str, str], optional
        Key(s) used to lookup the neighbor graphs for the different modalities in the
        corresponding :py:attr:`anndata.AnnData.obsp`.
    spatial_neighbors_key: str, optional
        Key used to lookup the spatial neighbors graph in :py:attr:`mudata.MuData.obsp`.
    partition_kwargs: None | collections.abc.Mapping[str, dict[str, typing.Any]], optional
        Keyword arguments for the modality and spatial partitions.
    random_state : int, optional
        Random seed.
    kwargs
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

    mdata.obs[key_added] = multiplex_leiden(
        *modal_connectivities,
        spatial_connectivities,
        resolutions=value_or_mapping_to_list(resolution, modalities),
        use_weights=use_weights,
        n_iterations=n_iterations,
        layer_weights=value_or_mapping_to_list(layer_weights, modalities),
        random_state=random_state,
        **kwargs,
    )

    mdata.obs[key_added] = mdata.obs[key_added].astype("category")
