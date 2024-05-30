from collections.abc import Callable

import scanpy as sc
from anndata import AnnData

from ._multiplex_leiden import spatialleiden


def _search_resolution(
    fn: Callable[[float], int],
    ncluster: int,
    start: float = 1,
    step: float = 0.1,
    n_iterations: int = 15,
) -> float:
    # adapted from SpaGCN.search_res (https://github.com/jianhuupenn/SpaGCN)
    res = start
    old_ncluster = fn(res)
    iter = 1
    while old_ncluster != ncluster:
        old_sign = 1 if (old_ncluster < ncluster) else -1
        new_ncluster = fn(res + step * old_sign)
        if new_ncluster == ncluster:
            res = res + step * old_sign
            # print(f"Recommended res = {res:.2f}")
            return res
        new_sign = 1 if (new_ncluster < ncluster) else -1
        if new_sign == old_sign:
            res = res + step * old_sign
            # print(f"Res changed to {res:.2f}")
            old_ncluster = new_ncluster
        else:
            step = step / 2
            # print(f"Step changed to {step:.2f}")
        if iter > n_iterations:
            # print("Exact resolution not found")
            # print(f"Recommended res =  {res:.2f}")
            return res
        iter += 1
    # print(f"Recommended res = {res:.2f}")
    return res


def search_resolution_latent(
    adata: AnnData,
    ncluster: int,
    *,
    start: float = 1,
    step: float = 0.1,
    n_iterations: int = 15,
    **kwargs,
) -> float:
    """
    Search the resolution to obtain `n` clusters using Leiden clustering.

    Parameters
    ----------
    adata : anndata.AnnData
    ncluster : int
        Number of clusters.
    start : float, optional
        Starting point for resolution.
    step : float, optional
        Increment if cluster number is incorrect.
    n_iterations : int, optional
        Maximum number of iterations before stopping. If correct number of clusters is
        obtained it will stop early.
    kwargs
        Other keyword arguments are passed to :py:func:`scanpy.tl.leiden`.

    Returns
    -------
    float
        Target resolution.
    """

    def ncluster4res_leiden(resolution: float) -> int:
        sc.tl.leiden(adata, resolution=resolution, **kwargs)
        return adata.obs[key_added].cat.categories.size

    key_added = kwargs.pop("key_added", "leiden")

    return _search_resolution(ncluster4res_leiden, ncluster, start, step, n_iterations)


def search_resolution_spatial(
    adata: AnnData,
    ncluster: int,
    *,
    start: float = 0.4,
    step: float = 0.1,
    n_iterations: int = 15,
    **kwargs,
) -> float:
    """
    Search the resolution of the spatial layer to obtain `n` clusters using
    SpatialLeiden clustering.

    The resolution of the latent space must be provided using ``kwargs``. If you want
    to search both resolutions refer to :py:func:`spatialleiden.search_resolution`.

    Parameters
    ----------
    adata : anndata.AnnData
    ncluster : int
        Number of clusters.
    start : float, optional
        Starting point for resolution.
    step : float, optional
        Increment if cluster number is incorrect.
    n_iterations : int, optional
        Maximum number of iterations before stopping. If correct number of clusters is
        obtained it will stop early.
    kwargs
        Other keyword arguments are passed to :py:func:`spatialleiden.spatialleiden`.

    Returns
    -------
    float
        Target resolution for the topological space.
    """

    def ncluster4res_spatialleiden(resolution: float) -> int:
        spatialleiden(
            adata,
            resolution=(resolution_user[0], resolution),
            key_added=key_added,
            **kwargs,
        )
        return adata.obs[key_added].cat.categories.size

    key_added = kwargs.pop("key_added", "spatialleiden")
    resolution_user = kwargs.pop("resolution", (1, 1))

    return _search_resolution(
        ncluster4res_spatialleiden, ncluster, start, step, n_iterations
    )


def search_resolution(
    adata: AnnData,
    ncluster: int,
    *,
    start: tuple[float, float] = (1.0, 0.4),
    step: float = 0.1,
    n_iterations: int = 15,
    latent_kwargs: dict | None = None,
    spatial_kwargs: dict | None = None,
) -> tuple[float, float]:
    """
    Search the resolutions of the topological and latent layer to obtain `n` clusters
    using SpatialLeiden clustering.

    Parameters
    ----------
    adata : anndata.AnnData
    ncluster : int
        Number of clusters.
    start : float, optional
        Starting point for resolution.
    step : float, optional
        Increment if cluster number is incorrect.
    n_iterations : int, optional
        Maximum number of iterations before stopping. If correct number of clusters is
        obtained it will stop early.
    latent_kwargs : dict | None, optional
        Keyword arguments passed to :py:func:`scanpy.tl.leiden`.
    spatial_kwargs : dict | None, optional
        Keyword arguments passed to :py:func:`spatialleiden.spatialleiden`.

    Returns
    -------
    tuple[float, float]
        Target resolution for the latent and topological space.
    """
    latent_kwargs = dict() if latent_kwargs is None else latent_kwargs
    spatial_kwargs = dict() if spatial_kwargs is None else spatial_kwargs

    resolution_latent = search_resolution_latent(
        adata,
        ncluster,
        start=start[0],
        step=step,
        n_iterations=n_iterations,
        **latent_kwargs,
    )

    spatial_kwargs["resolution"] = (resolution_latent, 1)
    resolution_spatial = search_resolution_spatial(
        adata,
        ncluster,
        start=start[1],
        step=step,
        n_iterations=n_iterations,
        **spatial_kwargs,
    )
    return (resolution_latent, resolution_spatial)
