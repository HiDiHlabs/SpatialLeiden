from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING
from warnings import warn

from ._multiplex_leiden import leiden, spatialleiden

if TYPE_CHECKING:
    from anndata import AnnData


def _search_resolution(
    fn: Callable[[float], int],
    n_clusters: int,
    start: float = 1,
    step: float = 0.1,
    n_iter: int = 15,
) -> float:
    if n_iter <= 2:
        raise ValueError("At least 2 iterations are necessary")
    resolution = start
    n = fn(resolution)
    i = 1
    while n != n_clusters:
        if i >= n_iter:
            warn(
                "Correct resolution not found. Consider increasing the number of "
                "iterations or adjusting the step size."
            )
            break
        sign = 1 if n < n_clusters else -1
        resolution += step * sign  # TODO what happens when approaching zero
        n = fn(resolution)
        new_sign = 1 if n < n_clusters else -1
        if new_sign != sign:
            step /= 2
        i += 1

    return resolution


def search_resolution_latent(
    adata: AnnData,
    n_clusters: int,
    *,
    start: float = 1,
    step: float = 0.1,
    n_iter: int = 15,
    **kwargs,
) -> float:
    """
    Search the resolution to obtain `n` clusters using Leiden clustering.

    Parameters
    ----------
    adata : anndata.AnnData
    n_clusters : int
        Number of clusters.
    start : float, optional
        Starting point for resolution.
    step : float, optional
        Increment if cluster number is incorrect.
    n_iter : int, optional
        Maximum number of iterations before stopping. If correct number of clusters is
        obtained it will stop early.
    kwargs
        Other keyword arguments are passed to :py:func:`spatialleiden.leiden`.

    Returns
    -------
    float
        Target resolution.
    """

    def ncluster4res_leiden(resolution: float) -> int:
        leiden(adata, resolution=resolution, key_added=key_added, **kwargs)
        return adata.obs[key_added].cat.categories.size

    key_added = kwargs.pop("key_added", "leiden")

    return _search_resolution(ncluster4res_leiden, n_clusters, start, step, n_iter)


def search_resolution_spatial(
    adata: AnnData,
    n_clusters: int,
    *,
    start: float = 0.4,
    step: float = 0.1,
    n_iter: int = 15,
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
    n_clusters : int
        Number of clusters.
    start : float, optional
        Starting point for resolution.
    step : float, optional
        Increment if cluster number is incorrect.
    n_iter : int, optional
        Maximum number of iterations before stopping. If correct number of clusters is
        obtained it will stop early.
    kwargs
        Other keyword arguments are passed to :py:func:`spatialleiden.spatialleiden`.

    Returns
    -------
    float
        Target resolution for the spatial layer.
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
        ncluster4res_spatialleiden, n_clusters, start, step, n_iter
    )


def search_resolution(
    adata: AnnData,
    n_clusters: int,
    *,
    start: tuple[float, float] = (1.0, 0.4),
    step: float = 0.1,
    n_iter: int = 15,
    latent_kwargs: dict | None = None,
    spatial_kwargs: dict | None = None,
) -> tuple[float, float]:
    """
    Search the resolutions of the spatial and latent space layer to obtain `n` clusters
    using SpatialLeiden clustering.

    Parameters
    ----------
    adata : anndata.AnnData
    n_clusters : int
        Number of clusters.
    start : tuple[float, float], optional
        Starting points for resolution, latent and spatial layer, respectively.
    step : float, optional
        Increment if cluster number is incorrect.
    n_iter : int, optional
        Maximum number of iterations before stopping. If correct number of clusters is
        obtained it will stop early.
    latent_kwargs : dict | None, optional
        Keyword arguments passed to :py:func:`spatialleiden.leiden`.
    spatial_kwargs : dict | None, optional
        Keyword arguments passed to :py:func:`spatialleiden.spatialleiden`.

    Returns
    -------
    tuple[float, float]
        Target resolution for the latent space and spatial layer.
    """
    latent_kwargs = dict() if latent_kwargs is None else latent_kwargs
    spatial_kwargs = dict() if spatial_kwargs is None else spatial_kwargs

    resolution_latent = search_resolution_latent(
        adata, n_clusters, start=start[0], step=step, n_iter=n_iter, **latent_kwargs
    )

    spatial_kwargs["resolution"] = (resolution_latent, 1)
    resolution_spatial = search_resolution_spatial(
        adata, n_clusters, start=start[1], step=step, n_iter=n_iter, **spatial_kwargs
    )
    return (resolution_latent, resolution_spatial)
