from scipy.sparse import csr_array, sparray, spmatrix


def distance2connectivity(distances: sparray | spmatrix) -> sparray | spmatrix:
    """
    Transforms distances to connectivites.

    A sparse distance matrix is transformed to connectivities by calculating
    :math:`1-d/d_{max}`.

    Parameters
    ----------
    distances : scipy.sparse.sparray | scipy.sparse.spmatrix
        Sparse matrix of pairwise distances.

    Returns
    -------
    scipy.sparse.sparray | scipy.sparse.spmatrix
    """
    if distances.format not in {"csr", "csc", "coo"}:
        connectivity = csr_array(distances)
    else:
        connectivity = distances.copy()
    connectivity.data = 1 - (connectivity.data / connectivity.data.max())
    return connectivity
