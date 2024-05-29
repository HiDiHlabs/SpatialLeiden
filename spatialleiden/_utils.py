from typing import TypeAlias

import numpy as np
from scipy.sparse import csr_array, sparray, spmatrix

_GraphArray: TypeAlias = sparray | spmatrix | np.ndarray


def distance2connectivity(distances: sparray | spmatrix) -> sparray | spmatrix:
    if distances.format not in {"csr", "csc", "coo"}:
        connectivity = csr_array(distances)
    else:
        connectivity = distances.copy()
    connectivity.data = 1 - (connectivity.data / connectivity.data.max())
    return connectivity
