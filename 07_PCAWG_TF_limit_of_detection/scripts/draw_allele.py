import numpy as np
from numba import jit

@jit(nopython=True, parallel=False)
# Input x: vafs* tfs
# Output x: an array of 1,0 based on the vafs* tfs
def draw_allele(x):
    for i in range(x.shape[0]):
        x[i] = np.random.binomial(n = 1, p = x[i])
    return x