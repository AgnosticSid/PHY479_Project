

import numpy as np


### Bipartite Nearest Neighbor Hopping Model
### ... a_i ~t1~ b_i+1/2 ~t2~ a_i+1 ~t1~ ...

# System size = N ab pairs = 2N

N = 4

def make_H(N, t1, t2, bc='periodic'):
    H = np.zeros((2*N, 2*N), dtype=np.double)

    o_diag1 = np.empty((2*N - 1))
    o_diag1[0::2] = -t1 * np.ones(N)
    o_diag1[1::2] = -t2 * np.ones(N-1)

    H += np.diag(o_diag1, k=-1) + np.diag(o_diag1, k=1)
    H[0, 2*N - 1] = -t2 if bc=='periodic' else 0
    H[2*N - 1, 0] = -t2

    return H

def eigenvalues(H):
    return np.sort(np.unique(np.linalg.eig(H)[0]))

def energies(N, t1, t2):
    e1 = lambda k: np.abs(t1 * np.exp(1j*k) + t2)
    e2 = lambda k: -np.abs(t1 * np.exp(1j*k) + t2)
    k = (2 * np.arange(4) * np.pi / 4).astype(float)
    return np.sort(np.unique(np.concatenate((e1(k), e2(k)), axis=0)))

# Periodic BC's

eigs = eigenvalues(make_H(N, 1.0, 0.5))
eigs2 = energies(N, 1.0, 0.5)

print(make_H(N, 1.0, 0.5))
print(eigs, eigs2)

