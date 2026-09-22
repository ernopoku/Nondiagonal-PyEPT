"""Independent determinant algebra, only for tiny-system validation."""
import numpy as np
from itertools import combinations

def apply(det,operators):
    phase=1
    for p,creator in reversed(operators):
        if bool(det>>p&1)==creator:return None,0
        phase*=(-1)**((det&((1<<p)-1)).bit_count());det^=1<<p
    return det,phase

def dagger(op):return tuple((p,not c) for p,c in reversed(op))

def operator_matrix(n,op):
    matrix=np.zeros((2**n,2**n))
    for det in range(2**n):
        d,s=apply(det,op)
        if s:matrix[d,det]=s
    return matrix

def physical_hamiltonian(ints):
    n,o=ints.n,ints.nocc;g=ints.g('p','p','p','p')
    h=np.diag(ints.energy)-np.einsum('piqi->pq',g[:,:o,:,:o])
    mat=np.zeros((2**n,2**n))
    for p in range(n):
        for q in range(n):
            mat+=h[p,q]*operator_matrix(n,((p,1),(q,0)))
    for p,q in combinations(range(n),2):
        for r,s in combinations(range(n),2):
            mat+=g[p,q,r,s]*operator_matrix(n,((p,1),(q,1),(s,0),(r,0)))
    return mat

def reference_states(ints,t):
    ref=np.zeros(2**ints.n);ref[(1<<ints.nocc)-1]=1
    corr=np.zeros_like(ref)
    for i,j in combinations(range(ints.nocc),2):
        for a,b in combinations(range(ints.nocc,ints.n),2):
            corr+=t[i,j,a-ints.nocc,b-ints.nocc]*operator_matrix(ints.n,((a,1),(b,1),(j,0),(i,0)))@ref
    return ref,corr
