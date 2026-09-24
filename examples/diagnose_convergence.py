"""Reproduce H0-preconditioned Davidson iterations without early exit.

Diagnostic only: energy-change convergence is intentionally NOT used to return
scientific results. Use EPT.kernel() for residual-controlled results.
"""
import argparse,json
from pathlib import Path
import numpy as np
from pyscf import gto,scf
from nondiagonal_ept import EPT
from validate_reference import GEOMETRIES

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--molecule',choices=list(GEOMETRIES),default='H2O')
    parser.add_argument('--target',type=int,default=3,help='Original zero-based spatial MO')
    parser.add_argument('--method',choices=['NRP3','NRQ3','NRL3'],default='NRL3')
    parser.add_argument('--cycles',type=int,default=10)
    parser.add_argument('--output',default='convergence_trace.json')
    args=parser.parse_args()
    frozen=1 if args.molecule in ('HF','H2O') else 2
    mf=scf.RHF(gto.M(atom=GEOMETRIES[args.molecule],basis='cc-pvtz',verbose=0)).run(conv_tol=1e-13)
    h=EPT(mf,args.method,frozen=frozen).hamiltonian
    original=h.ints.original_mos[h.ints.spatial[h.simple]]
    found=np.flatnonzero(original==args.target)
    if len(found)!=1:raise ValueError('Invalid or frozen target.')
    index=int(found[0]);q=np.zeros(h.shape[0]);q[index]=1
    basis=[];actions=[];trace=[];previous=None;last=None
    diag=np.r_[h.ints.energy[h.simple],h.dh[tuple(h.ip.T)],h.dp[tuple(h.ea.T)]]
    for iteration in range(1,args.cycles+1):
        basis.append(q);actions.append(h.matvec(q));b=np.array(basis).T;ab=np.array(actions).T
        projected=b.T@ab;values,u=np.linalg.eigh(projected)
        if previous is None:k=np.argmax(abs(b[index]@u)**2)
        else:k=np.argmax(abs(previous@b@u)**2)
        energy=values[k];x=b@u[:,k];r=ab@u[:,k]-energy*x
        row={'iteration':iteration,'energy':float(energy),
             'strength':float(np.linalg.norm(x[:h.ns])**2),'residual':float(np.linalg.norm(r)),
             'energy_change':None if last is None else float(energy-last),
             'last_projected_row':projected[-1].tolist()}
        trace.append(row);print(row,flush=True)
        if np.linalg.norm(r)<1e-10:break
        den=energy-diag;den=np.where(abs(den)<1e-5,1e-5,den);q=r/den
        for _ in range(2):q-=b@(b.T@q)
        q/=np.linalg.norm(q);last=energy;previous=x
    Path(args.output).write_text(json.dumps(trace,indent=2)+'\n')
if __name__=='__main__':main()
