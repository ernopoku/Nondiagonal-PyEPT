"""PySCF BCCD preparation for BD-T1; preserves the caller's RHF object."""
import numpy as np
from .integrals import from_pyscf
from .solver import ConvergenceError

def prepare(mf,frozen,max_memory_mb,tol):
    from pyscf import cc
    from pyscf.cc.bccd import bccd_kernel_
    from pyscf.cc.addons import spatial2spin
    # Validate reference before launching the more expensive BCCD iterations.
    from_pyscf(mf,frozen,max_memory_mb)
    copied=mf.copy()
    copied.mo_coeff=mf.mo_coeff.copy();copied.mo_energy=mf.mo_energy.copy()
    bd=cc.CCSD(copied,frozen=frozen)
    bd.conv_tol=min(1e-11,tol*.01);bd.conv_tol_normt=min(1e-9,tol*.1)
    bd.kernel()
    if not bd.converged:raise ConvergenceError('Initial CCSD did not converge for BCCD.')
    bd=bccd_kernel_(bd,conv_tol_normu=tol,max_cycle=50,canonicalization=True,verbose=0)
    # PySCF reinitializes the object during semicanonicalization: explicitly
    # reconverge and check both CC residuals and the singles norm afterward.
    bd.conv_tol=min(1e-11,tol*.01);bd.conv_tol_normt=min(1e-9,tol*.1)
    bd.kernel(t1=bd.t1,t2=bd.t2)
    if not bd.converged or np.linalg.norm(bd.t1)>max(tol*3,1e-8):
        raise ConvergenceError(f'Brueckner reference failed convergence; |t1|={np.linalg.norm(bd.t1):.3e}.')
    ref=bd._scf;ref.mo_coeff=bd.mo_coeff
    density=ref.make_rdm1();fock=ref.get_fock(dm=density)
    ref.mo_energy=np.diag(ref.mo_coeff.T@fock@ref.mo_coeff).copy()
    ints=from_pyscf(ref,frozen,max_memory_mb,_brueckner=True)
    return ints,spatial2spin(bd.t2),bd
