"""Non-diagonal electron propagators using a PySCF RHF reference."""
from .integrals import from_pyscf
from .blocks import Hamiltonian
from .methods import METHODS,method_spec
from .solver import davidson,make_pole,self_energy,ConvergenceError,HARTREE_TO_EV
from .sector import SectorEPT, run_sector_methods, SECTOR_METHODS
import numpy as np

class EPT:
    def __init__(self,mf,method='NRL3',*,frozen=0,sector='ip',spin=0,max_memory_mb=2000,brueckner_tol=1e-8,
                 static_tol=1e-10,static_max_cycle=None,static_space="sector"):
        spec=method_spec(method,sector)
        if spec.name != 'nD-NRL3' and static_space != 'sector':
            raise ValueError('static_space is only supported for nD-NRL3.')
        if static_max_cycle is None:
            static_max_cycle = 5000 if spec.name == 'nD-NRL3' else 500
        if spec.name.upper()=='BD-T1':
            from .brueckner import prepare
            self.integrals,t,self.brueckner=prepare(mf,frozen,max_memory_mb,brueckner_tol)
            self.hamiltonian=Hamiltonian(self.integrals,method,sector,spin,doubles=t)
        elif spec.name == 'nD-NRL3':
            from .non_dyson import StaticNRL3
            from pyscf.lib import logger
            self.integrals=from_pyscf(mf,frozen,max_memory_mb)
            self.hamiltonian=StaticNRL3(self.integrals,sector,spin,
                                           static_tol=static_tol,
                                           static_max_cycle=static_max_cycle,
                                           max_memory_mb=max_memory_mb, static_space=static_space,
                                           _log=logger.new_logger(mf))
        else:
            self.integrals=from_pyscf(mf,frozen,max_memory_mb)
            self.hamiltonian=Hamiltonian(self.integrals,method,sector,spin,
                                         static_tol=static_tol,static_max_cycle=static_max_cycle)
        self.method=self.hamiltonian.spec.name
        self.results=[]

    def kernel(self,targets=None,*,tol=1e-9,max_cycle=150,max_space=40):
        """Primary poles for zero-based ORIGINAL spatial MO indices.

        Default: active occupied orbitals for IP, active virtual orbitals for EA.
        nD-NRL3 defaults to sector-specific simple and triple spaces; use
        static_space="full" explicitly for the legacy full simple space.
        """
        h=self.hamiltonian
        original=self.integrals.original_mos[self.integrals.spatial[h.simple]]
        if targets is None:
            choose=h.simple<self.integrals.nocc
            if h.sector=='ea':choose=~choose
            targets=original[choose]
        results=[]
        for target in targets:
            found=np.flatnonzero(original==target)
            if not len(found):raise ValueError(f"MO {target} is frozen, invalid, or outside this method's {h.sector.upper()} simple space.")
            omega,x,res,it=davidson(h,int(found[0]),tol,max_cycle,max_space)
            if any(abs(p.vector@x)>1-1e-7 for p in results):
                raise ConvergenceError('Two targets converged to the same pole; use dense_spectrum to resolve strongly mixed states.')
            results.append(make_pole(h,omega,x,res,it,int(target)))
        self.results=results
        return results

    def dense_spectrum(self,max_dimension=2500):
        """All poles, including satellites; intended for small-system verification."""
        h=self.hamiltonian;matrix=h.dense(max_dimension)
        if np.max(np.abs(matrix-matrix.T))>1e-9:raise ArithmeticError('Non-Hermitian Hamiltonian.')
        energies,vectors=np.linalg.eigh(matrix)
        return [make_pole(h,e,x,np.linalg.norm(matrix@x-e*x),1) for e,x in zip(energies,vectors.T)]

    def self_energy(self,energy,derivative=False,**kwargs):
        return self_energy(self.hamiltonian,energy,derivative,**kwargs)

def run_methods(mf, methods, *, targets=None, tol=1e-9, max_cycle=150,
                max_space=40, **options):
    """Run several methods on one converged SCF reference.

    Return an insertion-ordered dictionary mapping canonical method names to
    lists of Pole objects. EPT options (frozen, sector, spin, etc.) are shared.
    Methods run sequentially; each builds its own propagator intermediates.
    No SCF calculation is repeated. Any convergence failure raises immediately.
    """
    if isinstance(methods, str):
        raise TypeError("methods must be a sequence, e.g. ['NRL3', 'NRQ3']")
    names = [method_spec(name, options.get('sector', 'ip')).name for name in methods]
    if not names:
        raise ValueError('Request at least one method.')
    if len(set(names)) != len(names):
        raise ValueError('Duplicate methods or aliases are not allowed.')
    targets = None if targets is None else list(targets)
    results = {}
    for name in names:
        calculation = EPT(mf, name, **options)
        results[name] = calculation.kernel(targets, tol=tol, max_cycle=max_cycle,
                                           max_space=max_space)
    return results


__all__=['run_methods','EPT','Hamiltonian','METHODS','ConvergenceError','HARTREE_TO_EV',
         'SectorEPT','run_sector_methods','SECTOR_METHODS']
