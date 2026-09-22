"""Real canonical RHF -> active, antisymmetrized Dirac spin integrals."""
from dataclasses import dataclass
import numpy as np

@dataclass
class Integrals:
    energy: np.ndarray
    spatial_eri: np.ndarray
    spatial: np.ndarray
    spin: np.ndarray
    nocc: int
    coefficients: np.ndarray
    original_mos: np.ndarray
    overlap: np.ndarray
    fock: np.ndarray | None = None

    @property
    def n(self): return len(self.energy)
    @property
    def nvir(self): return self.n - self.nocc

    def indices(self, key):
        if isinstance(key, str):
            return {'o': np.arange(self.nocc), 'v': np.arange(self.nocc, self.n),
                    'p': np.arange(self.n)}[key]
        return np.atleast_1d(key).astype(int)

    def g(self, p, q, r, s):
        """<pq||rs> = (pr|qs) delta_spin - (ps|qr) delta_spin.

        All selectors retain their axes, including one-element selectors.
        Only requested blocks are expanded into spin orbitals.
        """
        p,q,r,s = [self.indices(k) for k in (p,q,r,s)]
        P,Q,R,S = [self.spatial[k] for k in (p,q,r,s)]
        sp,sq,sr,ss = [self.spin[k] for k in (p,q,r,s)]
        direct = self.spatial_eri[np.ix_(P,R,Q,S)].transpose(0,2,1,3)
        exchange = self.spatial_eri[np.ix_(P,S,Q,R)].transpose(0,2,3,1)
        direct *= (sp[:,None,None,None]==sr[None,None,:,None]) & (sq[None,:,None,None]==ss[None,None,None,:])
        exchange *= (sp[:,None,None,None]==ss[None,None,None,:]) & (sq[None,:,None,None]==sr[None,None,:,None])
        return direct-exchange


    def virtual_ladder(self, tensor):
        """1/2 <ab||cd> X_icd using antisymmetry and spatial ERIs.

        The exchange contraction equals the direct contraction for an
        antisymmetric X. Interleaved alpha/beta ordering is required.
        """
        occupied_spatial=self.nocc//2
        eri=self.spatial_eri[occupied_spatial:,occupied_spatial:,
                             occupied_spatial:,occupied_spatial:]
        out=np.empty_like(tensor)
        for sa in (0,1):
            for sb in (0,1):
                out[:,sa::2,sb::2]=np.einsum('acbd,icd->iab',eri,
                                            tensor[:,sa::2,sb::2],optimize=True)
        return out


def from_pyscf(mf, frozen=0, max_memory_mb=2000, *, _brueckner=False):
    """Transform a converged, real, canonical closed-shell RHF reference.

    frozen: count of lowest occupied spatial MOs, or explicit spatial MO indices
    (occupied or virtual). Frozen orbitals retain their mean-field contribution.
    """
    from pyscf import ao2mo, scf
    if not isinstance(mf, scf.hf.RHF) or isinstance(mf, scf.rohf.ROHF) or hasattr(mf, 'xc'):
        raise TypeError('A molecular closed-shell RHF reference is required; UHF/ROHF/DFT are not supported.')
    if not mf.converged: raise ValueError('SCF did not converge.')
    if getattr(mf,'with_df',None) is not None:
        raise TypeError('Use conventional RHF; density-fitted SCF references are not supported.')
    if max_memory_mb<=0:raise ValueError('max_memory_mb must be positive.')
    if not np.all(np.isfinite(mf.mo_coeff)) or not np.all(np.isfinite(mf.mo_energy)):
        raise ValueError('Non-finite reference orbitals or energies.')
    if np.iscomplexobj(mf.mo_coeff): raise ValueError('Complex orbitals are not supported.')
    occ = np.asarray(mf.mo_occ)
    if not np.all(np.isin(occ, [0,2])): raise ValueError('Integer closed-shell occupations are required.')
    occupied = np.flatnonzero(occ==2)
    if isinstance(frozen, (int,np.integer)):
        if frozen<0 or frozen>len(occupied): raise ValueError('Invalid frozen occupied count.')
        excluded = set(occupied[:frozen])
    else:
        excluded = set(frozen)
        if any(not isinstance(i,(int,np.integer)) or i<0 or i>=len(occ) for i in excluded):
            raise ValueError('Frozen indices must be valid zero-based spatial MO indices.')
    active = np.array([i for i in occupied if i not in excluded] +
                      [i for i in np.flatnonzero(occ==0) if i not in excluded])
    no = sum(occ[active]==2)
    if no==0 or no==len(active): raise ValueError('At least one active occupied and virtual spatial orbital are needed.')
    c = np.asarray(mf.mo_coeff[:,active])
    f = c.T @ mf.get_fock() @ c
    check=f - np.diag(mf.mo_energy[active])
    if _brueckner:
        check[:no,no:]=0;check[no:,:no]=0
    if np.max(np.abs(check))>2e-6:
        raise ValueError('Reference orbitals are not canonical RHF orbitals.')
    # Transformation workspace + restored spatial ERIs; no full spin n^4 array.
    nspin=2*len(active);ospin=2*no;vspin=nspin-ospin
    estimate = (4*len(active)**4 + 6*nspin*ospin*vspin*vspin
                + 4*nspin*ospin*ospin*vspin)*8/1e6
    if estimate>max_memory_mb:
        raise MemoryError(f'Estimated AO-to-MO workspace {estimate:.0f} MB exceeds max_memory_mb.')
    eri = ao2mo.restore(1, ao2mo.kernel(mf.mol,c,max_memory=max_memory_mb),len(active))
    return Integrals(np.repeat(mf.mo_energy[active],2), eri,
                     np.repeat(np.arange(len(active)),2), np.tile([0,1],len(active)),
                     2*no, c, active, mf.get_ovlp(), np.kron(f,np.eye(2)))
