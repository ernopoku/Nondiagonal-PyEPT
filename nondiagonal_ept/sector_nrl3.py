"""Experimental third-order canonical sector reduction of the NRL3 matrix.

This is an auxiliary-matrix ISR construction, NOT a separately derived
correlated-ground-state ADC formulation. See docs/SECTOR_ISR.md for the
derivation, the orthogonal transition-moment completion and limitations.
"""
import warnings
import numpy as np
from scipy.linalg import expm, eigh
from .blocks import Hamiltonian
from .integrals import from_pyscf
from .solver import Pole, ConvergenceError


def comm(a, b):
    return a@b-b@a


def canonical_reduction(orders, ip_indices, *, gap_tol=1e-8):
    """Block diagonalize H0+lambda H1+lambda^2 H2+lambda^3 H3 to order 3.

    H0 is diagonal. The generator S has only IP/EA off-block entries and
    U=exp(S); K=U.T H U. Return the strict cubic K and orthogonal U at lambda=1.
    """
    h0, h1, h2, h3 = orders
    n = len(h0)
    ip = np.asarray(ip_indices, dtype=int)
    ea = np.setdiff1d(np.arange(n), ip)
    if not np.isfinite(gap_tol) or gap_tol <= 0:
        raise ValueError('gap_tol must be finite and positive.')
    if (not len(ip) or not len(ea) or len(np.unique(ip)) != len(ip)
            or np.any(ip < 0) or np.any(ip >= n)):
        raise ValueError('IP indices must define a nonempty proper partition.')
    for h in orders:
        if h.shape != (n, n) or not np.all(np.isfinite(h)) or not np.allclose(h, h.T, atol=1e-12, rtol=0):
            raise ValueError('Order matrices must be real finite symmetric square matrices.')
    if np.max(np.abs(h0-np.diag(np.diag(h0)))) > 1e-12:
        raise ValueError('H0 must be diagonal.')
    eps = np.diag(h0)
    denominator = eps[ip, None]-eps[None, ea]
    gap = float(np.min(np.abs(denominator)))
    if gap < gap_tol:
        raise ValueError(f'IP/EA zeroth-order sectors overlap: gap {gap:.3e} < {gap_tol:.3e} Hartree.')

    def generator(r):
        s = np.zeros_like(h0)
        s[np.ix_(ip, ea)] = -r[np.ix_(ip, ea)]/denominator
        s[np.ix_(ea, ip)] = -s[np.ix_(ip, ea)].T
        return s

    s1 = generator(h1)
    c01 = comm(h0, s1)
    r2 = h2+comm(h1, s1)+.5*comm(c01, s1)
    s2 = generator(r2)
    c02 = comm(h0, s2)
    r3 = (h3+comm(h2, s1)+comm(h1, s2)
          +.5*(comm(c01, s2)+comm(c02, s1)+comm(comm(h1, s1), s1))
          +comm(comm(c01, s1), s1)/6)
    s3 = generator(r3)
    k = h0+h1+r2+r3
    k[np.ix_(ip, ea)] = 0
    k[np.ix_(ea, ip)] = 0
    if np.max(np.abs(k-k.T)) > 1e-10:
        raise ArithmeticError('Canonical sector matrix is not Hermitian.')
    # Exponential completion preserves positivity and combined spectral weight
    # exactly. It agrees with the cubic transition expansion through order 3.
    rotation = expm(s1+s2+s3)
    return (k+k.T)/2, rotation, (s1, s2, s3), gap


def nrl3_orders(ham, max_dimension):
    """Extract formal orders by scaling fluctuation integrals at fixed HF eps.

    Bare vertices and triple interactions are order 1; half-weight doubles
    vertices order 2; the linear singles static self-energy is order 3.
    """
    full = ham.dense(max_dimension)
    e0 = np.r_[ham.ints.energy[ham.simple], ham.dh[tuple(ham.ip.T)],
               ham.dp[tuple(ham.ea.T)]]
    h0 = np.diag(e0)
    h2 = np.zeros_like(full)
    bare_h = ham.ints.g('p', 'v', 'o', 'o').transpose(0, 2, 3, 1)
    bare_p = ham.ints.g('p', 'o', 'v', 'v')
    bare_h = bare_h[(ham.simple[:, None],)+tuple(ham.ip.T)]
    bare_p = bare_p[(ham.simple[:, None],)+tuple(ham.ea.T)]
    h2[:ham.ns, ham.ns:] = np.hstack([ham.bh-bare_h, ham.bp-bare_p])
    h2[ham.ns:, :ham.ns] = h2[:ham.ns, ham.ns:].T
    h3 = np.zeros_like(full)
    h3[:ham.ns, :ham.ns] = ham.a-np.diag(e0[:ham.ns])
    return h0, full-h0-h2-h3, h2, h3


class NRL3SectorISR:
    method = 'NRL3-ISR(3)'

    def __init__(self, mf, *, frozen=0, sector='ip', spin=0,
                 max_memory_mb=2000, max_dimension=1200, gap_tol=1e-8):
        if sector not in ('ip', 'ea'):
            raise ValueError("sector must be 'ip' or 'ea'.")
        if spin not in (0, 1):
            raise ValueError('spin must be 0 or 1.')
        self.integrals = from_pyscf(mf, frozen, max_memory_mb)
        self._build(self.integrals, sector, spin, max_memory_mb, max_dimension, gap_tol)

    @classmethod
    def from_integrals(cls, ints, *, sector='ip', spin=0, max_memory_mb=2000,
                       max_dimension=1200, gap_tol=1e-8):
        """Low-level constructor for perturbative validation."""
        obj = cls.__new__(cls)
        obj.integrals = ints
        obj._build(ints, sector, spin, max_memory_mb, max_dimension, gap_tol)
        return obj

    def _build(self, ints, sector, spin, memory, limit, gap_tol):
        no, nv = ints.nocc//2, ints.nvir//2
        dimension = no+nv+(no*(no-1)//2+no**2)*nv+no*(nv*(nv-1)//2+nv**2)
        if dimension > limit:
            raise MemoryError(f'NRL3-ISR(3) full auxiliary dimension {dimension} exceeds {limit}; '
                              'this experimental implementation is dense.')
        estimate = 48*dimension**2*8/1e6+ints.spatial_eri.nbytes/1e6
        if not np.isfinite(memory) or estimate > memory:
            raise MemoryError(f'Estimated dense ISR workspace {estimate:.0f} MB exceeds max_memory_mb.')
        ham = Hamiltonian(ints, 'NRL3', sector, spin)
        ip = np.r_[np.flatnonzero(ham.simple < ints.nocc),
                   np.arange(ham.ns, ham.ns+ham.nh)]
        ea = np.setdiff1d(np.arange(dimension), ip)
        orders = nrl3_orders(ham, limit)
        k, rotation, generators, gap = canonical_reduction(orders, ip, gap_tol=gap_tol)
        keep = ip if sector == 'ip' else ea
        self.matrix = k[np.ix_(keep, keep)]
        self.transition = rotation[:ham.ns, keep]
        self.sector, self.spin = sector, spin
        self.original_mos, self.coefficients = ints.original_mos, ints.coefficients
        self.overlap = ints.overlap
        self.dimension = len(keep)
        self.results = []
        rotated = rotation.T@sum(orders)@rotation
        self.diagnostics = {
            'full_dimension': int(dimension), 'sector_dimension': len(keep),
            'minimum_denominator_hartree': gap,
            'generator_frobenius_norms': [float(np.linalg.norm(s)) for s in generators],
            'remaining_cross_sector_norm_hartree': float(np.linalg.norm(rotated[np.ix_(ip, ea)])),
            'omitted_higher_order_norm_hartree': float(np.linalg.norm(rotated-k)),
        }
        warnings.warn('NRL3-ISR(3) is an experimental auxiliary-matrix sector reduction, '
                      'not a published or independently validated non-Dyson NRL3 formulation.',
                      UserWarning, stacklevel=3)

    def kernel(self, nroots=1, *, tol=1e-9, max_cycle=150, max_space=40):
        """IP roots in decreasing signed energy, EA roots in increasing energy.

        Dense eigensolve; max_cycle and max_space are accepted for the common
        sector API but do not control this solver. target=None, iterations=1.
        """
        self.results = []
        if not isinstance(nroots, (int, np.integer)) or not 1 <= nroots <= self.dimension:
            raise ValueError(f'nroots must be between 1 and {self.dimension}.')
        if not np.isfinite(tol) or tol <= 0:
            raise ValueError('tol must be finite and positive.')
        selected = [self.dimension-nroots, self.dimension-1] if self.sector == 'ip' else [0, nroots-1]
        energies, vectors = eigh(self.matrix, subset_by_index=selected)
        if self.sector == 'ip':
            energies, vectors = energies[::-1], vectors[:, ::-1]
        results = []
        for energy, vector in zip(energies, vectors.T):
            residual = float(np.linalg.norm(self.matrix@vector-energy*vector))
            if not np.isfinite(residual) or residual > tol:
                raise ConvergenceError(f'NRL3-ISR(3) residual {residual:.3e} exceeds {tol:.3e}.')
            mo = self.transition@vector
            if mo[np.argmax(np.abs(mo))] < 0:
                mo, vector = -mo, -vector
            ps = float(mo@mo)
            ao = self.coefficients@mo
            if not np.isclose(ao@self.overlap@ao, ps, atol=1e-9, rtol=1e-9):
                raise ArithmeticError('AO Dyson norm does not equal PS.')
            results.append(Pole(float(energy), ps, mo, ao, residual, 1, None,
                                self.sector, vector))
        self.results = results
        return results
