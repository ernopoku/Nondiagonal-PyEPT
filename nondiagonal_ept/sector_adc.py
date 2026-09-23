"""Restricted non-Dyson ADC(3), using PySCF's IP/EA implementation."""
import re
import numpy as np
from .integrals import reference_space
from .solver import Pole, ConvergenceError


class NonDysonADC3:
    method = 'nD-ADC(3)'

    def __init__(self, mf, *, frozen=0, sector='ip', spin=0,
                 max_memory_mb=2000):
        import pyscf
        from pyscf import adc
        if sector not in ('ip', 'ea'):
            raise ValueError("sector must be 'ip' or 'ea'.")
        if spin not in (0, 1):
            raise ValueError('spin must be 0 or 1.')
        # Frozen-space/transition-moment behavior is verified with this API.
        version = tuple(map(int, re.findall(r'\d+', pyscf.__version__)[:2]))
        if version < (2, 14):
            raise ImportError('nD-ADC(3) requires PySCF >= 2.14.')
        active, c, _, no = reference_space(mf, frozen, max_memory_mb)
        if np.any(np.diff(active) < 0):
            raise ValueError('PySCF ADC requires occupied MOs before virtual MOs.')
        excluded = sorted(set(range(len(mf.mo_occ)))-set(active))
        self.backend = adc.ADC(mf, frozen=excluded or None)
        self.backend.method = 'adc(3)'
        self.backend.method_type = sector
        self.backend.compute_properties = True
        self.backend.approx_trans_moments = False
        self.backend.max_memory = max_memory_mb
        self.backend.verbose = 0  # Package examples print PS, using per-spin factors.
        self.backend.if_heri_eris = True
        setattr(self.backend, sector+'_adc', self._solve)
        self.sector, self.spin = sector, spin
        self.original_mos, self.coefficients = active, c
        self.overlap = mf.get_ovlp()
        nv = len(active)-no
        self.dimension = no+no*no*nv if sector == 'ip' else nv+no*nv*nv
        self.results = []

    def _solve(self, nroots=1, guess=None, eris=None):
        from pyscf import lib
        from pyscf.adc import radc_ip, radc_ea
        es = (radc_ip.RADCIP if self.sector == 'ip' else radc_ea.RADCEA)(self.backend)
        action, diagonal = es.gen_matvec(es.get_imds(eris), eris)
        guess = es.get_init_guess(nroots, diagonal, ascending=True)
        # PySCF's default linear-dependence cutoff may discard corrections
        # before a requested sub-1e-7 residual is attained. Lower it explicitly.
        converged, es.E, vectors = lib.linalg_helper.davidson_nosym1(
            lambda xs: [action(x) for x in xs], guess, diagonal,
            nroots=nroots, tol=es.conv_tol, tol_residual=es.tol_residual,
            lindep=min(1e-24, es.tol_residual**2*.01),
            max_cycle=es.max_cycle, max_space=es.max_space,
            max_memory=es.max_memory, verbose=0)
        if not np.all(converged):
            raise ConvergenceError('ADC Davidson did not converge all requested roots.')
        es.U = np.asarray(vectors).T.copy()
        es.P, es.X = es.get_properties(nroots)
        self._action = action
        return es.E, es.U, es.P, es.X, es

    def kernel(self, nroots=1, *, tol=1e-9, max_cycle=150, max_space=40):
        """Lowest IP costs or lowest attachment energies; nroots is not MO indices.

        PySCF spin-summed factors are divided by two. Its transition amplitudes
        already have the desired one-spin normalization and are not rescaled.
        iterations=0 means the backend does not expose an iteration count.
        """
        self.results = []
        if not isinstance(nroots, (int, np.integer)) or not 1 <= nroots <= self.dimension:
            raise ValueError(f'nroots must be between 1 and {self.dimension}.')
        if not np.isfinite(tol) or tol <= 0 or max_cycle < 1 or max_space < max(4, nroots):
            raise ValueError('Require finite tol > 0, max_cycle >= 1 and max_space >= max(4,nroots).')
        adc = self.backend
        adc.conv_tol = min(tol, 1e-10)
        adc.tol_residual = tol*.1
        adc.max_cycle, adc.max_space = max_cycle, max_space
        energy, vectors, factors, moments, eris = adc.kernel(nroots=nroots)
        action = self._action
        if len(energy) != nroots:
            raise ConvergenceError('ADC returned fewer roots than requested.')
        results = []
        for e, x, factor, mo in zip(energy, vectors.T, factors, moments.T):
            if any(np.iscomplexobj(a) for a in (e, x, mo)):
                raise ConvergenceError('ADC returned a complex eigenpair.')
            residual = np.linalg.norm(action(x)-e*x)/np.linalg.norm(x)
            if not np.isfinite(residual) or residual > tol:
                raise ConvergenceError(f'ADC eigenpair residual {residual:.3e} exceeds {tol:.3e}.')
            mo, x = mo.copy(), x.copy()
            if mo[np.argmax(np.abs(mo))] < 0:
                mo, x = -mo, -x
            ps = float(mo@mo)
            ao = self.coefficients@mo
            if not np.isfinite(ps) or not np.all(np.isfinite(ao)):
                raise ArithmeticError('Non-finite ADC transition moments.')
            if not np.isclose(ps, factor/2, atol=1e-10, rtol=1e-9):
                raise ArithmeticError('Unexpected PySCF ADC spectroscopic-factor convention.')
            if not np.isclose(ao@self.overlap@ao, ps, atol=1e-9, rtol=1e-9):
                raise ArithmeticError('AO Dyson norm does not equal PS.')
            omega = -float(e) if self.sector == 'ip' else float(e)
            results.append(Pole(omega, ps, mo, ao, float(residual), 0,
                                None, self.sector, x))
        self.results = results
        return results
