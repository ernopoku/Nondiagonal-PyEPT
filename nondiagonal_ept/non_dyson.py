"""Static opposite-sector extension of NRL3 (not a separate ADC ISR)."""
from dataclasses import replace
import numpy as np
from scipy.sparse.linalg import LinearOperator, minres
from .blocks import Hamiltonian
from .solver import ConvergenceError


class StaticNRL3(Hamiltonian):
    """Freeze S_opp[p,q] as (S_opp(eps_p)[p,q]+S_opp(eps_q)[p,q])/2.

    Retain the full simple-orbital space and one triple manifold. The opposite
    manifold contributes to A once during construction, not during Davidson.
    """
    def __init__(self, ints, sector='ip', spin=0):
        super().__init__(ints, 'NRL3', sector, spin)
        ip_sector = sector == 'ip'
        opposite = self.bp if ip_sector else self.bh
        dimension = self.np if ip_sector else self.nh
        eps = ints.energy[self.simple]
        # Capture the original dimensions before reducing the secular space.
        nh, np_ = self.nh, self.np
        zeros_h, zeros_p = np.zeros(nh), np.zeros(np_)

        def opposite_action(x):
            h, p = self.triple_action(zeros_h, x) if ip_sector else self.triple_action(x, zeros_p)
            return p if ip_sector else h

        solutions = np.zeros_like(opposite)
        self.static_solve_residuals = np.zeros(self.ns)
        for row, energy in enumerate(eps):
            rhs = opposite[row]
            if dimension == 0 or np.linalg.norm(rhs) == 0:
                continue
            op = LinearOperator((dimension, dimension),
                                matvec=lambda x, e=energy: e*x-opposite_action(x),
                                dtype=float)
            solution, info = minres(op, rhs, rtol=1e-12,
                                   maxiter=max(200, 3*dimension))
            # MINRES stopping is relative to an estimated operator norm. Refine
            # against the actual residual so that ill-scaled systems cannot pass
            # merely because the solver reports success.
            for _ in range(2):
                error = rhs-op@solution
                if np.linalg.norm(error) <= 1e-10*max(1., np.linalg.norm(rhs)):
                    break
                correction, correction_info = minres(op, error, rtol=1e-13,
                                                     maxiter=max(200, 3*dimension))
                solution += correction
                if correction_info != 0:
                    info = correction_info
            residual = np.linalg.norm(op@solution-rhs)
            if info != 0 or not np.all(np.isfinite(solution)) or residual > 1e-10*max(1., np.linalg.norm(rhs)):
                raise ConvergenceError(
                    f'nD-NRL3 static resolvent failed for simple orbital {row}: '
                    f'info={info}, residual={residual:.3e}. '
                    'An HF energy may be near an opposite-sector pole.')
            solutions[row] = solution
            self.static_solve_residuals[row] = residual
        raw = solutions@opposite.T
        self.static_opposite = (raw+raw.T)/2
        self.a += self.static_opposite
        self.spec = replace(self.spec, name='nD-NRL3',
                            ip_interaction=ip_sector, ea_interaction=not ip_sector)
        if ip_sector:
            self.np = 0
            self.ea = np.empty((0,3), dtype=int)
            self.bp = np.zeros((self.ns,0))
            self.govov = None
        else:
            self.nh = 0
            self.ip = np.empty((0,3), dtype=int)
            self.bh = np.zeros((self.ns,0))
            self.goooo = self.gvoov = None
        self.shape = (self.ns+self.nh+self.np,)*2
        self._diag = None
