"""Static opposite-sector extension of NRL3 (not a separate ADC ISR)."""
from dataclasses import replace
from time import perf_counter
import numpy as np
from scipy.sparse.linalg import LinearOperator, minres
from .blocks import Hamiltonian
from .solver import ConvergenceError


class _ResidualConverged(Exception):
    """Internal early exit after checking the physical (unpreconditioned) residual."""
    def __init__(self, solution):
        self.solution = solution.copy()


class StaticNRL3(Hamiltonian):
    """Freeze S_opp[p,q] as (S_opp(eps_p)[p,q]+S_opp(eps_q)[p,q])/2.

    Retain the full simple-orbital space and one triple manifold. The opposite
    manifold contributes to A once during construction, not during Davidson.
    """
    def __init__(self, ints, sector='ip', spin=0, *, static_tol=1e-10,
                 static_max_cycle=5000, max_memory_mb=2000, _log=None):
        if not np.isfinite(static_tol) or static_tol <= 0:
            raise ValueError('static_tol must be finite and positive.')
        if not isinstance(static_max_cycle, (int, np.integer)) or static_max_cycle < 1:
            raise ValueError('static_max_cycle must be a positive integer.')
        if not np.isfinite(max_memory_mb) or max_memory_mb <= 0:
            raise ValueError('max_memory_mb must be finite and positive.')
        started = perf_counter()
        super().__init__(ints, 'NRL3', sector, spin)
        ip_sector = sector == 'ip'
        opposite = self.bp if ip_sector else self.bh
        dimension = self.np if ip_sector else self.nh
        eps = ints.energy[self.simple]
        # Apply only the eliminated block. The parent triple_action evaluates
        # both sectors, even when one input is identically zero.
        opposite_action, diagonal, cache_bytes = _opposite_operator(
            self, ip_sector, max_memory_mb)
        self.static_diagnostics = dict(solver='preconditioned MINRES',
                                       iterations=[], matvecs=0, reference_matvecs=0,
                                       ladder_cache_bytes=int(cache_bytes))
        self.static_solve_iterations = np.zeros(self.ns, dtype=int)

        zero_h, zero_p = np.zeros(self.nh), np.zeros(self.np)
        def reference_action(vector):
            self.static_diagnostics['matvecs'] += 1
            self.static_diagnostics['reference_matvecs'] += 1
            h, p = (self.triple_action(zero_h, vector) if ip_sector
                    else self.triple_action(vector, zero_p))
            return p if ip_sector else h

        solutions = np.zeros_like(opposite)
        self.static_solve_residuals = np.zeros(self.ns)
        solve_started = perf_counter()
        if _log is not None:
            _log.info('nD-NRL3: constructing static correction (%d orbital shifts, dimension %d)',
                      self.ns, dimension)
        for row, energy in enumerate(eps):
            rhs = opposite[row]
            if dimension == 0 or np.linalg.norm(rhs) == 0:
                continue
            def action(x, e=energy):
                self.static_diagnostics['matvecs'] += 1
                return e*x-opposite_action(x)
            op = LinearOperator((dimension, dimension), matvec=action, dtype=float)
            # MINRES requires a positive-definite preconditioner, including
            # when the shifted auxiliary matrix itself is indefinite.
            gap = np.abs(energy-diagonal)
            # Avoid a nearly zero diagonal creating enormous weights in an
            # indefinite system. This changes only the preconditioner, never
            # the actual shifted matrix or its required residual.
            floor = max(1e-8, .01*np.max(gap))
            inverse = 1/np.maximum(gap, floor)
            preconditioner = LinearOperator(op.shape, matvec=lambda x: inverse*x,
                                           dtype=float)
            threshold = static_tol*max(1., np.linalg.norm(rhs))
            def solve(right_hand_side, rtol):
                iterations = 0
                def check_iteration(x):
                    nonlocal iterations
                    iterations += 1
                    self.static_solve_iterations[row] += 1
                    # Refinement RHSs can be tiny. A fixed relative tolerance
                    # otherwise oversolves them by many orders of magnitude.
                    # Test the actual residual in the original units instead.
                    if iterations % 8 == 0:
                        error = right_hand_side-op@x
                        if np.linalg.norm(error) <= .5*threshold:
                            raise _ResidualConverged(x)
                try:
                    return minres(op, right_hand_side, rtol=rtol,
                                  M=preconditioner, maxiter=static_max_cycle,
                                  callback=check_iteration)
                except _ResidualConverged as done:
                    return done.solution, 0
            solution, info = solve(rhs, min(1e-12, static_tol*.01))
            # MINRES stopping is relative to an estimated operator norm. Refine
            # against the actual residual so that ill-scaled systems cannot pass
            # merely because the solver reports success.
            for _ in range(2):
                error = rhs-(energy*solution-reference_action(solution))
                if np.linalg.norm(error) <= static_tol*max(1., np.linalg.norm(rhs)):
                    break
                correction, correction_info = solve(error, min(1e-13, static_tol*.001))
                solution += correction
                if correction_info != 0:
                    info = correction_info
            # Certify with the original unpacked parent action. Pair symmetry
            # identities can differ at roundoff level for transformed ERIs;
            # any resulting discrepancy is repaired by the refinement above.
            residual = np.linalg.norm(energy*solution-reference_action(solution)-rhs)
            if info < 0 or not np.all(np.isfinite(solution)) or not np.isfinite(residual) or residual > static_tol*max(1., np.linalg.norm(rhs)):
                raise ConvergenceError(
                    f'nD-NRL3 static resolvent failed for simple orbital {row}: '
                    f'info={info}, residual={residual:.3e}. '
                    'Increase static_max_cycle if the iteration limit was reached; '
                    'an HF energy may also be near an opposite-sector pole.')
            solutions[row] = solution
            self.static_solve_residuals[row] = residual
            if _log is not None:
                _log.info('  static orbital %d/%d: %d iterations, residual %.3e',
                          row+1, self.ns, self.static_solve_iterations[row], residual)
        self.static_diagnostics.update(
            iterations=self.static_solve_iterations.tolist(),
            residuals=self.static_solve_residuals.tolist(),
            solve_seconds=perf_counter()-solve_started)
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
        self.static_diagnostics['setup_seconds'] = perf_counter()-started


def _opposite_operator(ham, particles, max_memory_mb):
    """Return an exact single-sector action, its diagonal, and cache size.

    Cache the spatial virtual ladder in BLAS layout when it fits a small
    fraction of the memory budget. Never construct a dense auxiliary matrix
    or a spin-orbital four-virtual integral tensor.
    """
    ints = ham.ints
    o, v = ints.nocc, ints.nvir
    if not particles:
        i, j, a = ham.ip.T
        diagonal = (ham.dh[i,j,a] - ham.goooo[i,j,i,j]
                    + ham.gvoov[a,i,a,i] + ham.gvoov[a,j,a,j])
        def action(vector):
            x = np.zeros((o,o,v))
            x[i,j,a] = vector
            x[j,i,a] = -vector
            z = ham.dh*x - .5*np.einsum('ijkl,kla->ija',ham.goooo,x,optimize=True)
            ring = np.einsum('biak,kjb->ija',ham.gvoov,x,optimize=True)
            z += ring-ring.swapaxes(0,1)
            return z[i,j,a]
        return action, diagonal, 0
    i, a, b = ham.ea.T
    aa, bb = ints.spatial[o+a], ints.spatial[o+b]
    same = ints.spin[o+a] == ints.spin[o+b]
    eri = ints.spatial_eri
    diagonal = (ham.dp[i,a,b] + eri[aa,aa,bb,bb]
                - same*eri[aa,bb,bb,aa]
                - ham.govov[i,a,i,a] - ham.govov[i,b,i,b])
    spatial = eri[o//2:,o//2:,o//2:,o//2:]
    ladder, cache_bytes = _pair_ladder(spatial, max_memory_mb)
    spin = ham.spin
    other = 1-spin
    # Only two independent spin blocks can be nonzero in this conserved
    # creator-spin sector. The mixed block determines its transpose partner.
    selections = ((slice(spin,None,2),slice(spin,None,2),slice(spin,None,2)),
                  (slice(other,None,2),slice(spin,None,2),slice(other,None,2)))
    # Ring contractions conserve spin(i)-spin(a). Pack their nonzero pair
    # blocks once, rather than multiplying the full mostly-zero spin tensor.
    ring_matrix = ham.govov.transpose(2,1,0,3).reshape(o*v,o*v)
    pair_i, pair_a = np.indices((o,v))
    pair_spin = ints.spin[pair_i]-ints.spin[o+pair_a]
    ring_blocks = []
    for sb in (spin, other):
        pairs = np.flatnonzero(pair_spin.ravel() == sb-spin)
        columns = np.arange(sb,v,2)
        ring_blocks.append((pairs,columns,ring_matrix[np.ix_(pairs,pairs)]))
    def action(vector):
        y = np.zeros((o,v,v))
        y[i,a,b] = vector
        y[i,b,a] = -vector
        w = ham.dp*y
        values = ladder(y[selections[0]], y[selections[1]])
        w[selections[0]] += values[0]
        w[selections[1]] += values[1]
        w[other::2,other::2,spin::2] -= values[1].transpose(0,2,1)
        flat = y.reshape(o*v,v)
        ring = np.zeros_like(flat)
        for pairs,columns,matrix in ring_blocks:
            ring[np.ix_(pairs,columns)] = matrix@flat[np.ix_(pairs,columns)]
        ring = ring.reshape(o,v,v)
        w -= ring-ring.swapaxes(1,2)
        return w[i,a,b]
    return action, diagonal, cache_bytes



def _pair_ladder(spatial, max_memory_mb):
    """Exact spatial ladder in normalized symmetric/antisymmetric pair bases.

    Pair exchange commutes with (ac|bd), so these two blocks do not couple.
    The same-spin tensor is antisymmetric; the mixed-spin tensor uses both.
    The transformations here are internal and do not change the external
    unique-configuration normalization used by Hamiltonian.
    """
    nv = spatial.shape[0]
    aa, ab = np.triu_indices(nv, 1)
    sa, sb = np.triu_indices(nv)
    symmetric_weights = np.where(sa == sb, 1/np.sqrt(2.), 1.)
    cache_bytes = (len(aa)**2 + len(sa)**2)*spatial.dtype.itemsize
    limit = min(1024e6, .1*max_memory_mb*1e6)

    def matrix(a, b, first, last, symmetric):
        x, y = a[first:last,None], b[first:last,None]
        c, d = a[None,:], b[None,:]
        direct = spatial[x,c,y,d]
        exchange = spatial[x,d,y,c]
        if symmetric:
            direct += exchange
            direct *= symmetric_weights[first:last,None]
            direct *= symmetric_weights[None,:]
        else:
            direct -= exchange
        return direct

    cached = cache_bytes <= limit
    anti_matrix = matrix(aa,ab,0,len(aa),False) if cached else None
    sym_matrix = matrix(sa,sb,0,len(sa),True) if cached else None

    def multiply(rows, a, b, packed, symmetric):
        if packed is not None:
            return rows@packed.T
        if len(a) == 0:
            return np.empty_like(rows)
        # Three slab-sized arrays can be live during kernel construction.
        slab = max(1, min(len(a), int(limit/(3*len(a)*spatial.dtype.itemsize))))
        result = np.empty_like(rows)
        for first in range(0,len(a),slab):
            last = min(first+slab,len(a))
            result[:,first:last] = rows@matrix(a,b,first,last,symmetric).T
        return result

    inverse_sqrt2 = 1/np.sqrt(2.)
    def apply(same_spin, mixed_spin):
        occupied = same_spin.shape[0]
        same = np.sqrt(2.)*same_spin[:,aa,ab]
        mixed_anti = (mixed_spin[:,aa,ab]-mixed_spin[:,ab,aa])*inverse_sqrt2
        mixed_sym = ((mixed_spin[:,sa,sb]+mixed_spin[:,sb,sa])
                     *np.where(sa==sb,.5,inverse_sqrt2))
        anti = multiply(np.concatenate((same,mixed_anti)),aa,ab,anti_matrix,False)
        sym = multiply(mixed_sym,sa,sb,sym_matrix,True)
        result_same = np.zeros_like(same_spin)
        result_same[:,aa,ab] = anti[:occupied]*inverse_sqrt2
        result_same[:,ab,aa] = -anti[:occupied]*inverse_sqrt2
        result_mixed = np.zeros_like(mixed_spin)
        symmetric_part = sym*np.where(sa==sb,1.,inverse_sqrt2)
        result_mixed[:,sa,sb] = symmetric_part
        result_mixed[:,sb,sa] = symmetric_part
        result_mixed[:,aa,ab] += anti[occupied:]*inverse_sqrt2
        result_mixed[:,ab,aa] -= anti[occupied:]*inverse_sqrt2
        return result_same, result_mixed
    return apply, cache_bytes if cached else 0
