"""Schirmer--Angonoa Dyson-expansion static self-energy for Dyson ADC(3).

Dynamic ADC(3) vertices/blocks remain fixed at the canonical RHF reference.
Solve Sigma = W[Q + L(Sigma)], where Q = contour(G0 M G0), and L is the
occupied--virtual free-propagator response. No density trace rescaling is used.
See docs/DYSON_ADC.md for equations, conventions and validation scope.
"""
from copy import copy
from dataclasses import replace
import numpy as np
from scipy.sparse.linalg import LinearOperator, minres
from .solver import ConvergenceError


def coulomb_exchange(eri, density):
    """One-spin RHF density -> spin-diagonal static self-energy, 2J-K."""
    return (2*np.einsum('pqrs,rs->pq', eri, density, optimize=True)
            -np.einsum('prqs,rs->pq', eri, density, optimize=True))


def _triple_operator(ham, holes):
    # Avoid contracting the unused manifold, particularly its virtual ladder.
    side = copy(ham)
    if holes:
        side.np = 0
        side.ea = np.empty((0, 3), dtype=int)
        side.spec = replace(side.spec, ea_interaction=False)
        return lambda x: side.triple_action(x, np.empty(0))[0]
    side.nh = 0
    side.ip = np.empty((0, 3), dtype=int)
    side.spec = replace(side.spec, ip_interaction=False)
    return lambda x: side.triple_action(np.empty(0), x)[1]


def _checked_solve(action, rhs, tol, max_cycle, label):
    n = rhs.size
    if not n or np.linalg.norm(rhs) == 0:
        return np.zeros_like(rhs), 0.
    operator = LinearOperator((n, n), matvec=action, dtype=float)
    solution, info = minres(operator, rhs, rtol=tol*.01, maxiter=max_cycle)
    residual = float(np.linalg.norm(action(solution)-rhs))
    # Actual equation residual, not just the iterative solver's status flag.
    if (info != 0 or not np.all(np.isfinite(solution)) or not np.isfinite(residual)
            or residual > tol*max(1., np.linalg.norm(rhs))):
        raise ConvergenceError(f'ADC(3) {label} failed: info={info}, residual={residual:.3e}. '
                               'Check the reference/gaps or increase static_max_cycle.')
    return solution, residual


def dynamic_density(ham, *, tol=1e-10, max_cycle=500):
    """Contour integral Q of G0 M G0, in the active spatial MO basis.

    Xp_i = (eps_i - Dp)^(-1) Bp_i and Xh_a = (eps_a - Dh)^(-1) Bh_a.
    Only opposite-manifold resolvents are needed; no triple diagonalization.
    """
    no = ham.ints.nocc//2
    eps = ham.ints.energy[ham.simple]
    gap = eps[:no, None]-eps[None, no:]
    if np.any(gap >= -1e-10):
        raise ValueError('ADC(3) requires occupied HF energies below virtual energies with a nonzero gap.')
    holes, particles = _triple_operator(ham, True), _triple_operator(ham, False)
    xp = np.zeros((no, ham.np))
    xh = np.zeros((ham.ns-no, ham.nh))
    residuals = []
    for i, energy in enumerate(eps[:no]):
        xp[i], res = _checked_solve(lambda x, e=energy: e*x-particles(x),
                                   ham.bp[i], tol, max_cycle, f'2ph resolvent at MO {i}')
        residuals.append(res)
    for a, energy in enumerate(eps[no:]):
        xh[a], res = _checked_solve(lambda x, e=energy: e*x-holes(x),
                                   ham.bh[no+a], tol, max_cycle, f'2hp resolvent at MO {no+a}')
        residuals.append(res)
    q = np.zeros((ham.ns, ham.ns))
    q[:no, :no] = -xp@xp.T
    q[no:, no:] = xh@xh.T
    q[:no, no:] = (xp@ham.bp[no:].T+ham.bh[:no]@xh.T)/gap
    q[no:, :no] = q[:no, no:].T
    return q, residuals


def dem_static(ham, *, tol=1e-10, max_cycle=500):
    """Return the DEM static matrix and convergence diagnostics.

    The unknown t_ia = Sigma_ia/(eps_i-eps_a) solves a symmetric response
    equation. Once solved, all static blocks follow by one 2J-K contraction.
    Frozen orbitals carry their mean-field contribution only.
    """
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError('static_tol must be finite and positive.')
    if not isinstance(max_cycle, (int, np.integer)) or max_cycle < 1:
        raise ValueError('static_max_cycle must be a positive integer.')
    q, residuals = dynamic_density(ham, tol=tol, max_cycle=max_cycle)
    no, n = ham.ints.nocc//2, ham.ns
    eps = ham.ints.energy[ham.simple]
    gap = eps[:no, None]-eps[None, no:]
    eri = ham.ints.spatial_eri
    b = coulomb_exchange(eri, q)

    def response_density(t):
        rho = np.zeros((n, n))
        rho[:no, no:] = t.reshape(gap.shape)
        rho[no:, :no] = rho[:no, no:].T
        return rho

    def action(t):
        return (gap*t.reshape(gap.shape)
                -coulomb_exchange(eri, response_density(t))[:no, no:]).ravel()

    t, response_residual = _checked_solve(action, b[:no, no:].ravel(), tol,
                                         max_cycle, 'static response')
    density = q+response_density(t)
    sigma = coulomb_exchange(eri, density)
    # Evaluate the defining equation in ALL blocks, not only the solved block.
    check = sigma-coulomb_exchange(eri, q+response_density(sigma[:no, no:]/gap))
    residual = float(np.linalg.norm(check))
    if not np.all(np.isfinite(sigma)) or residual > tol*max(1., np.linalg.norm(sigma)):
        raise ConvergenceError(f'ADC(3) static fixed-point residual {residual:.3e} exceeds tolerance.')
    if np.max(np.abs(sigma-sigma.T)) > tol:
        raise ArithmeticError('ADC(3) static self-energy is not Hermitian.')
    return sigma, {
        'scheme': 'Schirmer-Angonoa DEM',
        'static_residual': residual,
        'response_residual': response_residual,
        'resolvent_residuals': residuals,
        'correlation_density_trace': float(2*np.trace(density)),
        'dynamic_density': q,
        'correlation_density': density,
    }
