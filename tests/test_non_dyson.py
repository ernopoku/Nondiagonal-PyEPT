"""Independent dense-resolvent checks for the static NRL3 extension."""
from dataclasses import replace
import numpy as np
import pytest
from test_theory import ints
from nondiagonal_ept.blocks import Hamiltonian
from nondiagonal_ept.non_dyson import StaticNRL3
from nondiagonal_ept.solver import self_energy


@pytest.mark.parametrize('sector', ['ip', 'ea'])
def test_static_resolvent_and_retained_spectrum(ints, sector):
    parent = Hamiltonian(ints, 'NRL3', sector)
    full = parent.dense()
    nd = StaticNRL3(ints, sector)
    simple = np.arange(parent.ns)
    holes = np.arange(parent.ns, parent.ns+parent.nh)
    particles = np.arange(parent.ns+parent.nh, parent.shape[0])
    retained, eliminated = (holes, particles) if sector == 'ip' else (particles, holes)
    b = full[np.ix_(simple, eliminated)]
    d = full[np.ix_(eliminated, eliminated)]
    eps = ints.energy[parent.simple]
    row_values = np.array([b[p]@np.linalg.solve(e*np.eye(len(d))-d,b.T)
                           for p,e in enumerate(eps)])
    static = (row_values+row_values.T)/2
    keep = np.r_[simple, retained]
    expected = full[np.ix_(keep,keep)].copy()
    expected[:parent.ns,:parent.ns] += static
    actual = nd.dense()
    np.testing.assert_allclose(nd.static_opposite, static, atol=1e-12)
    np.testing.assert_allclose(actual, expected, atol=1e-12)
    np.testing.assert_allclose(actual, actual.T, atol=1e-13)
    np.testing.assert_allclose(nd.diagonal(), np.diag(actual), atol=1e-13)
    assert nd.shape[0] == parent.shape[0]-len(eliminated)
    # Each diagonal exactly reproduces the frozen-at-own-HF-energy prescription.
    for p,e in enumerate(eps):
        np.testing.assert_allclose(static[p,p],b[p]@np.linalg.solve(e*np.eye(len(d))-d,b[p]),atol=1e-13)
    energy = -.45 if sector == 'ip' else .15
    derivative = self_energy(nd,energy,derivative=True)
    numerical = (self_energy(nd,energy+1e-5)-self_energy(nd,energy-1e-5))/2e-5
    np.testing.assert_allclose(derivative,numerical,atol=1e-8)
    # Removing the frozen sector's energy derivative changes the residues.
    assert np.linalg.eigvalsh(derivative).max() < 1e-10


@pytest.mark.parametrize('sector', ['ip', 'ea'])
def test_noninteracting_limit(ints, sector):
    zero = replace(ints, spatial_eri=np.zeros_like(ints.spatial_eri))
    nd = StaticNRL3(zero,sector)
    np.testing.assert_allclose(nd.static_opposite,0,atol=0)
    np.testing.assert_allclose(nd.a,np.diag(zero.energy[nd.simple]),atol=0)


def test_molecular_api_and_observables():
    from pyscf import gto,scf
    from nondiagonal_ept import EPT,run_methods
    mf=scf.RHF(gto.M(atom='O 0 0 0; H 0 -.7586 .5043; H 0 .7586 .5043',basis='sto-3g',verbose=0)).run(conv_tol=1e-12)
    for sector,target in [('ip',4),('ea',5)]:
        ep=EPT(mf,'nD-NRL3',frozen=1,sector=sector)
        import json
        json.dumps(ep.hamiltonian.static_diagnostics)
        p=ep.kernel([target])[0]
        dense=ep.dense_spectrum()
        q=min(dense,key=lambda q:abs(q.energy-p.energy))
        np.testing.assert_allclose([p.energy,p.strength],[q.energy,q.strength],atol=1e-9)
        assert p.residual<1e-9
        np.testing.assert_allclose(p.dyson_ao@mf.get_ovlp()@p.dyson_ao,p.strength,atol=1e-12)
        assert 0<=p.strength<=1
        np.testing.assert_allclose(sum(q.strength for q in dense),ep.hamiltonian.ns,atol=1e-10)
        beta=EPT(mf,'nD-NRL3',frozen=1,sector=sector,spin=1).kernel([target])[0]
        np.testing.assert_allclose(beta.energy,p.energy,atol=1e-10)
    results=run_methods(mf,['NRL3','nD-NRL3'],frozen=1,targets=[4])
    assert list(results)==['NRL3','nD-NRL3']
    with pytest.raises(ValueError):Hamiltonian(ep.integrals,'nD-NRL3')

@pytest.mark.parametrize('sector,index', [('ip',0),('ea',2)])
def test_weak_coupling_difference_starts_at_fourth_order(sector,index):
    from nondiagonal_ept.integrals import Integrals
    from nondiagonal_ept.solver import davidson
    rng=np.random.default_rng(31);n,no=4,2
    z=rng.normal(size=(n,n,n,n))*.08
    z=(z+z.swapaxes(0,1)+z.swapaxes(2,3)+z.swapaxes(0,1).swapaxes(2,3))/4
    z=(z+z.transpose(2,3,0,1))/2
    ints=Integrals(np.repeat([-1.4,-.8,.3,.9],2),z,np.repeat(np.arange(n),2),
                   np.tile([0,1],n),2*no,np.eye(n),np.arange(n),np.eye(n))
    errors=[]
    for scale in [.5,.25,.125]:
        scaled=replace(ints,spatial_eri=scale*z)
        dyson=Hamiltonian(scaled,'NRL3',sector)
        nd=StaticNRL3(scaled,sector)
        errors.append(abs(davidson(dyson,index,tol=1e-13)[0]-davidson(nd,index,tol=1e-13)[0]))
    # Fourth-order leading changes decrease by approximately 16 on halving V.
    assert all(12 < errors[i]/errors[i+1] < 22 for i in range(2))


def test_static_solver_failure_is_not_silently_accepted(ints,monkeypatch):
    import nondiagonal_ept.non_dyson as module
    from nondiagonal_ept import ConvergenceError
    monkeypatch.setattr(module,'minres',lambda op,rhs,**kw:(np.zeros_like(rhs),0))
    with pytest.raises(ConvergenceError,match='static resolvent failed'):
        StaticNRL3(ints,'ip')


@pytest.mark.parametrize('particles', [False, True])
@pytest.mark.parametrize('spin', [0, 1])
@pytest.mark.parametrize('cache_mb', [0, 2000])
def test_fast_opposite_action_against_parent(ints, particles, spin, cache_mb):
    from nondiagonal_ept.non_dyson import _opposite_operator
    parent = Hamiltonian(ints, 'NRL3', spin=spin)
    action, diagonal, cache = _opposite_operator(parent, particles, cache_mb)
    dim = parent.np if particles else parent.nh
    offset = parent.ns + (parent.nh if particles else 0)
    dense = parent.dense()[offset:offset+dim, offset:offset+dim]
    vector = np.random.default_rng(4).normal(size=dim)
    np.testing.assert_allclose(action(vector), dense@vector, atol=1e-13)
    np.testing.assert_allclose(diagonal, np.diag(dense), atol=1e-13)
    if particles and cache_mb:
        assert cache > 0
    else:
        assert cache == 0


def test_static_iteration_limit_with_good_residual_is_accepted(ints, monkeypatch):
    import nondiagonal_ept.non_dyson as module
    def exact_but_limit(op, rhs, **kwargs):
        dense = np.column_stack([op@x for x in np.eye(len(rhs))])
        return np.linalg.solve(dense, rhs), kwargs['maxiter']
    monkeypatch.setattr(module, 'minres', exact_but_limit)
    nd = StaticNRL3(ints)
    assert max(nd.static_solve_residuals) < 1e-12


@pytest.mark.parametrize('options', [dict(static_tol=0), dict(static_tol=float('nan')),
                                     dict(static_max_cycle=0), dict(static_max_cycle=1.5)])
def test_static_invalid_controls(ints, options):
    with pytest.raises(ValueError):
        StaticNRL3(ints, **options)


def test_indefinite_shift_uses_unchanged_resolvent(ints):
    # A high virtual HF sampling energy lies inside the opposite-sector
    # spectrum. Preconditioning must not replace or shift that resolvent.
    energies = ints.energy.copy()
    energies[-2:] = 3.1
    modified = replace(ints, energy=energies)
    parent = Hamiltonian(modified, 'NRL3')
    full = parent.dense()
    offset = parent.ns+parent.nh
    d = full[offset:,offset:]
    b = parent.bp
    eps = energies[parent.simple]
    assert np.min(np.linalg.eigvalsh(eps[-1]*np.eye(len(d))-d)) < 0
    assert np.max(np.linalg.eigvalsh(eps[-1]*np.eye(len(d))-d)) > 0
    raw = np.array([np.linalg.solve(e*np.eye(len(d))-d,b[p])@b.T
                    for p,e in enumerate(eps)])
    nd = StaticNRL3(modified)
    np.testing.assert_allclose(nd.static_opposite, (raw+raw.T)/2, atol=1e-11, rtol=1e-10)


def test_static_controls_and_progress_are_reported(ints, monkeypatch):
    import nondiagonal_ept.non_dyson as module
    calls = []
    original = module.minres
    def observed(op, rhs, **kwargs):
        calls.append(kwargs)
        return original(op, rhs, **kwargs)
    class Log:
        def __init__(self): self.lines = []
        def info(self, message, *args): self.lines.append(message % args)
    monkeypatch.setattr(module, 'minres', observed)
    log = Log()
    nd = StaticNRL3(ints, static_max_cycle=123, _log=log)
    assert calls and all(c['maxiter'] == 123 for c in calls)
    assert all(c['M'] is not None for c in calls)
    assert nd.static_diagnostics['matvecs'] > 0
    assert len(nd.static_diagnostics['iterations']) == nd.ns
    assert nd.static_diagnostics['setup_seconds'] >= nd.static_diagnostics['solve_seconds']
    assert len(log.lines) == nd.ns+1


@pytest.mark.parametrize('nv', [1, 3, 5])
@pytest.mark.parametrize('memory_mb', [0, 2000])
def test_pair_ladder_against_full_spatial_contraction(nv, memory_mb):
    from nondiagonal_ept.non_dyson import _pair_ladder
    rng = np.random.default_rng(83)
    eri = rng.normal(size=(nv,)*4)
    eri = (eri+eri.swapaxes(0,1)+eri.swapaxes(2,3)
           +eri.swapaxes(0,1).swapaxes(2,3))/4
    eri = (eri+eri.transpose(2,3,0,1))/2
    same = rng.normal(size=(2,nv,nv))
    same -= same.transpose(0,2,1)
    mixed = rng.normal(size=same.shape)
    apply, cache_bytes = _pair_ladder(eri, memory_mb)
    actual = apply(same,mixed)
    for tensor,result in zip((same,mixed),actual):
        expected = np.einsum('acbd,icd->iab',eri,tensor,optimize=False)
        np.testing.assert_allclose(result,expected,atol=1e-13,rtol=1e-13)


def test_compressed_action_roundoff_is_checked_against_parent(ints, monkeypatch):
    import nondiagonal_ept.non_dyson as module
    expected = StaticNRL3(ints).static_opposite
    original = module._opposite_operator
    def perturbed(*args):
        action, diagonal, cache = original(*args)
        # Exaggerate compression/roundoff error to exercise true-parent
        # residual refinement rather than checking the fast operator itself.
        return lambda x: action(x)+1e-7*x, diagonal, cache
    monkeypatch.setattr(module, '_opposite_operator', perturbed)
    nd = StaticNRL3(ints)
    np.testing.assert_allclose(nd.static_opposite,expected,atol=1e-11,rtol=1e-10)
    assert max(nd.static_solve_residuals) <= 1e-10
    assert nd.static_diagnostics['reference_matvecs'] > nd.ns
