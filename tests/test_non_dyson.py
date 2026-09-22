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
