import numpy as np
import pytest
from pyscf import gto,scf,mp
from nondiagonal_ept import EPT,ConvergenceError
from nondiagonal_ept.blocks import amplitudes

@pytest.fixture(scope='module')
def mf():
    return scf.RHF(gto.M(atom='O 0 0 0; H 0 -.7586 .5043; H 0 .7586 .5043',basis='sto-3g',verbose=0)).run(conv_tol=1e-12)

@pytest.mark.parametrize('method',['ND2','2ph-TDA','NR2','NRP3','NRQ3','NRL3','3+'])
def test_molecular_poles_and_sum_rule(mf,method):
    ep=EPT(mf,method,frozen=1);spectrum=ep.dense_spectrum()
    assert abs(sum(p.strength for p in spectrum)-ep.hamiltonian.ns)<1e-10
    p=ep.kernel([4],tol=1e-10)[0]
    dense=min(spectrum,key=lambda q:abs(q.energy-p.energy))
    assert abs(dense.energy-p.energy)<1e-10
    assert abs(dense.strength-p.strength)<1e-8
    assert 0<=p.strength<=1
    np.testing.assert_allclose(np.linalg.norm(p.dyson_mo)**2,p.strength,atol=1e-13)
    assert p.residual<1e-10


def test_nd2_self_energy_and_derivative(mf):
    ep=EPT(mf,'ND2',frozen=1);h=ep.hamiltonian;e=-.3
    b=np.hstack([h.bh,h.bp]);diag=h.diagonal()[h.ns:]
    expected=(b/(e-diag))@b.T
    np.testing.assert_allclose(ep.self_energy(e),expected,atol=1e-10)
    analytic=ep.self_energy(e,derivative=True)
    numerical=(ep.self_energy(e+1e-5)-ep.self_energy(e-1e-5))/2e-5
    np.testing.assert_allclose(analytic,numerical,atol=1e-8)
    p=ep.kernel([4])[0];u=p.dyson_mo/np.sqrt(p.strength)
    z=1/(1-u@ep.self_energy(p.energy,derivative=True)@u)
    np.testing.assert_allclose(z,p.strength,atol=1e-9)


def test_mp2_conventions(mf):
    ep=EPT(mf,'ND2',frozen=1);t,_=amplitudes(ep.integrals)
    emp2=.25*np.einsum('ijab,ijab',t,ep.integrals.g('o','o','v','v'))
    reference=mp.MP2(mf,frozen=1).run().e_corr
    np.testing.assert_allclose(emp2,reference,atol=2e-12)


def test_guards(mf):
    with pytest.raises(ValueError):EPT(mf,'ADC(3)')
    with pytest.raises(ValueError):EPT(mf,'ND2',frozen=5)
    with pytest.raises(ValueError):EPT(mf,'ND2',frozen=1).kernel([0])
    with pytest.raises(ConvergenceError):EPT(mf,'NRL3',frozen=1).kernel([4],max_cycle=1)
    unconverged=mf.copy();unconverged.converged=False
    with pytest.raises(ValueError):EPT(unconverged)


def test_ea_and_mo_phase_invariance(mf):
    a=EPT(mf,'NRQ3',frozen=1,sector='ea')
    flipped=mf.copy();flipped.mo_coeff=mf.mo_coeff*np.array([-1,1,-1,1,-1,1,-1])
    b=EPT(flipped,'NRQ3',frozen=1,sector='ea')
    np.testing.assert_allclose([p.energy for p in a.dense_spectrum()], [p.energy for p in b.dense_spectrum()],atol=1e-12)
    p=a.kernel([5])[0]
    assert p.residual<1e-9


def test_brueckner_path_preserves_reference():
    mf=scf.RHF(gto.M(atom='H 0 0 0; H 0 0 .74',basis='6-31g',verbose=0)).run(conv_tol=1e-12)
    before=mf.mo_coeff.copy();ep=EPT(mf,'BD-T1');h=ep.hamiltonian;d=h.dense()
    np.testing.assert_allclose(mf.mo_coeff,before,atol=0)
    assert np.linalg.norm(ep.brueckner.t1)<3e-8
    np.testing.assert_allclose(d,d.T,atol=1e-12)
    np.testing.assert_allclose(np.diag(d),h.diagonal(),atol=1e-12)
    p=ep.kernel([0])[0];assert p.residual<1e-9


def test_multiple_methods_match_individual(mf):
    from nondiagonal_ept import run_methods
    methods = ['NRL3', 'NRQ3', 'NRP3']
    results = run_methods(mf, methods, frozen=1, targets=iter([4, 2]))
    assert list(results) == methods
    for method in methods:
        expected = EPT(mf, method, frozen=1).kernel([4, 2])
        assert [p.target for p in results[method]] == [4, 2]
        np.testing.assert_allclose([p.energy for p in results[method]],
                                   [p.energy for p in expected], atol=1e-12)
        assert all(p.residual < 1e-9 for p in results[method])
    for invalid in [[], ['ND2', 'ADC(2)']]:
        with pytest.raises(ValueError): run_methods(mf, invalid)
    with pytest.raises(TypeError): run_methods(mf, 'NRL3')
    with pytest.raises(ConvergenceError):
        run_methods(mf, ['NRL3'], frozen=1, targets=[4], max_cycle=1)
