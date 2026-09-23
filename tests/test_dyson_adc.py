"""Independent contour, response-equation and FCI limits for Dyson ADC."""
from dataclasses import replace
import numpy as np
import pytest
from test_theory import ints
from oracle import physical_hamiltonian, operator_matrix
from nondiagonal_ept import EPT, run_methods, ConvergenceError
from nondiagonal_ept.blocks import Hamiltonian
from nondiagonal_ept.dyson_adc import dynamic_density, coulomb_exchange


def test_adc2_alias_is_existing_dyson_second_order(ints):
    for name in ('ADC(2)', 'Dyson-ADC(2)', 'ADC(2)-Dyson'):
        np.testing.assert_array_equal(Hamiltonian(ints,name).dense(), Hamiltonian(ints,'ND2').dense())


@pytest.mark.parametrize('method', ['3+', 'ND2'])
def test_dynamic_density_by_independent_complex_contour(ints,method):
    h = Hamiltonian(ints, method)
    full = h.dense()
    b = full[:h.ns, h.ns:]
    d = full[h.ns:, h.ns:]
    values, vectors = np.linalg.eigh(d)
    couplings = b@vectors
    eps = ints.energy[h.simple]
    negative = np.r_[eps[eps < 0], values[values < 0]]
    positive = np.r_[eps[eps > 0], values[values > 0]]
    left, right = min(negative)-.5, (max(negative)+min(positive))/2
    center, radius = (left+right)/2, (right-left)/2
    integral = np.zeros((h.ns, h.ns), dtype=complex)
    for theta in 2*np.pi*(np.arange(1024)+.5)/1024:
        dz = radius*np.exp(1j*theta)
        z = center+dz
        m = (couplings/(z-values))@couplings.T
        integral += dz*m/(z-eps[:,None])/(z-eps[None,:])/1024
    q, _ = dynamic_density(h)
    np.testing.assert_allclose(q, integral, atol=2e-12)


@pytest.mark.parametrize('method,base', [('ADC(3)','3+'), ('ADC(2)-DEM','ND2')])
def test_static_response_against_full_spin_linear_system(ints,method,base):
    h = Hamiltonian(ints, method)
    q = h.static_diagnostics['dynamic_density']
    n = ints.n
    g = ints.g('p','p','p','p')
    occ = np.arange(n) < ints.nocc
    eps = ints.energy
    response = np.zeros((n,n))
    unlike = occ[:,None] != occ[None,:]
    response[unlike] = ((occ[:,None].astype(float)-occ[None,:]) /
                       np.where(unlike,eps[:,None]-eps[None,:],1))[unlike]
    # Expand BOTH spin copies; no 2J-K code in this oracle.
    qspin = np.kron(q,np.eye(2))
    w = lambda rho: np.einsum('prqs,rs->pq',g,rho)
    eye = np.eye(n*n)
    coefficient = np.column_stack([(col.reshape(n,n)-w(response*col.reshape(n,n))).ravel()
                                  for col in eye])
    expected = np.linalg.solve(coefficient,w(qspin).ravel()).reshape(n,n)
    actual = h.a-np.diag(eps[h.simple])
    np.testing.assert_allclose(actual,expected[np.ix_(h.simple,h.simple)],atol=2e-12)
    np.testing.assert_allclose(coulomb_exchange(ints.spatial_eri,q),w(qspin)[::2,::2],atol=1e-13)
    strict = Hamiltonian(ints,base)
    np.testing.assert_array_equal(h.bh,strict.bh)
    np.testing.assert_array_equal(h.bp,strict.bp)
    np.testing.assert_array_equal(h.dense()[h.ns:,h.ns:],strict.dense()[h.ns:,h.ns:])
    assert np.linalg.norm(h.a-strict.a) > 1e-8
    assert h.static_diagnostics['static_residual'] < 1e-10


def test_static_is_correct_through_fourth_order_against_fci(ints):
    determinants = [i for i in range(2**ints.n) if i.bit_count() == ints.nocc]
    errors, differences = [], []
    for scale in (.4, .2, .1):
        scaled = replace(ints,spatial_eri=scale*ints.spatial_eri)
        physical = physical_hamiltonian(scaled)[np.ix_(determinants,determinants)]
        _, vectors = np.linalg.eigh(physical)
        psi = vectors[:,0]
        density = np.zeros((ints.n,ints.n))
        for r in range(ints.n):
            for s in range(ints.n):
                op = operator_matrix(ints.n,((r,1),(s,0)))[np.ix_(determinants,determinants)]
                density[r,s] = psi@op@psi
        density -= np.diag((np.arange(ints.n) < ints.nocc).astype(float))
        exact = np.einsum('prqs,rs->pq',scaled.g('p','p','p','p'),density)[::2,::2]
        adc = Hamiltonian(scaled,'ADC(3)')
        sigma = adc.a-np.diag(scaled.energy[adc.simple])
        errors.append(np.linalg.norm(sigma-exact))
        differences.append(np.linalg.norm(adc.a-Hamiltonian(scaled,'3+').a))
    # The DEM static remainder begins at fifth order; its change to 3+ at fourth.
    assert all(24 < errors[i]/errors[i+1] < 42 for i in range(2)), errors
    assert all(12 < differences[i]/differences[i+1] < 22 for i in range(2)), differences


@pytest.mark.parametrize('method', ['ADC(3)', 'ADC(2)-DEM'])
def test_noninteracting_limit_and_failed_resolvent(ints,monkeypatch,method):
    zero = replace(ints,spatial_eri=np.zeros_like(ints.spatial_eri))
    h = Hamiltonian(zero,method)
    np.testing.assert_array_equal(h.a,np.diag(zero.energy[h.simple]))
    assert h.static_diagnostics['static_residual'] == 0
    import nondiagonal_ept.dyson_adc as module
    monkeypatch.setattr(module,'minres',lambda op,rhs,**kw:(np.zeros_like(rhs),0))
    with pytest.raises(ConvergenceError,match='resolvent'):
        Hamiltonian(ints,method)


@pytest.fixture(scope='module')
def mf():
    from pyscf import gto, scf
    return scf.RHF(gto.M(atom='O 0 0 0; H 0 -.7586 .5043; H 0 .7586 .5043',
                        basis='sto-3g',verbose=0)).run(conv_tol=1e-12)


@pytest.mark.parametrize('sector,target', [('ip',4),('ea',5)])
@pytest.mark.parametrize('method', ['ADC(3)', 'ADC(2)-DEM'])
def test_molecular_dyson_adc_api(mf,sector,target,method):
    before = mf.mo_coeff.copy()
    ep = EPT(mf,method,frozen=[0,6],sector=sector)
    p = ep.kernel([target])[0]
    dense = ep.dense_spectrum()
    q = min(dense,key=lambda q:abs(q.energy-p.energy))
    np.testing.assert_allclose([p.energy,p.strength],[q.energy,q.strength],atol=1e-9)
    assert p.residual < 1e-9
    np.testing.assert_allclose(p.dyson_ao@mf.get_ovlp()@p.dyson_ao,p.strength,atol=1e-12)
    np.testing.assert_allclose(sum(q.strength for q in dense),ep.hamiltonian.ns,atol=1e-12)
    u = p.normalized_dyson_mo
    z = 1/(1-u@ep.self_energy(p.energy,derivative=True)@u)
    np.testing.assert_allclose(z,p.strength,atol=1e-9)
    changed = mf.copy()
    changed.mo_coeff = mf.mo_coeff*np.array([-1,1,-1,1,-1,1,-1])
    other = EPT(changed,method,frozen=[0,6],sector=sector,spin=1).kernel([target])[0]
    np.testing.assert_allclose([p.energy,p.strength],[other.energy,other.strength],atol=1e-10)
    np.testing.assert_array_equal(mf.mo_coeff,before)


def test_methods_and_static_controls(mf):
    result = run_methods(mf,['ADC(2)','ADC(3)','3+'],frozen=1,targets=[4])
    assert list(result) == ['ND2','ADC(3)','3+']
    for name in ['ADC(3)-DEM','Dyson-ADC(3)','ADC(3)-Dyson']:
        p = EPT(mf,name,frozen=1).kernel([4])[0]
        np.testing.assert_allclose(p.energy,result['ADC(3)'][0].energy,atol=1e-12)
    with pytest.raises(ConvergenceError): EPT(mf,'ADC(3)',frozen=1,static_max_cycle=1)
    for options in ({'static_tol':0}, {'static_tol':float('nan')}, {'static_max_cycle':0}):
        with pytest.raises(ValueError): EPT(mf,'ADC(3)',**options)
    with pytest.raises(ValueError): run_methods(mf,['ADC(3)','ADC(3)-DEM'])


def test_adc2_dem_leading_static_order_and_fixed_dynamic_blocks(ints):
    # The second-order density gives the complete third-order static term.
    from nondiagonal_ept.blocks import amplitudes, static_self_energy
    errors, sizes = [], []
    for scale in (.2, .1, .05):
        scaled = replace(ints, spatial_eri=scale*ints.spatial_eri)
        ham = Hamiltonian(scaled, 'ADC(2)-DEM')
        t, singles = amplitudes(scaled)
        third = static_self_energy(scaled,t,singles,'quadratic')[::2,::2]
        sigma = ham.a-np.diag(scaled.energy[ham.simple])
        errors.append(np.linalg.norm(sigma-third))
        sizes.append(np.linalg.norm(sigma))
    assert all(13 < errors[i]/errors[i+1] < 19 for i in range(2)), errors
    assert all(7 < sizes[i]/sizes[i+1] < 9 for i in range(2)), sizes


def test_adc2_dem_method_selection_and_controls(mf):
    result = run_methods(mf,['ND2','ADC(2)-DEM','ADC(3)-DEM'],frozen=1,targets=[4])
    assert list(result) == ['ND2','ADC(2)-DEM','ADC(3)']
    assert abs(result['ND2'][0].energy-result['ADC(2)-DEM'][0].energy) > 1e-6
    for options in ({'static_tol':0}, {'static_tol':float('nan')}, {'static_max_cycle':0}):
        with pytest.raises(ValueError): EPT(mf,'ADC(2)-DEM',**options)
    with pytest.raises(ConvergenceError): EPT(mf,'ADC(2)-DEM',static_max_cycle=1)
