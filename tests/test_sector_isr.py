"""Sector spectra, transition amplitudes and independent perturbation checks."""
from dataclasses import replace
from contextlib import nullcontext
import numpy as np
import pytest
from test_theory import ints
from nondiagonal_ept import SectorEPT, run_sector_methods, ConvergenceError
from nondiagonal_ept.blocks import Hamiltonian
from nondiagonal_ept.sector_nrl3 import canonical_reduction, nrl3_orders, NRL3SectorISR


@pytest.fixture(scope='module')
def mf():
    from pyscf import gto, scf
    return scf.RHF(gto.M(atom='O 0 0 0; H 0 -.7586 .5043; H 0 .7586 .5043',
                        basis='sto-3g', verbose=0)).run(conv_tol=1e-12)


@pytest.mark.parametrize('sector', ['ip', 'ea'])
@pytest.mark.parametrize('frozen', [1, [0, 6]])
def test_adc_against_direct_pyscf(mf, sector, frozen):
    from pyscf import adc
    direct = adc.ADC(mf, frozen=frozen)
    direct.method, direct.method_type = 'adc(3)', sector
    direct.approx_trans_moments = False
    direct.conv_tol, direct.tol_residual = 1e-12, 1e-8
    e, _, factors, amplitudes = direct.kernel(nroots=2)
    calc = SectorEPT(mf, 'nD-ADC(3)', frozen=frozen, sector=sector)
    poles = calc.kernel(nroots=2)
    np.testing.assert_allclose([p.energy for p in poles], -e if sector == 'ip' else e, atol=1e-10)
    np.testing.assert_allclose([p.strength for p in poles], factors/2, atol=2e-8)
    for j, p in enumerate(poles):
        # Overall signs of Dyson orbitals are arbitrary.
        np.testing.assert_allclose(np.outer(p.dyson_mo, p.dyson_mo),
                                   np.outer(amplitudes[:, j], amplitudes[:, j]), atol=2e-8)
        np.testing.assert_allclose(p.dyson_ao@mf.get_ovlp()@p.dyson_ao, p.strength, atol=1e-12)
        assert p.residual < 1e-9 and p.target is None
        assert p.binding_energy == -p.energy
    # Both spectra are unbound for attachment in this minimal-basis example.
    assert all(p.binding_energy > 0 if sector == 'ip' else p.binding_energy < 0 for p in poles)


def test_adc_failure_and_guards(mf):
    calc = SectorEPT(mf, frozen=1)
    with pytest.raises(ConvergenceError):
        calc.kernel(nroots=2, max_cycle=1)
    assert calc.results == []
    for n in (0, 1.5, calc.dimension+1):
        with pytest.raises(ValueError): calc.kernel(nroots=n)
    for kwargs in ({'tol': float('nan')}, {'tol': 0}, {'max_space': 1}):
        with pytest.raises(ValueError): calc.kernel(**kwargs)
    with pytest.raises(TypeError): calc.kernel(targets=[4])
    for opts in ({'sector':'ee'}, {'spin':2}, {'frozen':5}):
        with pytest.raises(ValueError): SectorEPT(mf, **opts)
    with pytest.raises(ValueError): run_sector_methods(mf, [])
    with pytest.raises(ValueError): run_sector_methods(mf, ['nD-ADC(3)', 'ND-ADC(3)'])


def test_nrl3_order_extraction(ints):
    ham = Hamiltonian(ints, 'NRL3')
    orders = nrl3_orders(ham, 1000)
    for coupling in (0., .2, .5, 1.):
        actual = Hamiltonian(replace(ints, spatial_eri=coupling*ints.spatial_eri), 'NRL3').dense()
        expected = sum(coupling**j*h for j, h in enumerate(orders))
        np.testing.assert_allclose(actual, expected, atol=1e-13)


def test_canonical_reduction_against_exact_exponentials():
    # Independent dense random example: no package molecular block equations.
    rng = np.random.default_rng(908)
    orders = [np.diag([-1.8, -.7, .3, 1.2])]
    for _ in range(3):
        a = rng.normal(size=(4, 4))*.1
        orders.append((a+a.T)/2)
    errors, cross_errors, moment_errors, resolvent_errors = [], [], [], []
    for lam in (.4, .2, .1):
        hs = [lam**j*h for j, h in enumerate(orders)]
        k, u, (s1, s2, s3), gap = canonical_reduction(hs, [0, 1])
        rotated = u.T@sum(hs)@u
        errors.append(np.linalg.norm(rotated-k))
        cross_errors.append(np.linalg.norm(rotated[:2, 2:]))
        # Independently expand exp(S) by perturbative degree through third order.
        polynomial = (np.eye(4)+s1+s2+s3+.5*(s1@s1+s1@s2+s2@s1)+s1@s1@s1/6)
        moment_errors.append(np.linalg.norm(u-polynomial))
        z = .1+1.2j
        resolvent_errors.append(np.linalg.norm(
            np.linalg.inv(z*np.eye(4)-sum(hs))-u@np.linalg.inv(z*np.eye(4)-k)@u.T))
        np.testing.assert_allclose(u.T@u, np.eye(4), atol=1e-13)
        np.testing.assert_allclose(k, k.T, atol=1e-13)
        np.testing.assert_allclose(k[:2, 2:], 0, atol=0)
        assert gap == 1.
    # Each discarded term and the full propagator error start at fourth order.
    for es in (errors, cross_errors, moment_errors, resolvent_errors):
        assert all(12 < es[i]/es[i+1] < 20 for i in range(2)), es


@pytest.mark.parametrize('sector', ['ip', 'ea'])
def test_nrl3_primary_poles_and_dyson_amplitudes_converge_at_order_four(ints, sector):
    energy_errors, amplitude_errors = [], []
    target = ints.nocc//2-1 if sector == 'ip' else ints.nocc//2
    for lam in (.5, .25, .125):
        scaled = replace(ints, spatial_eri=lam*ints.spatial_eri)
        full = Hamiltonian(scaled, 'NRL3')
        e, x = np.linalg.eigh(full.dense())
        root = np.argmax(x[target]**2)
        with pytest.warns(UserWarning, match='experimental'):
            reduced = NRL3SectorISR.from_integrals(scaled, sector=sector)
        pole = min(reduced.kernel(reduced.dimension), key=lambda p: abs(p.energy-e[root]))
        energy_errors.append(abs(pole.energy-e[root]))
        ref = x[:full.ns, root]
        amplitude_errors.append(np.linalg.norm(np.outer(pole.dyson_mo, pole.dyson_mo)-np.outer(ref, ref)))
    for es in (energy_errors, amplitude_errors):
        assert all(10 < es[i]/es[i+1] < 24 for i in range(2)), es


def test_sector_completeness_and_noninteracting_limit(ints):
    for scale in (0., 1.):
        scaled = replace(ints, spatial_eri=scale*ints.spatial_eri)
        spectrum = []
        for sector in ('ip', 'ea'):
            with pytest.warns(UserWarning):
                calc = NRL3SectorISR.from_integrals(scaled, sector=sector)
            spectrum.extend(calc.kernel(calc.dimension))
        # Combined one-spin orbital completeness, not an electron-number rule.
        moments = sum(np.outer(p.dyson_mo, p.dyson_mo) for p in spectrum)
        np.testing.assert_allclose(moments, np.eye(len(ints.coefficients)), atol=1e-12)
        assert all(-1e-13 <= p.strength <= 1+1e-13 for p in spectrum)
        if scale == 0:
            primary = [p for p in spectrum if p.strength > .5]
            np.testing.assert_allclose(sorted(p.energy for p in primary), ints.energy[::2], atol=1e-12)
            np.testing.assert_allclose([p.strength for p in primary], 1, atol=1e-12)


def test_molecular_sector_api_phase_spin_and_guards(mf):
    before = mf.mo_coeff.copy()
    with pytest.warns(UserWarning):
        results = run_sector_methods(mf, frozen=1, nroots=2)
    assert list(results) == ['nD-ADC(3)', 'NRL3-ISR(3)']
    flipped = mf.copy()
    flipped.mo_coeff = mf.mo_coeff*np.array([-1, 1, -1, 1, -1, 1, -1])
    for name, poles in results.items():
        with pytest.warns(UserWarning) if name == 'NRL3-ISR(3)' else nullcontext():
            changed = SectorEPT(flipped, name, frozen=1, spin=1).kernel(nroots=2)
        np.testing.assert_allclose([p.energy for p in poles], [p.energy for p in changed], atol=1e-10)
        for p, q in zip(poles, changed):
            np.testing.assert_allclose(np.outer(p.dyson_ao,p.dyson_ao), np.outer(q.dyson_ao,q.dyson_ao), atol=1e-9)
    np.testing.assert_array_equal(mf.mo_coeff, before)
    with pytest.raises(MemoryError): SectorEPT(mf, 'NRL3-ISR(3)', max_dimension=1)
    with pytest.raises(MemoryError): SectorEPT(mf, 'NRL3-ISR(3)', max_memory_mb=.5)
    h0 = np.diag([-1., 0., 0., 1.])
    with pytest.raises(ValueError, match='overlap'):
        canonical_reduction([h0, h0*0, h0*0, h0*0], [0, 1])
