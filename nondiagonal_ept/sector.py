"""Separate-sector solvers; root counts, unlike EPT's original MO targets."""

SECTOR_METHODS = ('nD-ADC(3)', 'NRL3-ISR(3)')


def sector_method(name):
    names = {key.upper(): key for key in SECTOR_METHODS}
    try:
        return names[name.upper().replace(' ', '')]
    except (KeyError, AttributeError):
        raise ValueError(f'Choose one of {SECTOR_METHODS}.') from None


def SectorEPT(mf, method='nD-ADC(3)', **options):
    """Construct a separate IP/EA solver. Call kernel(nroots=...), not targets.

    nD-ADC(3): established PySCF restricted ADC(3), PySCF >= 2.14 required.
    NRL3-ISR(3): experimental dense canonical reduction of the NRL3 matrix.
    """
    method = sector_method(method)
    if method == 'nD-ADC(3)':
        from .sector_adc import NonDysonADC3
        return NonDysonADC3(mf, **options)
    from .sector_nrl3 import NRL3SectorISR
    return NRL3SectorISR(mf, **options)


def run_sector_methods(mf, methods=SECTOR_METHODS, *, nroots=1, tol=1e-9,
                       max_cycle=150, max_space=40, **options):
    """Run sector methods sequentially on the same RHF reference."""
    if isinstance(methods, str):
        raise TypeError('methods must be a sequence of method names.')
    names = [sector_method(name) for name in methods]
    if not names or len(set(names)) != len(names):
        raise ValueError('Request at least one method, without duplicates.')
    return {name: SectorEPT(mf, name, **options).kernel(
        nroots=nroots, tol=tol, max_cycle=max_cycle, max_space=max_space)
        for name in names}
