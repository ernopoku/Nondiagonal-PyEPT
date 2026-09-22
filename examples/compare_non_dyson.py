"""Compare Dyson NRL3 with its static opposite-sector extension."""
import json
from pyscf import gto,scf
from nondiagonal_ept import EPT

cases=[('HF','H 0 0 0; F 0 0 .9168',1,[4,2],[5]),
       ('H2O','O 0 0 0; H 0 -.7586 .5043; H 0 .7586 .5043',1,[4,3],[5]),
       ('N2','N 0 0 0; N 0 0 1.1136',2,[6,4,3],[7])]
report=[]
for molecule,atom,frozen,ip,ea in cases:
    mf=scf.RHF(gto.M(atom=atom,unit='Angstrom',basis='cc-pvdz',verbose=0)).run(conv_tol=1e-12)
    for sector,targets in [('ip',ip),('ea',ea)]:
        methods={m:EPT(mf,m,frozen=frozen,sector=sector) for m in ['NRL3','nD-NRL3']}
        results={m:e.kernel(targets) for m,e in methods.items()}
        for i,target in enumerate(targets):
            a,b=(results[m][i] for m in methods)
            report.append(dict(molecule=molecule,atom=atom,basis='cc-pvdz',frozen=frozen,
                sector=sector,target=target,dyson_ev=a.binding_energy_ev,
                non_dyson_ev=b.binding_energy_ev,difference_ev=b.binding_energy_ev-a.binding_energy_ev,
                dyson_PS=a.strength,non_dyson_PS=b.strength,residual=b.residual,
                dimensions={m:e.hamiltonian.shape[0] for m,e in methods.items()},
                static_solve_residual=max(methods['nD-NRL3'].hamiltonian.static_solve_residuals)))
print(json.dumps(report,indent=2))
