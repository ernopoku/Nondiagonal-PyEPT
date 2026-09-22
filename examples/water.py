from pyscf import gto,scf,lib
from nondiagonal_ept import EPT
lib.num_threads(1)
mol=gto.M(atom='O 0 0 0; H 0 -.7586 .5043; H 0 .7586 .5043',basis='sto-3g',verbose=0)
mf=scf.RHF(mol).run(conv_tol=1e-12)
for method in ('ND2','2ph-TDA','NR2','NRP3','NRQ3','NRL3','3+'):
    ep=EPT(mf,method,frozen=1)
    for pole in ep.kernel(targets=[4,3,2]):
        print(f'{method:8s} MO {pole.target:2d} IP {pole.binding_energy_ev:12.6f} eV  Z {pole.strength:.8f} residual {pole.residual:.2e}')
