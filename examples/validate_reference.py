"""Compare reference poles; prints a validation report.

python examples/validate_reference.py --molecules HF H2O N2 --output report.json
The molecular geometries, spherical cc-pVTZ basis, and frozen cores match the
provided .com files. Target indices refer to the original PySCF MOs.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from pyscf import gto,scf
from nondiagonal_ept import EPT

GEOMETRIES={
 'HF':'H 0 0 0; F 0 0 .9178',
 'N2':'N 0 0 0; N 0 0 1.1136',
 'F2':'F 0 0 0; F 0 0 1.4057',
 'CO':'C 0 0 0; O 0 0 1.1372',
 'H2O': [('O',(0,0,0)),('H',(0,0,.9606)),('H',(.9606*np.sin(np.deg2rad(103.4614)),0,.9606*np.cos(np.deg2rad(103.4614))))],
}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--molecules',nargs='+',default=['HF'],choices=list(GEOMETRIES))
 parser.add_argument('--output',default='reference_comparison.json');args=parser.parse_args()
 data=json.loads((Path(__file__).resolve().parents[1]/'tests/reference_poles.json').read_text())
 report=[]
 for name in args.molecules:
  mf=scf.RHF(gto.M(atom=GEOMETRIES[name],basis='cc-pvtz',verbose=0)).run(conv_tol=1e-13)
  for method,sector in [('NRP3','ip'),('NRQ3','ip'),('NRL3','ip'),('NRL3','ea'),('NRP3','ea'),('NRQ3','ea')]:
   refs=[r for r in data if r['molecule']==name and r['method']==method and r['sector']==sector]
   start=time.time();ep=EPT(mf,method,sector=sector,frozen=refs[0]['frozen'])
   poles=ep.kernel(targets=[r['original_mo'] for r in refs],tol=2e-9)
   for ref,p in zip(refs,poles):
    result=dict(ref,python_energy=p.energy,python_strength=p.strength,
                delta_energy=p.energy-ref['energy'],delta_strength=p.strength-ref['strength'],
                residual=p.residual,iterations=p.iterations,python_scf=mf.e_tot,seconds=time.time()-start)
    report.append(result)
    print(name,method,sector,p.target,f'dE={result["delta_energy"]:.3e} dZ={result["delta_strength"]:.3e}',flush=True)
   Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
