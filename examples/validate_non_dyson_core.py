"""Cross-basis IP/EA comparisons for HF, H2O, N2, and CO; includes HF QZ tests."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='8'
import json,time,platform,hashlib,inspect
from pathlib import Path
import numpy as np
import scipy,pyscf
from pyscf import gto,scf,lib
from nondiagonal_ept import EPT
from nondiagonal_ept.non_dyson import StaticNRL3
from nondiagonal_ept.solver import davidson,make_pole
lib.num_threads(8)
ROOT=Path.cwd() / 'nd-core-validation'; ROOT.mkdir(exist_ok=True); records=[]
cases={
'HF':('F 0 0 0; H 0 0 .9168',1,[4,2],[5]),
'H2O':('O 0 0 0; H 0 -.7586 .5043; H 0 .7586 .5043',1,[4,3],[5]),
'N2':('N 0 0 0; N 0 0 1.1136',2,[6,4,3],[7]),
'CO':('C 0 0 0; O 0 0 1.128',2,[6,5],[7]),
}
for molecule,(atom,frozen,ipt,eat) in cases.items():
 bases=['cc-pvdz','aug-cc-pvdz','cc-pvtz','aug-cc-pvtz']
 if molecule=='HF':bases+=['cc-pvqz','aug-cc-pvqz']
 for basis in bases:
  print('START',molecule,basis,flush=True)
  mol=gto.M(atom=atom,basis=basis,verbose=0,max_memory=16000)
  mf=scf.RHF(mol);mf.conv_tol=1e-12;mf.kernel();assert mf.converged
  for sector,targets in [('ip',ipt),('ea',eat)]:
   t=time.perf_counter();record=dict(molecule=molecule,atom=atom,basis=basis,frozen=frozen,sector=sector,targets=targets)
   record['non_dyson_source_sha256']=hashlib.sha256(Path(inspect.getfile(StaticNRL3)).read_bytes()).hexdigest()
   try:
    nr=EPT(mf,'NRL3',frozen=frozen,sector=sector,max_memory_mb=16000)
    ref=nr.kernel(targets,tol=1e-9,max_cycle=300,max_space=60)
    record['NRL3_seconds']=time.perf_counter()-t;t=time.perf_counter()
    nd=StaticNRL3(nr.integrals,sector=sector,static_tol=1e-10,max_memory_mb=16000,static_space="sector")
    record['nd_setup_seconds_excluding_integral_transform']=time.perf_counter()-t
    original=nr.integrals.original_mos[nr.integrals.spatial[nd.simple]]
    record['static_diagnostics']=nd.static_diagnostics;record['poles']=[]
    for p in ref:
     index=int(np.flatnonzero(original==p.target)[0]);e,x,r,it=davidson(nd,index,tol=1e-9,max_cycle=300,max_space=60)
     q=make_pole(nd,e,x,r,it,p.target)
     record['poles'].append(dict(target=p.target,NRL3_eV=p.binding_energy_ev,nD_NRL3_eV=q.binding_energy_ev,difference_eV=q.binding_energy_ev-p.binding_energy_ev,NRL3_PS=p.strength,nD_PS=q.strength,NRL3_residual=p.residual,nD_residual=q.residual,dyson_overlap_squared=float((p.normalized_dyson_mo@q.normalized_dyson_mo)**2)))
    record['dimensions']={'NRL3':nr.hamiltonian.shape[0],'nD-NRL3':nd.shape[0]}
    del nr,nd
   except Exception as exc:
    record['error']=repr(exc)
   records.append(record)
   (ROOT/'results.json').write_text(json.dumps(records,indent=2)+'\n')
   print(json.dumps({k:v for k,v in record.items() if k not in ['static_diagnostics','atom']}),flush=True)
(ROOT/'environment.json').write_text(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pyscf=pyscf.__version__),indent=2)+'\n')

if any('error' in record for record in records):
 raise SystemExit(1)
