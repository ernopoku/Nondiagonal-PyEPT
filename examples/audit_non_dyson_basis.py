import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
 os.environ[key]='8'
import json,copy,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import LinearOperator,minres
from pyscf import gto,scf,lib
from nondiagonal_ept import EPT
from nondiagonal_ept.non_dyson import StaticNRL3,_opposite_operator
from nondiagonal_ept.solver import davidson,make_pole,HARTREE_TO_EV
lib.num_threads(8)
root=Path.cwd() / 'nd_nrl3_basis_audit_results'
root.mkdir(exist_ok=True)
for basis in ['cc-pvdz','aug-cc-pvdz','cc-pvtz','aug-cc-pvtz']:
 print('START',basis,flush=True); start=time.perf_counter()
 mol=gto.M(atom='F 0 0 0; H 0 0 .9168',unit='Angstrom',basis=basis,verbose=0,max_memory=16000)
 mf=scf.RHF(mol);mf.conv_tol=1e-12;mf.kernel()
 nr=EPT(mf,'NRL3',frozen=1,sector='ip',max_memory_mb=16000)
 poles=nr.kernel(targets=[4,2],tol=1e-10)
 h=StaticNRL3(nr.integrals,static_tol=1e-11,max_memory_mb=16000,static_space="full")
 ns=h.ns;no=nr.integrals.nocc//2
 def getp(ham,target):
  idx=int(np.flatnonzero(nr.integrals.original_mos==target)[0])
  e,x,r,it=davidson(ham,idx,tol=1e-10,max_cycle=300,max_space=60)
  return make_pole(ham,e,x,r,it,target)
 nd=[getp(h,p.target) for p in poles]
 action,diag,_=_opposite_operator(nr.hamiltonian,True,16000)
 def solve_at(e,rhs):
  op=LinearOperator((len(diag),)*2,matvec=lambda x:e*x-action(x),dtype=float)
  inv=1/np.maximum(abs(e-diag),1e-6)
  sol,info=minres(op,rhs,M=LinearOperator(op.shape,matvec=lambda x:inv*x,dtype=float),rtol=1e-13,maxiter=20000)
  res=float(np.linalg.norm(op@sol-rhs))
  assert res<1e-9,(info,res)
  return sol
 data={'basis':basis,'ns':ns,'static_max_residual':float(max(h.static_solve_residuals)), 'K_block_norms':{},'poles':[]}
 K=h.static_opposite
 for label,block in [('oo',K[:no,:no]),('ov',K[:no,no:]),('vv',K[no:,no:])]:
  data['K_block_norms'][label]=float(np.linalg.norm(block))
 for p,q in zip(poles,nd):
  s=q.vector[:ns]; rhs=nr.hamiltonian.bp.T@s
  dynamic=nr.hamiltonian.bp@solve_at(q.energy,rhs)
  d={'target':p.target,'NRL3_eV':p.binding_energy_ev,'nD_NRL3_eV':q.binding_energy_ev,'NRL3_PS':p.strength,'nD_PS':q.strength,'dyson_overlap_squared':float((p.normalized_dyson_mo@q.normalized_dyson_mo)**2),'nd_residual':q.residual,'static_minus_dynamic_expectation_eV':float(s@(K@s-dynamic)*HARTREE_TO_EV),'static_minus_dynamic_action_norm':float(np.linalg.norm(K@s-dynamic)),'virtual_simple_weight':float(s[no:]@s[no:])}
  data['poles'].append(d)
 # Controlled diagnostic: replace ONLY static occupied-virtual entries by
 # values sampled at the occupied orbital energy, not the virtual energy.
 K_occ=K.copy()
 for i in range(no):
  row=nr.hamiltonian.bp@solve_at(nr.integrals.energy[h.simple[i]],nr.hamiltonian.bp[i])
  K_occ[i,no:]=row[no:];K_occ[no:,i]=row[no:]
 probe=copy.copy(h);probe.a=nr.hamiltonian.a+K_occ;probe._diag=None
 data['diagnostic_only_occupied_sampled_ov_eV']=[getp(probe,p.target).binding_energy_ev for p in poles]
 data['max_abs_K_ov']=float(np.max(abs(K[:no,no:])))
 data['max_abs_occupied_sampled_ov']=float(np.max(abs(K_occ[:no,no:])))
 np.savez(root/(basis+'.npz'),K=K,K_occ=K_occ,eps=nr.integrals.energy[h.simple],a=nr.hamiltonian.a)
 data['seconds']=time.perf_counter()-start
 (root/(basis+'.json')).write_text(json.dumps(data,indent=2)+'\n')
 print(json.dumps(data),flush=True)
