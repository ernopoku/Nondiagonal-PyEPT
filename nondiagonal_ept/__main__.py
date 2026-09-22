"""Command-line interface: python -m nondiagonal_ept input.json -o results.json."""
import argparse,json
from pathlib import Path
from pyscf import gto,scf
from . import EPT

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('-o','--output',type=Path,default=Path('ept_results.json'))
    args=parser.parse_args();config=json.loads(args.input.read_text())
    allowed={'atom','basis','charge','unit','method','frozen','sector','targets','spin',
             'scf_tol','tol','max_cycle','max_space','max_memory_mb','brueckner_tol'}
    unknown=set(config)-allowed
    if unknown:parser.error(f'Unknown input fields: {sorted(unknown)}')
    mol=gto.M(atom=config['atom'],basis=config.get('basis','cc-pvdz'),
              charge=config.get('charge',0),unit=config.get('unit','Angstrom'),spin=0,verbose=0)
    mf=scf.RHF(mol).run(conv_tol=config.get('scf_tol',1e-12))
    options={k:config[k] for k in ('frozen','sector','spin','max_memory_mb','brueckner_tol') if k in config}
    ep=EPT(mf,config.get('method','NRL3'),**options)
    solver={k:config[k] for k in ('tol','max_cycle','max_space') if k in config}
    poles=ep.kernel(config.get('targets'),**solver)
    result={'method':ep.method,'sector':ep.hamiltonian.sector,'spin':ep.hamiltonian.spin,
            'hf_energy_hartree':float(mf.e_tot),'active_original_mos':ep.integrals.original_mos.tolist(),
            'hamiltonian_dimension':ep.hamiltonian.shape[0],'poles':[]}
    for p in poles:
        result['poles'].append({'original_mo':p.target,'omega_hartree':p.energy,
            'binding_energy_ev':p.binding_energy_ev,'pole_strength':p.strength,
            'residual_norm':p.residual,'iterations':p.iterations,
            'dyson_active_mo':p.dyson_mo.tolist(),'dyson_ao':p.dyson_ao.tolist()})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(f'Saved {len(poles)} converged poles to {args.output}')

if __name__=='__main__':main()
