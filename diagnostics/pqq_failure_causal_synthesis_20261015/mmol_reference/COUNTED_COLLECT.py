"""Focused copy of pinned PedE collect(): manifest atom count replaces192 only.
Scalar/population parsers and external actual-TRAH correction remain unchanged.
Its eight real ExaA outputs are replayed in the Slurm preflight.
"""
import json
from pathlib import Path
import re
from types import SimpleNamespace

def pin(path):
    import hashlib
    path=Path(path).resolve();h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(path),'sha256':h.hexdigest()}

def collect(old, m, mp, pf):
    numeric, trah_macro, population, NUMBER = old.numeric, old.trah_macro, old.population, old.NUMBER
    inventory = {r['task_id']:r for r in pf['electronic_inventory']}; rows = []
    for t in m['tasks']:
        row = {k:t[k] for k in ('task_id','sample','metal','medium','charge','multiplicity')}
        out = Path(t['output_path']); receipt = out.with_name(out.name+'.execution.json')
        if not out.exists() or not receipt.exists(): row['status']='not_completed'; rows.append(row); continue
        try:
            r = json.loads(receipt.read_text()); text = out.read_text(errors='replace'); row['execution']=pin(receipt); row['output']=pin(out)
            energy = re.findall(r'FINAL SINGLE POINT ENERGY\s+('+NUMBER+r')',text)
            nel = re.findall(r'Number of Electrons\s+NEL\s+\.+\s+(\d+)',text)
            row['energy_hartree'] = numeric(energy[-1]) if energy else None
            row['printed_explicit_electrons'] = int(nel[-1]) if nel else None
            row['scf_summary'] = [{'criterion':a.strip(),'value':numeric(b),'tolerance':numeric(c)} for a,b,c in
                re.findall(r'^\s*Last (.*?)\s+\.{3}\s*('+NUMBER+r')\s+Tolerance\s*:\s*('+NUMBER+r')',text,re.M)]
            macro = trah_macro(text)
            row['final_solver'] = 'TRAH' if macro else 'DIIS/SOSCF (see raw iteration table)'
            row['final_trah_macro'] = {'iteration':int(macro[-1][0]),'energy':numeric(macro[-1][1]),'orbital_residual':numeric(macro[-1][2])} if macro else None
            row['scf_policy_lines'] = [x.strip() for x in text.splitlines() if any(s in x for s in ('ConvCheckMode','TolE ','TolG ','TolErr ','ConvForced'))]
            row['raw_solver_lines'] = [x for x in text.splitlines() if re.match(r'^\s*\d+\s+',x) and ('MAcro' in x or re.match(r'^\s*\d+\s+-\d+\.\d+\s+',x))]
            row['populations'] = {kind:population(text,title,next(s['atom_count'] for s in m['sources'] if s['sample']==t['sample'])) for kind,title in [('Mulliken','MULLIKEN ATOMIC CHARGES'),('Loewdin','LOEWDIN ATOMIC CHARGES')]}
            source = next(s for s in m['sources'] if s['sample']==t['sample'])
            mapping = json.loads(Path(source['source_preparation']['path']).read_text())['mapping']['source_to_qm']
            groups = {'metal':[0], 'PQQ':[x['qm_index'] for x in mapping if x['kind']=='opaque_cofactor'],
                      'caps':[x['qm_index'] for x in mapping if x['kind']=='sigma_link_H'],
                      'protein':[x['qm_index'] for x in mapping if x['kind']=='source']}
            row['group_indices'] = groups; row['group_charges'] = {}
            for kind, values in row['populations'].items():
                if values is not None: row['group_charges'][kind] = {g:sum(values[i]['charge'] for i in indices) for g,indices in groups.items()}
            row['CPCM_dielectric_hartree'] = [numeric(x) for x in re.findall(r'CPCM Dielectric\s*:\s*('+NUMBER+r')',text)]
            row['dipole_lines'] = [x.strip() for x in text.splitlines() if 'dipole' in x.lower() or 'Debye' in x]
            row['basis_ecp_lines'] = [x.strip() for x in text.splitlines() if any(s in x for s in ('def2-mTZVPP','replacing 46','Def2-ECP','ECP PARAMETER'))]
            row['wall_memory_lines'] = [x.strip() for x in text.splitlines() if 'TOTAL RUN TIME' in x or 'Maximum memory' in x]
            retained = [p for p in out.parent.glob('endpoint.runtime*') if p.is_file() and (p.suffix in ('.gbw','.densities','.json') or 'property' in p.name)]
            row['retained_artifacts'] = [pin(p) for p in sorted(retained)]
            checks = {'native_receipt_success':r['returncode']==0 and r['normal_termination'] and r['scf_converged'],
                'manifest_identity':r['manifest']['sha256']==pin(mp)['sha256'],
                'output_identity':r['artifacts']['output']['sha256']==pin(out)['sha256'],
                'version':r['orca_version']=='6.1.1', 'energy_present':row['energy_hartree'] is not None,
                'electron_count':row['printed_explicit_electrons']==inventory[t['task_id']]['explicit_electrons'],
                'populations_complete':all(v is not None for v in row['populations'].values()),
                'gbw_retained':any(p.suffix=='.gbw' for p in retained),
                'native_basis':'def2-mTZVPP' in text,
                'La_ECP':t['metal']!='La' or bool(re.search(r'replacing\s+46\s+core',text))}
            if macro:
                tol = next((v['tolerance'] for v in row['scf_summary'] if v['criterion']=='Orbital Gradient'),None)
                checks['final_TRAH_gradient'] = tol is not None and row['final_trah_macro']['orbital_residual'] <= tol
            row['checks']=checks; row['status']='complete_diagnostic' if all(checks.values()) else 'failed_collection'
            row['density_response_qualification']='unestablished; solver summaries retained, no same-state repeat or derivative qualification'
        except Exception as e: row['status']='collection_error'; row['error']=f'{type(e).__name__}: {e}'
        rows.append(row)
    return {'schema':'pede_fixed_origin_dft_collection.v1','manifest':pin(mp),'cells':rows,
            'complete_cells':sum(r['status']=='complete_diagnostic' for r in rows),
            'declared_cells':len(m['tasks']),'scientific_scope':'fixed-origin solvent response only; no fit, labels, forces or production promotion'}

def corrected_collect(old, corrected, manifest, manifest_path, preflight):
    proxy=SimpleNamespace(collect=lambda m,p,pf:collect(old,m,p,pf))
    return corrected.collect(proxy,manifest,manifest_path,preflight)
