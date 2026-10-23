#!/usr/bin/env python3
"""Prepare only: literal near/far common-basis trimers and unchanged model inputs."""
import argparse
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import shutil
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
W = ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
D = ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'
GIDS = tuple(f'aquo{n}_central_{m}' for n in (8,9) for m in ('La','Dy'))


def require(ok, message):
    if not ok: raise ValueError(message)


def pin(path):
    p=Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(pin(p['path'])==p, 'Changed pin: '+p['path'])
    return Path(p['path'])


def read(p): return json.loads(verify(p).read_text())


def module(p, name):
    spec=importlib.util.spec_from_file_location(name, verify(p))
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def write(p, value):
    require(not p.exists(), 'Immutable output exists')
    p.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    return pin(p)


def default_config():
    comparison=json.loads((W/'redistribution_comparison_v1/CONFIG.json').read_text())
    admission=pin(W/'reference_chemistry_v1/common_basis_decomposition_v1/collection_v1/ROOT_ADMISSION.json')
    a=read(admission); common=read(a['result'])
    return {'implementation':pin(__file__),'declaration':pin(HERE/'TRIMER_DECLARATION.md'),
        'common_adapter':pin(D/'reference_chemistry_v1/common_basis_decomposition_v1/REFERENCE.py'),
        'common_panel':common['panel'],'common_collection':a['result'],
        'common_admission':admission,'common_review':a['independent_review'],
        'model_result':comparison['model_result'],'model_review':comparison['model_independent_review'],
        'model_adapter':read(comparison['model_config'])['implementation'],
        'quantum_moments':comparison['quantum_result'],'quantum_moments_review':comparison['quantum_review'],
        'comparison_result':pin(W/'redistribution_comparison_v1/actual_v1/RESULT.json'),
        'comparison_review':pin(W/'redistribution_comparison_v1/independent_arithmetic_review_v1/REVIEW.json'),
        'maximum_new_reference_programs':16,'maximum_new_model_programs':16,
        'reference_MPI_ranks':4,'model_workers':4,'native_release':None}


def selection(g):
    n=g['water_count']; x=g['positions_A']; pairs=[]
    require(len(x)==1+3*n and all(math.isfinite(v) for row in x for v in row), 'Invalid parent geometry')
    for i,j in itertools.combinations(range(1,n+1),2):
        a=x[3*i-2]; z=x[3*j-2]
        distance=math.dist(a,z)
        u=[v-c for v,c in zip(a,x[0])]; v=[v-c for v,c in zip(z,x[0])]
        cosine=math.fsum(a*b for a,b in zip(u,v))/(math.sqrt(math.fsum(a*a for a in u))*math.sqrt(math.fsum(a*a for a in v)))
        pairs.append({'water_indices':[i,j], 'O_O_distance_A':distance,
                      'M_O_O_angle_degrees':math.degrees(math.acos(max(-1.,min(1.,cosine))))})
    return {label:min(pairs,key=lambda p:(sign*round(p['O_O_distance_A'],9),*p['water_indices']))
            for label,sign in [('near',1),('far',-1)]}


def definition(common, g, selected, label, role):
    require(role in ('metal_two_water','two_water'), 'Invalid trimer role')
    i,j=selected['water_indices']; require(1<=i<j<=g['water_count'],'Invalid water pair')
    water=[a for k in (i,j) for a in (3*k-2,3*k-1,3*k)]
    real=([0]+water) if role=='metal_two_water' else water
    c=common.definition(g,i,'metal_water' if role=='metal_two_water' else 'water_only')
    c.pop('water_index'); c.update(cell_id=f'pbe0__aquo{g["water_count"]}__{g["metal"]}__{label}__{role}',
        role=role,water_indices=[i,j],selection=selected,real_parent_atom_indices=real,
        water_parent_atom_indices=water,ghost_mask=[a not in real for a in range(len(g['symbols']))])
    z=[0 if c['ghost_mask'][a] else (11 if a==0 else common.b.Z[s]) for a,s in enumerate(c['symbols'])]
    c.update(expected_effective_nuclear_charges=z,expected_electrons=sum(z)-c['charge'],
             expected_correlated_electrons=sum(z)-c['charge'])
    require(c['expected_electrons']==(28 if role=='metal_two_water' else 20),'Electron count differs')
    return c


def proof(masks, terms, n):
    atoms=[sum(c for k,c in terms.items() if a in masks[k]) for a in range(n)]
    pairs=[sum(c for k,c in terms.items() if a in masks[k] and b in masks[k])
           for a,b in itertools.combinations(range(n),2)]
    require(atoms==[0]*n and all(v==0 for v in pairs),'One-/two-body cancellation failed')
    return {'atom_coefficients':atoms,'all_distinct_atom_pair_coefficients_zero':True,
            'atom_pairs_checked':len(pairs),'absolute_site_coefficient_sum':sum(abs(c)*len(masks[k]) for k,c in terms.items())}


def prepare(cfg, config_path, out):
    tick=time.monotonic(); cpu=time.process_time()
    require(cfg['implementation']==pin(__file__) and cfg['declaration']==pin(HERE/'TRIMER_DECLARATION.md'),'Source/declaration differs')
    require(cfg['maximum_new_reference_programs']==cfg['maximum_new_model_programs']==16 and
            cfg['reference_MPI_ranks']==cfg['model_workers']==4 and cfg['native_release'] is None,'Preparation scope differs')
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}: verify(v)
    common=module(cfg['common_adapter'],'trimer_common_helper'); model=module(cfg['model_adapter'],'trimer_model_helper')
    panel=read(cfg['common_panel']); collection=read(cfg['common_collection']); ca=read(cfg['common_admission'])
    mr=read(cfg['model_result']); mm=read(mr['manifest']); mi=read(cfg['model_review'])
    qm=read(cfg['quantum_moments']); qi=read(cfg['quantum_moments_review'])
    require(ca['result']==cfg['common_collection'] and ca['independent_review']==cfg['common_review'],'Reference admission differs')
    require(collection['new_complete']==68 and collection['reused_complete']==12 and collection['error'] is None,'Common-basis result incomplete')
    require(mi['input_result']==cfg['model_result'] and mi['status'].startswith('passed') and len(mr['rows'])==78 and all(r['status']=='complete' for r in mr['rows'].values()),'Model raw admission differs')
    require(qi['actual_result']==cfg['quantum_moments'] and qi['status'].startswith('passed'),'Quantum dipole admission differs')
    require(mm['panel']==panel['pair_panel'] and panel['adapter']==cfg['common_adapter'],'Parent panel differs')
    for name,path in common.HELPERS.items(): require(panel[name]==pin(path),'Referenced helper differs')
    for assets in panel['assets'].values():
        for p in assets.values(): verify(p)
    for n in (8,9):
        require(panel['geometries'][f'aquo{n}_central_La']['positions_A']==panel['geometries'][f'aquo{n}_central_Dy']['positions_A'],'Paired coordinate arrays differ')
    out.mkdir(parents=True,exist_ok=False); (out/'quantum_inputs').mkdir(); (out/'model_inputs').mkdir()
    result={'status':'prepared_unreleased_trimer_discriminator','implementation':cfg['implementation'],
        'configuration':pin(config_path),'declaration':cfg['declaration'],'sources':cfg,
        'candidate_result':mm['candidate_result'],'selected_parameters':mm['selected_parameters'],
        'reference_assets':panel['assets'],'orca':panel['orca'],
        'model_analyze':mr['config']['analyze'],'reference_cells':[], 'model_rows':{}, 'cycles':{},
        'selections':{},'reused_reference':{},'reused_model':{},
        'maximum_new_reference_programs':16,'maximum_new_model_programs':16,
        'new_native_calls':0,'native_release':None,'no_execution_entry_point':True}
    checks=0
    for gid in GIDS:
        g=panel['geometries'][gid]; chosen=selection(g); result['selections'][gid]=chosen
        parent=mm['rows'][gid]; xyz=[r.split() for r in verify(parent['input']['xyz']).read_text().splitlines()[1:] if r.strip()]
        require([[float(v) for v in r[2:5]] for r in xyz]==g['positions_A'],'Model/reference coordinates differ')
        for label,s in chosen.items():
            i,j=s['water_indices']; new={}
            for role in ('metal_two_water','two_water'):
                c=definition(common,g,s,label,role); qid=c['cell_id']; mid=f'{gid}__{label}__{role}'; new[role]=(qid,mid)
                path=out/'quantum_inputs'/(qid+'.inp'); path.write_text(common.deck(panel,c,4,3000)); c['template_input']=pin(path)
                text=path.read_text(); lines=text.split('* xyz ',1)[1].splitlines(); q,m=map(int,lines[0].split()); rr=[r.split() for r in lines[1:-1]]
                require(q==c['charge'] and m==1 and [r[0].endswith(':') for r in rr]==c['ghost_mask'] and [[float(v) for v in r[1:]] for r in rr]==g['positions_A'],'Rendered full-basis geometry differs')
                require(('NewECP' in text)==bool(c['expected_ecp_core']) and ('DelECP '+g['metal'] in text)==c['ghost_mask'][0] and 'EnGrad' not in text and 'D4' not in text,'ECP/electronic-SP recipe differs')
                result['reference_cells'].append(c)
                indices=c['real_parent_atom_indices']; folder=out/'model_inputs'/mid; folder.mkdir()
                (folder/'system.xyz').write_text(model.extract(xyz,indices,mid))
                for key,name in [('key','system.key'),('parameters','amoeba09.prm')]:shutil.copyfile(verify(parent['input'][key]),folder/name)
                row={k:parent[k] for k in ('metal','parent_geometry_id','source_input')}
                row.update(id=mid,kind=role,atoms=len(indices),parent_atom_indices=indices,
                    input={k:pin(folder/name) for k,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]},
                    charge=c['charge'],energy_replay_kcal_mol=None,energy_replay_authority=None)
                for key in ('types','names','native_names','masses_amu','frame_signatures'):row[key]=[parent[key][a] for a in indices]
                row['COM_A']=[math.fsum(parent['masses_amu'][a]*g['positions_A'][a][k] for a in indices)/sum(row['masses_amu']) for k in range(3)]
                require(row['input']['key']['sha256']==parent['input']['key']['sha256'] and row['input']['parameters']['sha256']==parent['input']['parameters']['sha256'],'Model terms changed')
                require(all(set(f['axis_parent_indices']).issubset(indices) and set(f['bond_parent_indices']).issubset(indices) for f in row['frame_signatures']),'Fragment cuts local frame/bond')
                result['model_rows'][mid]=row; checks+=5
            qprefix=f'pbe0__aquo{g["water_count"]}__{g["metal"]}'
            qt={new['metal_two_water'][0]:1,new['two_water'][0]:-1,
                f'{qprefix}__water_{i:02d}__metal_water':-1,f'{qprefix}__water_{j:02d}__metal_water':-1,
                f'{qprefix}__ion_ghost_shell':1,f'{qprefix}__water_{i:02d}__water_only':1,f'{qprefix}__water_{j:02d}__water_only':1}
            mt={new['metal_two_water'][1]:1,new['two_water'][1]:-1,
                f'{gid}__dimer_{i:02d}':-1,f'{gid}__dimer_{j:02d}':-1,'ion_'+g['metal']:1,
                f'{gid}__water_{i:02d}':1,f'{gid}__water_{j:02d}':1}
            for cid in qt:
                if cid not in (new['metal_two_water'][0],new['two_water'][0]):
                    record=collection['new_fragments'].get(cid,collection['reused_terms'].get(cid)); require(record is not None and record['status']=='complete' and cid in qm['rows'],'Missing reused reference')
                    result['reused_reference'][cid]={'record':record,'moment_row':qm['rows'][cid]}
            for cid in mt:
                if cid not in result['model_rows']:result['reused_model'][cid]={'definition':mm['rows'][cid],'record':mr['rows'][cid]}
            masks={cid:(result['model_rows'][cid] if cid in result['model_rows'] else mm['rows'][cid])['parent_atom_indices'] for cid in mt}
            p=proof(masks,mt,len(xyz)); require(p['absolute_site_coefficient_sum']==28,'Printing bound denominator differs')
            require(sum(qt[cid]*(3 if ('metal_two_water' in cid or '__metal_water' in cid or '__ion_ghost_shell' in cid) else 0) for cid in qt)==0,'Quantum cycle charge differs')
            result['cycles'][gid+'__'+label]={'parent_geometry_id':gid,'selection':s,'reference_terms':qt,'model_terms':mt,
                'permanent_and_pair_cancellation':p,'induced_print_bound_Debye_per_component':14e-6,
                'M_print_bound_Debye_per_component':.0035,'energy_print_bound_kcal_mol':3.5e-8}
            checks+=4
    require(len(result['reference_cells'])==len(result['model_rows'])==16 and len(result['cycles'])==8,'Finite denominator differs')
    manifest=write(out/'MANIFEST.json',result)
    receipt={'status':'prepared_unreleased','implementation':cfg['implementation'],'configuration':pin(config_path),'manifest':manifest,
             'checks_passed':checks+1,'new_native_calls':0,'reference_cells':16,'model_cells':16,'cycles':8,
             'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'local_threads':1}}
    print(json.dumps({'receipt':write(out/'RESULT.json',receipt),'manifest':manifest,'new_native_calls':0}))


if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--config'); ap.add_argument('--out',type=Path); ap.add_argument('--write-default-config',type=Path)
    a=ap.parse_args()
    if a.write_default_config:write(a.write_default_config,default_config())
    else:
        require(a.config is not None and a.out is not None,'Configuration/output required')
        prepare(json.loads(Path(a.config).read_text()),a.config,a.out.resolve())
