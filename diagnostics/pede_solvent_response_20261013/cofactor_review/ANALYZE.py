#!/usr/bin/env python3
"""Pinned PQQ chemical-graph/reference audit; analysis requires Slurm."""
from pathlib import Path
import argparse, ast, hashlib, json, os, resource, socket, subprocess, time, traceback

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
NAME_MAP = {'C1':'C8','C2':'C9','C3':'C9A','C4':'C6A','C5':'C5','C6':'C4','C7':'C3A','C8':'C1A','C9':'C2','C10':'C3','C11':'C2X','C12':'C7','C13':'C7X','C14':'C9X','N1':'N1','N2':'N6','O1':'O5','O2':'O4','O3':'O2A','O4':'O2B','O5':'O7A','O6':'O7B','O7':'O9A','O8':'O9B'}
THREADS = ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def read(p): return json.loads(Path(p).read_text())
def pin(p):
    p=Path(p).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    if pin(p['path']) != p: raise ValueError('changed input '+p['path'])
    return Path(p['path'])
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f: json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def freeze(workspace):
    out=Path(workspace).resolve();out.mkdir(parents=True,exist_ok=True)
    geo_manifest=ROOT/'workspaces/pede_failure_diagnosis_20261013/geometry_v1/manifest.json'
    geo_result=ROOT/'workspaces/pede_failure_diagnosis_20261013/geometry_local_v2/execution/RESULT.json'
    recovered=ROOT/'workspaces/ca_coverage_recovery_20261012/two_ca_solvent_c5q88_v1/manifest.json'
    controls=ROOT/'workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3/manifest.json'
    by_source={c['source_case_id']:c for p in (recovered,controls) for c in read(p)['cases'] if c['source_case_id'].startswith('q88jh5-')}
    cases=[]
    for c in read(geo_manifest)['cases']:
        s=by_source.get(c['case_id']);cases.append({**{k:c[k] for k in ('case_id','family','sample','canonical','raw','pqq','existing_admission','core')},
          'source_preparation':s['source_preparation'] if s else None,'maps':s['maps'] if s else None,
          'endpoints':{z:s['matrix'][z]['origin']['xyz'] for z in ('Ca','La')} if s else None})
    holdout=ROOT/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/prepared/holdout_preparation.json'
    refs=[{'id':c['pdb_id'],'raw':c['raw_source'],'expected_PQQ_instances':2,'state_note':('primary_paper_interprets_semiquinone_not_an_oxidized_geometry_standard' if c['pdb_id']=='1H4I' else 'redox_not_independently_assigned_here_PEG_present_in_full_crystal_omitted_in_dry_core')} for c in read(holdout)['targets']]
    flg=ROOT/'workspaces/ca_coverage_recovery_20261012/refold/run_v3_h200/prepared/1FLG_AUDIT.json'
    refs.append({'id':'1FLG','raw':read(flg)['reference'],'expected_PQQ_instances':2,'state_note':'ExaA_structural_reference_not_PedE_or_a_redox_assignment'})
    prep=ROOT/'workspaces/ca_coverage_recovery_20261012/refold/run_v3_h200/prepared/PREPARED.json'
    t=next(c for c in read(prep)['targets'] if c['accession']=='q88jh5')
    templates={'original_Ca':t['source_input'],'new_two_Ca':t['input'],'original_La':pin('/groups/banfield/projects/multienv/corkscrew/supplementary_structures/Ca-binding-controls/La-alone/q88jh5-pqq-la/q88jh5-pqq-la_data.json')}
    schema=ROOT/'workspaces/motion_envelope_20260923/pilot34_searches_v1/implementation/pqq_microstates.py'
    values={n.targets[0].id:ast.literal_eval(n.value) for n in ast.parse(schema.read_text()).body if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id in ('CCD_BONDS','PROTENIX_BONDS')}
    inputs=[pin(p) for p in (geo_manifest,geo_result,recovered,controls,holdout,flg,prep,schema,out/'PQQ_CCD.cif')]+list(templates.values())
    for c in cases:
        inputs += [c['raw']]+([c['core'],c['source_preparation'],*c['maps'].values(),*c['endpoints'].values()] if c['endpoints'] else [])
    inputs += [r['raw'] for r in refs]
    for p in inputs:verify(p)
    m={'protocol':'PedE_PQQ_saved_chemical_graph_geometry_v1','driver':pin(__file__),'plan':pin(HERE/'PLAN.md'),'python':pin('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python'),
      'inputs':list({p['path']:p for p in inputs}.values()),'geometry_result':pin(geo_result),'cases':cases,'crystals':refs,'template_inputs':templates,
      'CCD':pin(out/'PQQ_CCD.cif'),'CCD_URL':'https://files.rcsb.org/ligands/download/PQQ.cif','name_map':NAME_MAP,'schema_bonds':values,
      'tolerances':{'source_local_coordinate_A':1e-6,'paired_coordinate_A':1e-12,'constructed_H_length_A':2e-6},
      'output':str(out/'execution'),'new_molecular_calls':0,'new_coordinates':0,'numerical_threads':1}
    write(out/'manifest.json',m);print(json.dumps({'manifest':pin(out/'manifest.json'),'admitted_origins':9,'raw_sources':10,'experimental_sources':3,'executed':False}))
def xyz(p):
    lines=verify(p).read_text().splitlines();rows=[x.split() for x in lines[2:] if x.strip()]
    assert len(rows)==int(lines[0]);return [r[0] for r in rows],np.array([[float(v) for v in r[1:]] for r in rows])
def edges(bonds): return {tuple(sorted(b)) for b in bonds}
def canonical_metrics(saved,name_map,orders):
    out=[]
    for b in saved['bond_lengths']:
        a,z=b['atoms'];mapped=(name_map.get(a,a),name_map.get(z,z));key=tuple(sorted(mapped))
        out.append({'source_atoms':[a,z],'CCD_atoms':list(mapped),'distance_A':b['distance_A'],**orders[key]})
    carbonyl={a:next(b['distance_A'] for b in out if set(b['CCD_atoms'])==set(pair)) for a,pair in {'metal_facing':['C5','O5'],'other':['C4','O4']}.items()}
    acids={c:sorted(b['distance_A'] for b in out if set(b['CCD_atoms']) in ({c,a},{c,z})) for c,a,z in [('C2X','O2A','O2B'),('C7X','O7A','O7B'),('C9X','O9A','O9B')]}
    return {'bonds':out,'quinone_CO_A':carbonyl,'carboxyl_CO_sorted_A':acids,'ring_plane_RMS_A':saved['ring_plane_RMS_A'],'ring_plane_max_A':saved['ring_plane_max_A']}
def crystal_metrics(spec,bonds,orders):
    st=gemmi.read_structure(str(verify(spec['raw'])));assert len(st)==1;rows=[]
    for chain in st[0]:
        for res in chain:
            if res.name!='PQQ':continue
            atom={a.name:a for a in res if a.element.name!='H'}
            if len(atom)!=len(list(a for a in res if a.element.name!='H')):raise ValueError('ambiguous PQQ altloc')
            assert set(atom)=={a for b in bonds for a in b}
            p={k:np.array([a.pos.x,a.pos.y,a.pos.z]) for k,a in atom.items()}
            ring=sorted(k for k,a in atom.items() if a.element.name in ('C','N') and k not in ('C2X','C7X','C9X'))
            r=np.array([p[k] for k in ring]);_,_,vh=np.linalg.svd(r-r.mean(axis=0),full_matrices=False);dev=(r-r.mean(axis=0))@vh[-1]
            saved={'bond_lengths':[{'atoms':[a,b],'distance_A':float(np.linalg.norm(p[a]-p[b]))} for a,b in bonds],
                   'ring_plane_RMS_A':float(np.sqrt(np.mean(dev**2))),'ring_plane_max_A':float(np.max(np.abs(dev)))}
            rows.append({'source':spec['id'],'raw':spec['raw'],'chain':chain.name,'resnum':res.seqid.num,'icode':res.seqid.icode.strip(),
              'state_note':spec['state_note'],'atom_occupancy_B':[{'atom':k,'occupancy':a.occ,'B_iso':a.b_iso,'altloc':a.altloc.strip('\x00 ')} for k,a in atom.items()],
              **canonical_metrics(saved,{},orders)})
    assert len(rows)==spec['expected_PQQ_instances'];return rows
def local_audit(c,tolerance):
    if not c['endpoints']:return {'status':'unavailable','reason':'original preparation unavailable; raw metrics retained'}
    prep=read(verify(c['source_preparation']));core=read(verify(c['core']));f=next(f for f in core['qm_fragments'] if f['role']=='pqq_cofactor')
    assert f['formal_charge']==-3 and f['microstate_id']=='pqq_ox_3minus_v1'
    p=[v for v in prep['mapping']['source_to_qm'] if v['kind']=='opaque_cofactor' and ':PQQ' in v['fragment']]
    assert len(p)==27 and {v['name'] for v in p}==set(NAME_MAP)|{'H_N1','H_C3','H_C8'}
    sa,ca=xyz(c['endpoints']['Ca']);sl,la=xyz(c['endpoints']['La']);assert sa[0]=='Ca' and sl[0]=='La' and sa[1:]==sl[1:]
    paired=float(np.max(np.abs(ca-la)));assert paired<=tolerance['paired_coordinate_A']
    by_name={v['name']:v for v in p};errors=[]
    for v in p:
        i=v['qm_index'];assert sa[i]==v['element'];errors.append(float(np.max(np.abs(ca[i]-np.array(v['xyz_A'])))))
    assert max(errors)<=tolerance['source_local_coordinate_A']
    st=gemmi.read_structure(str(verify(c['raw'])));sel=c['pqq'];res=next(r for ch in st[0] if ch.name==sel['chain'] for r in ch if r.seqid.num==sel['resnum'] and r.seqid.icode.strip()==sel.get('icode',''))
    raw={a.name:np.array([a.pos.x,a.pos.y,a.pos.z]) for a in res if a.element.name!='H'}
    source_error=max(float(np.max(np.abs(raw[n]-ca[v['qm_index']]))) for n,v in by_name.items() if v['element']!='H');assert source_error<=tolerance['source_local_coordinate_A']
    h=[]
    for name,parent,length in [('H_N1','N1',1.01),('H_C3','C10',1.08),('H_C8','C1',1.08)]:
        actual=float(np.linalg.norm(ca[by_name[name]['qm_index']]-ca[by_name[parent]['qm_index']]));assert abs(actual-length)<=tolerance['constructed_H_length_A']
        h.append({'H':name,'parent':parent,'distance_A':actual,'declared_construction_A':length})
    caps=prep['mapping']['cut_bonds_and_caps'];bad=[v for v in caps if any(v[k].get('resname')=='PQQ' or v[k].get('canonical_resname')=='PQQ' for k in ('retained','omitted'))]
    assert not bad
    return {'status':'complete','PQQ_atoms':len(p),'PQQ_heavy_atoms':24,'PQQ_H':h,'formal_state':{k:f[k] for k in ('kind','microstate_id','formal_charge','atom_count')},
       'source_to_local_PQQ_heavy_max_error_A':source_error,'mapped_coordinate_max_error_A':max(errors),'paired_coordinate_max_error_A':paired,
       'PQQ_covalent_cut_count':len(bad),'all_caps':[{'qm_index':v['qm_index'],'kind':v['kind'],'retained':v['retained'],'omitted':v['omitted']} for v in caps],
       'electronic_localization':'not constrained or inferred from formal fragment state'}
def execute(manifest):
    job=os.environ.get('SLURM_JOB_ID','');assert job.isdigit(),'real Slurm required before scientific imports'
    allocation=subprocess.check_output(['scontrol','show','job',job,'-o'],text=True);assert 'JobId='+job+' ' in allocation
    for k in THREADS:os.environ[k]='1'
    m=read(manifest)
    for p in [m['driver'],m['plan'],m['python'],*m['inputs']]:verify(p)
    out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);t=time.monotonic();cpu=time.process_time()
    receipt={'manifest':pin(manifest),'job':job,'host':socket.gethostname(),'allocation_at_start':allocation,'new_molecular_calls':0,'new_coordinates':0};write(out/'START.json',receipt)
    try:
        global np,gemmi
        import numpy as np,gemmi
        bonds=m['schema_bonds']['CCD_BONDS'];generic=m['schema_bonds']['PROTENIX_BONDS'];assert edges(bonds)==edges([(NAME_MAP[a],NAME_MAP[b]) for a,b in generic])
        block=gemmi.cif.read_file(str(verify(m['CCD']))).sole_block();order={}
        for row in block.find('_chem_comp_bond.',['atom_id_1','atom_id_2','value_order','pdbx_aromatic_flag']):
            if row[0].startswith('H') or row[1].startswith('H'):continue
            order[tuple(sorted((row[0],row[1]))) ]={'CCD_bond_order':row[2],'CCD_aromatic':row[3]}
        assert set(order)==edges(bonds);assert order[('C5','O5')]['CCD_bond_order']=='DOUB' and order[('C4','O4')]['CCD_bond_order']=='DOUB'
        saved={r['case_id']:r for r in read(verify(m['geometry_result']))['rows']};rows=[]
        for c in m['cases']:
            g=saved[c['case_id']];entry={k:c[k] for k in ('case_id','family','sample','canonical','existing_admission')}
            entry.update(raw_metrics=canonical_metrics(g['raw']['PQQ'],NAME_MAP,order),local_metrics={k:canonical_metrics(v['PQQ'],NAME_MAP,order) for k,v in g.get('local',{}).items()})
            try:entry['preparation_audit']=local_audit(c,m['tolerances'])
            except Exception as exc:entry['preparation_audit']={'status':'failed','reason':str(exc),'traceback':traceback.format_exc()}
            rows.append(entry)
        refs=[]
        for c in m['crystals']:
            try:refs.extend(crystal_metrics(c,bonds,order))
            except Exception as exc:refs.append({'source':c['id'],'status':'unavailable','reason':str(exc),'traceback':traceback.format_exc()})
        templates={k:{'input':p,'ligands':[x['ligand'] for x in read(verify(p))['sequences'] if 'ligand' in x]} for k,p in m['template_inputs'].items()}
        sm=[next(l['smiles'] for l in v['ligands'] if l['id']=='C') for v in templates.values()];assert len(set(sm))==1
        result={'status':'complete','manifest':pin(manifest),'rows':rows,'crystals':refs,'template_inputs':templates,'same_source_PQQ_SMILES':True,
           'name_map':NAME_MAP,'graph_correspondence_verified':True,'CCD':m['CCD'],'CCD_summary':{k:block.find_value('_chem_comp.'+k) for k in ('formula','pdbx_formal_charge','pdbx_model_coordinates_db_code','pdbx_ideal_coordinates_missing_flag')},
           'geometry_result_reused':m['geometry_result'],'interpretation':'descriptive representation audit; no inferred redox, equilibrium state, repair or class promotion','new_molecular_calls':0}
        write(out/'RESULT.json',result);receipt['status']='complete';receipt['result']=pin(out/'RESULT.json')
    except Exception as exc:receipt.update(status='failed',reason=str(exc),traceback=traceback.format_exc());raise
    finally:
        receipt.update(wall_seconds=time.monotonic()-t,process_cpu_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);write(out/'RECEIPT.json',receipt)
    print(json.dumps({'status':'complete','result':receipt['result']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);s.add_parser('freeze').add_argument('--workspace',required=True);s.add_parser('execute').add_argument('--manifest',required=True);a=p.parse_args()
    freeze(a.workspace) if a.op=='freeze' else execute(a.manifest)
