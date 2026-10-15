#!/usr/bin/env python3
"""Metadata freeze; Slurm-only saved-graph inventory of canonical donor amides."""
import argparse, collections, functools, hashlib, json, os, re, shutil, socket, sys, time, traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/canonical_scope'
CURRENT=ROOT/'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1'
PROTOCOL='consumed_canonical25_existing_anchor_terminal_amide_scope_v1'
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def require(ok,why):
    if not ok:raise ValueError(why)
def check(p):require(pin(p['path'])==p,'Changed artifact: '+p['path']);return Path(p['path'])
@functools.lru_cache(None)
def json_path(path):return json.loads(Path(path).read_text())
def read(p):return json_path(str(check(p) if isinstance(p,dict) else Path(p)))
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def role(value):
    if isinstance(value,dict):return {'chain':value['chain'],'resnum':int(value['resnum']),'insertion_code':value.get('icode',value.get('insertion_code','')),'resname':value['resname']}
    m=re.fullmatch(r'([^:]+):([A-Z]{3})(-?\d+)([A-Za-z]?)',value)
    require(m is not None,'Unsupported original role selector')
    return dict(chain=m[1],resname=m[2],resnum=int(m[3]),insertion_code=m[4])
def same_res(a,r):return all(a.get(k,'')==r[k] for k in ('chain','resnum','insertion_code','resname'))
def freeze(output):
    native=pin(CURRENT/'manifest.json');canonical=pin(CURRENT/'CANONICAL_COLLECTION.json');n=read(native);cc=read(canonical);inv=read(n['inventory'])
    require(cc['native_manifest']==native and cc['inventory']==n['inventory'],'Current canonical authority differs')
    records=[r for r in inv['records'] if r['role']=='canonical_calibration'];require(len(records)==25 and [r['record_id'] for r in records]==[r['record_id'] for r in cc['rows']],'Exact canonical25 order differs')
    rows=[]
    for r in records:
        bm=read(r['baseline_manifest']);bc=bm['cases'][r['baseline_manifest_case_index']];require(bc['case_id']==r['case_id'],'Baseline case join differs')
        pool=read(bc['collection']);poolm=read(pool['manifest']);search=poolm['source_manifest'];sm=read(search);sc=next(c for c in sm['cases'] if c['case_id']==r['case_id'])
        parent=sc['union']['original_core']['parent'];core=read(parent);prep=read(r['source_preparation']);prot=core['protonation_manifest'];proton=read(prot)
        require(prep['source']==proton['output'],'Actual protonated source differs')
        selector=core['fixed_core']['requested_roles'].get('anchor_asparagine');tasks=[t for t in sm['tasks'] if t['case_id']==r['case_id']]
        require(len(tasks)==2 and {t['metal'] for t in tasks}=={'Ca','La'},'Paired archived task join differs')
        row={k:r[k] for k in ('record_id','case_id','source_case_id','selection_id','protein','maps','state_signature','source_preparation')}
        row.update(core_parent=parent,protonation=prot,protonated_source=proton['output'],old_pool=bc['collection'],old_pool_manifest=pool['manifest'],search_manifest=search,
            role_name='anchor_asparagine',original_role_selector=selector,role=role(selector) if selector is not None else None,
            included_role='anchor_asparagine' in core['fixed_core']['included_roles'],original_tasks=tasks,
            original_xyz={z:r['matrix'][z]['origin']['xyz'] for z in ('Ca','La')},six_saved_cells=r['matrix'])
        for p in (parent,prot,proton['output'],r['source_preparation'],r['state_signature'],*r['maps'].values(),*row['original_xyz'].values()):check(p)
        rows.append(row)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('INVENTORY.py','RUN.sbatch','PLAN.md'):shutil.copy2(HERE/name,out/name)
    m={'protocol':PROTOCOL,'implementation':pin(out/'INVENTORY.py'),'wrapper':pin(out/'RUN.sbatch'),'plan':pin(out/'PLAN.md'),
       'canonical_collection':canonical,'native_manifest':native,'inventory':n['inventory'],'existing_role_helper':pin(ROOT/'scripts/structure_informed_starts.py'),
       'python':pin('/usr/bin/python3'),'rows':rows,'declared_sources':25,'scope':'metadata and existing source graphs only; no coordinate construction or energy',
       'output':str(out/'execution'),'geometry_admission':'UNRUN','new_model_calls':0}
    write(out/'manifest.json',m);print(json.dumps({'manifest':pin(out/'manifest.json'),'role_metadata_counts':dict(collections.Counter(r['role']['resname'] if r['role'] else 'absent' for r in rows))}))

def inspect(r):
    out={k:r[k] for k in ('record_id','case_id','source_case_id','selection_id','protein','role','role_name','core_parent','source_preparation','state_signature','maps','search_manifest','protonated_source')}
    out.update(status='unsupported',reasons=[],geometry_admission='UNRUN',new_model_calls=0)
    if r['role'] is None or r['role']['resname'] not in ('ASN','GLN'):
        out.update(status='not_applicable',reasons=['Existing anchor role absent or not ASN/GLN; no remote replacement']);return out
    try:
        require(r['included_role'],'Original anchor role was not included');a=r['role'];rn=a['resname'];cp=read(r['source_preparation']);st=read(r['state_signature']);core=read(r['core_parent']);proton=read(r['protonation']);check(r['protonated_source'])
        gm={z:read(r['maps'][z]) for z in ('Ca','La')};require(gm['Ca']==gm['La'],'Paired physical mappings differ');g=gm['Ca']['context']
        require(cp['source']==proton['output']==r['protonated_source'],'Protonation/source join differs')
        out['source_completeness']={k:proton.get(k) for k in ('detected_missing_residues','repaired_missing_atom_count','repaired_missing_terminal_atom_count')}
        out['endpoint_charges']=st['endpoint_charges'];out['multiplicities']=st['multiplicities'];out['water_inventory']=st['water_inventory'];out['atom_count']=cp['new_atom_count']
        require(st['multiplicities']=={'Ca':1,'La':1},'Paired physical singlets required')
        label=f"{a['chain']}:{rn}{a['resnum']}{a['insertion_code']}";micro=[s for s in st['core_microstates'] if s[0]==label]
        require(len(micro)==1 and micro[0][1:]==[0,None],'Original neutral terminal-amide microstate unavailable');out['amide_microstate']=micro[0]
        qm={};physical={}
        for e in cp['mapping']['source_to_qm']:
            if 'source' in e and same_res(e['source'],a):
                name=e['source']['atom'];require(name not in qm,'Duplicate source-to-QM atom');qm[name]=e['qm_index']
        for i,e in enumerate(g['source_atom_metadata']):
            if e and same_res(e,a):
                name=e['atom'];require(name not in physical,'Duplicate physical source atom');physical[name]=i
        if rn=='ASN':axis=('CB','CG');terminal=('OD1','ND2','HD21','HD22');side=('CB','HB2','HB3','CG',*terminal);chi='chi2';bondnames=[('CB','CG'),('CG','OD1'),('CG','ND2'),('ND2','HD21'),('ND2','HD22'),('CB','HB2'),('CB','HB3')]
        else:axis=('CG','CD');terminal=('OE1','NE2','HE21','HE22');side=('CB','HB2','HB3','CG','HG2','HG3','CD',*terminal);chi='chi3';bondnames=[('CB','CG'),('CG','CD'),('CD','OE1'),('CD','NE2'),('NE2','HE21'),('NE2','HE22'),('CB','HB2'),('CB','HB3'),('CG','HG2'),('CG','HG3')]
        require(set(side)<=set(qm) and set(side)<=set(physical),'Incomplete named sidechain/source-to-QM map');out.update(qm_indices=qm,physical_indices=physical,terminal_atom_names=list(terminal),axis_atom_names=list(axis))
        bonds={frozenset(b) for b in g['bonds']};require(all(frozenset((physical[x],physical[y])) in bonds for x,y in bondnames),'Terminal-amide physical bond graph incomplete')
        pairs=[(h,p) for h,p in st['H_covalent_parents'] if h[:4]==[a['chain'],a['resnum'],a['insertion_code'],rn]]
        for h in terminal[2:]:
            parents=[p for hh,p in pairs if hh[4]==h];require(len(parents)==1 and parents[0][4]==terminal[1] and parents[0][:4]==[a['chain'],a['resnum'],a['insertion_code'],rn],'Amide hydrogen-parent identity differs')
        out['amide_H_parents']=[{'hydrogen':h,'parent':p} for h,p in pairs if h[4] in terminal[2:]]
        candidates=[v for v in g['modes'] if v['kind']=='sidechain_torsion' and v.get('axis_indices')==[physical[n] for n in axis]]
        require(len(candidates)==1 and candidates[0]['id'].endswith('/'+chi),'Unique existing terminal torsion unavailable');mode=candidates[0]
        require(set(mode['moving_indices'])=={physical[n] for n in (*terminal,axis[1])},'Terminal mode moves a different physical group')
        touched=[x for x in g['core_links'] if x[1]=='cap' and (x[2] in mode['moving_indices'] or x[3] in mode['moving_indices'])]
        require(not touched,'Terminal torsion moves a retained/omitted cap anchor')
        for n in side:require([qm[n],'source',physical[n]] in g['core_links'],'Physical-to-local source link differs: '+n)
        require(all(t['mapping']==r['maps'][t['metal']] and t['source_preparation']==r['source_preparation'] and t['xyz']==r['original_xyz'][t['metal']] for t in r['original_tasks']),'Archived optimizer task source join differs')
        require(r['original_tasks'][0]['active_mode_ids']==r['original_tasks'][1]['active_mode_ids'],'Paired original active modes differ')
        out.update(status='identity_and_graph_supported',terminal_mode=mode,terminal_mode_selected_by_metal={t['metal']:mode['id'] in t['active_mode_ids'] for t in r['original_tasks']},
            original_four_modes={t['metal']:t['active_mode_ids'] for t in r['original_tasks']},cap_anchor_intersections=touched,
            checks={'neutral_microstate':True,'complete_sidechain':True,'H_parents':True,'physical_bonds':True,'paired_maps':True,'exact_source_links':True,'terminal_mode':True,'no_moving_cap_anchor':True},
            limitation='Potential discrete flip only; no coordinates or full-source contact admission have been generated')
    except Exception as e:out['reasons'].append(str(e));out['traceback']=traceback.format_exc()
    return out

def execute(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before source graph interpretation')
    mp={'path':str(Path(path).resolve()),'sha256':sha};m=read(mp);require(pin(__file__)==m['implementation'] and pin(sys.executable)==m['python'],'Frozen implementation/interpreter differs')
    out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();receipt={'status':'failed','manifest':mp,'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'new_model_calls':0};rows=[]
    try:
        for k in ('wrapper','plan','canonical_collection','native_manifest','inventory','existing_role_helper'):check(m[k])
        require(m['protocol']==PROTOCOL and m['declared_sources']==25 and len(m['rows'])==25,'Declared population differs')
        require([r['record_id'] for r in m['rows']]==[r['record_id'] for r in read(m['canonical_collection'])['rows']],'Canonical identity/order differs')
        rows=[inspect(r) for r in m['rows']]
        result={'protocol':PROTOCOL,'manifest':mp,'canonical_collection':m['canonical_collection'],'declared_sources':25,'rows':rows,
            'status_counts':dict(collections.Counter(r['status'] for r in rows)),'new_model_calls':0,'coordinate_construction':False,'geometry_admission':'UNRUN','calibration_or_classification_changes':False}
        write(out/'RESULT.json',result);receipt['status']='complete'
        with (out/'SOURCES.tsv').open('x') as f:
            f.write('source_case_id\trole\tstatus\tterminal_mode\tselected_Ca\tselected_La\treason\n')
            for r in rows:f.write('\t'.join(map(str,[r['source_case_id'],r['role'],r['status'],r.get('terminal_mode',{}).get('id'),r.get('terminal_mode_selected_by_metal',{}).get('Ca'),r.get('terminal_mode_selected_by_metal',{}).get('La'),'; '.join(r['reasons'])]))+'\n')
    except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc());raise
    finally:
        receipt.update(elapsed_seconds=time.monotonic()-start,completed_source_rows=len(rows));write(out/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);a=s.add_parser('freeze');a.add_argument('--output',required=True);a=s.add_parser('execute');a.add_argument('--manifest',required=True);a.add_argument('--sha',required=True);a=p.parse_args()
    if a.op=='freeze':freeze(a.output)
    else:execute(a.manifest,a.sha)
