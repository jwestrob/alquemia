#!/usr/bin/env python3
"""Saved Asn contact trajectories; no energy, search or coordinate construction."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import shutil
import socket
import subprocess
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/contact_review'
BASE=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015'
CANDIDATES=('origin','amide_flip','relaxed_flip_Ca','relaxed_flip_La')
METALS=('Ca','La')

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def require(v,s):
    if not v:raise ValueError(s)
def verify(p):
    require(pin(p['path'])==p,'Changed artifact: '+p['path']);return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')

def freeze(collection,output,scout_check=False):
    """Metadata-only binding after the actual remaining-stage collection exists."""
    result=read(collection);c=read(result['config']);execution=Path(collection).with_name('EXECUTION.json');ex=read(execution)
    stage='scout' if scout_check else 'remaining'
    require(result['config']==pin(BASE/f'asn_relaxation_{stage}_v1/CONFIG.json')==ex['config'],'Wrong actual relaxation stage')
    require(result['stage']==c['stage']==stage and ex['status'] in ('complete','partial','failed'),'Terminal stage receipt required')
    p=read(c['preparation']);require(len(result['rows'])==len(p['rows'])==10 and [r['source_case_id'] for r in result['rows']]==c['all_source_ids'],'All ten source rows required')
    gp=BASE/'asn_geometry_response_v1/CONFIG.json';gc=read(gp);gr=BASE/'asn_geometry_response_v1/execution/RESULT.json'
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('ANALYZE.py','RUN.sbatch','METHOD.md'):shutil.copyfile(HERE/name,out/name)
    m={'protocol':'Asn_saved_origin_flip_relaxed_named_contact_review_v1','mode':'scout_check' if scout_check else 'final_analysis','stage':stage,'implementation':pin(out/'ANALYZE.py'),
       'wrapper':pin(out/'RUN.sbatch'),'method':pin(out/'METHOD.md'),'collection':pin(collection),'execution':pin(execution),
       'run_config':result['config'],'preparation':c['preparation'],'scout_collection':c['scout_collection'],
       'prior_geometry_config':pin(gp),'prior_geometry_result':pin(gr),'geometry_helper':gc['implementation']['ANALYZE.py'],
       'python':c['python'],'source_ids':c['all_source_ids'],'candidates':list(CANDIDATES),'new_energy_calls':0,
       'output':str(out/'execution')}
    write(out/'CONFIG.json',m);print(json.dumps(pin(out/'CONFIG.json')))

def source_graph(original,h):
    cp=read(original['source_preparation']);state=read(original['state_signature']);maps={z:read(original['maps'][z]) for z in METALS}
    require(maps['Ca']==maps['La'],'Paired source maps differ')
    entries={v['qm_index']:v for v in cp['mapping']['source_to_qm']};require(set(entries)==set(range(1,204)),'Wrong204-atom mapping')
    source_qm={h.key(v['source']):i for i,v in entries.items() if 'source' in v};kin=maps['Ca']['context']
    physical_qm={i:source_qm[h.key(s)] for i,s in enumerate(kin['source_atom_metadata']) if s and h.key(s) in source_qm}
    bonds={frozenset((physical_qm[i],physical_qm[j])) for i,j in kin['bonds'] if i in physical_qm and j in physical_qm}
    n=original['named_Asn_indices_zero_based']
    for name in ('OD1','ND2','HD21','HD22','H'):
        s=entries[n[name]]['source'];require((s['chain'],s['resnum'],s['resname'],s['atom'])==('A',290,'ASN',name),'Named Asn identity differs')
    gly,=[i for i,v in entries.items() if 'source' in v and tuple(v['source'].get(k) for k in ('chain','resnum','resname','atom','element'))==('A',292,'GLY','H','H')]
    parents={tuple(a):tuple(b) for a,b in state['H_covalent_parents']}
    for name in ('HD21','HD22'):
        require(parents[h.key(entries[n[name]]['source'])]==h.key(entries[n['ND2']]['source']) and frozenset((n[name],n['ND2'])) in bonds,'Amide H parent differs')
    require(state['endpoint_charges']=={'Ca':-4,'La':-3} and state['multiplicities']=={'Ca':1,'La':1},'Wrong chemical state')
    pairs={f'{name}__Gly292_H':(n[name],gly) for name in ('HD21','HD22')}
    pairs.update({f'{name}__Asn290_H':(n[name],n['H']) for name in ('HD21','HD22')})
    require(all(frozenset(v) not in bonds for v in pairs.values()),'Named diagnostic pair is covalently bonded')
    gate=read(original['geometry_admission']);require(gate['pass'],'Original fixed-flip geometry gate failed')
    return entries,bonds,n,pairs,set(gate['cap_indices']),state

def metrics(symbols,coords,origin,flip,entries,bonds,n,pairs,caps):
    moved={n[k] for k in ('OD1','ND2','HD21','HD22')};contacts=[]
    for i in sorted(moved):
        for j in range(204):
            if j in moved or j in caps or 'H' not in (symbols[i],symbols[j]) or frozenset((i,j)) in bonds:continue
            contacts.append({'moving_index':i,'other_index':j,'moving':entries[i],
                'other':entries.get(j,{'kind':'metal','qm_index':0}),'distance_A':math.dist(coords[i],coords[j]),
                'origin_A':math.dist(origin[i],origin[j]),'fixed_flip_A':math.dist(flip[i],flip[j]),
                'both_hydrogen':symbols[i]==symbols[j]=='H'})
    contacts.sort(key=lambda v:(v['distance_A'],v['moving_index'],v['other_index']))
    hh=[v for v in contacts if v['both_hydrogen']];short=[v for v in contacts if v['distance_A']<v['fixed_flip_A']]
    named={k:{'indices':list(ij),'identities':[entries[i] for i in ij],'distance_A':math.dist(coords[ij[0]],coords[ij[1]]),
        'origin_A':math.dist(origin[ij[0]],origin[ij[1]]),'fixed_flip_A':math.dist(flip[ij[0]],flip[ij[1]]),
        'change_from_fixed_flip_A':math.dist(coords[ij[0]],coords[ij[1]])-math.dist(flip[ij[0]],flip[ij[1]])} for k,ij in pairs.items()}
    return {'metal_donor_distances_A':{k:math.dist(coords[0],coords[n[k]]) for k in ('OD1','ND2')},
        'named_nonbonded_pairs':named,'amide_H_parent_distances_A':{k:math.dist(coords[n[k]],coords[n['ND2']]) for k in ('HD21','HD22')},
        'shortest_included_nonbonded_H_pair':contacts[0] if contacts else None,
        'shortest_included_nonbonded_HH_pair':hh[0] if hh else None,
        'shortest_pair_shortened_from_fixed_flip':short[0] if short else None}

def analyze(row,prepared,prior,h):
    original=prepared['original_row'];require(row['source_case_id']==prepared['source_case_id']==prior['source_case_id'],'Source join differs')
    for k in ('maps','state_signature','source_preparation'):require(original[k]==prior['inputs'][k],'Earlier geometry state differs')
    entries,bonds,n,pairs,caps,state=source_graph(original,h);parsed={};stages={}
    for name in CANDIDATES:
        cells={z:row['matrix'][z].get(name,{'status':'unavailable','reason':'missing candidate'}) for z in METALS}
        item={'status':'unavailable','cells':cells,'boundary_flag':None,'optimizer':None,'proposal_status':None,'reason':None}
        owner=name.removeprefix('relaxed_flip_') if name.startswith('relaxed_flip_') else None
        proposal_pin=row.get('proposals',{}).get(owner) if owner else None
        try:
            proposal=read(proposal_pin) if proposal_pin else None
            if owner:
                item.update(proposal_receipt=proposal_pin,proposal_status=proposal.get('status') if proposal else 'missing',
                    optimizer=proposal.get('optimizer') if proposal else None,boundary_flag=proposal.get('boundary_flag') if proposal else None,
                    final_geometry=proposal.get('final_geometry') if proposal else None)
            require(all(cells[z]['status']=='complete' for z in METALS),'Paired candidate unavailable')
            if owner:
                require(proposal and proposal['status']=='proposal_available' and proposal['optimizer']['success'],'No admitted converged proposal')
                require(proposal['proposal']['coordinate']==cells[owner]['xyz'],'Own proposal coordinate differs')
                require(all(cells[z]['proposal_receipt']==proposal_pin for z in METALS),'Cross-pool proposal identity differs')
                for z in METALS:require(read(cells[z]['full_source_admission'])['pass'],'Final full-source admission failed')
            else:
                for z in METALS:
                    expected=original['matrix'][z]['amide_flip']['original_xyz' if name=='origin' else 'xyz']
                    require(cells[z]['xyz']==expected,'Original/fixed flip coordinate pin differs')
            pair={z:h.xyz(cells[z]['xyz']) for z in METALS}
            require(pair['Ca'][1]==pair['La'][1] and pair['Ca'][0][1:]==pair['La'][0][1:] and pair['Ca'][0][0]=='Ca' and pair['La'][0][0]=='La','Paired geometry mismatch')
            require(all(math.isfinite(v) for xyz in pair['Ca'][1] for v in xyz),'Nonfinite saved geometry')
            item.update(status='available',paired_coordinates_exact=True);parsed[name]=pair['Ca']
        except Exception as e:item['reason']=str(e)
        stages[name]=item
    require('origin' in parsed and 'amide_flip' in parsed,'Saved origin/fixed-flip geometry unavailable')
    for name,(symbols,coords,_) in parsed.items():
        try:
            require(symbols==parsed['origin'][0],'Ordered atom identity differs')
            stages[name]['metrics']=metrics(symbols,coords,parsed['origin'][1],parsed['amide_flip'][1],entries,bonds,n,pairs,caps)
            if name in ('origin','amide_flip'):
                require(stages[name]['metrics']['metal_donor_distances_A']==prior['metal_donor_distances_A'][name],'Saved4626 donor-distance replay differs')
        except Exception as e:
            stages[name].update(status='unavailable',reason=str(e));stages[name].pop('metrics',None)
    return {'source_case_id':row['source_case_id'],'conditioning':row['conditioning'],'sample':row['sample'],
        'search_source_status':row['status'],'search_source_reason':row.get('reason'),'preparation_status':prepared['status'],
        'status':'complete' if all(v['status']=='available' for v in stages.values()) else 'partial','candidates':stages,
        'source_inputs':{k:original[k] for k in ('maps','state_signature','source_preparation','geometry_admission')},
        'direct_source_bonds_excluded':True,'synthetic_caps_excluded':True,'H_parents_verified':{'HD21':'ND2','HD22':'ND2'}}

def execute(path,sha):
    job=os.environ.get('SLURM_JOB_ID');require(job and job.isdigit(),'Actual Slurm required before numerical work')
    raw=subprocess.check_output(['scontrol','show','job','-o',job],text=True);fields=dict(x.split('=',1) for x in raw.split() if '=' in x)
    require(fields['JobState']=='RUNNING' and fields['NumCPUs']=='1' and fields['NumNodes']=='1','One allocated CPU required')
    mp={'path':str(Path(path).resolve()),'sha256':sha};m=read(mp);require(m['implementation']==pin(__file__) and m['python']==pin(sys.executable),'Wrong implementation/interpreter')
    for k in ('wrapper','method','execution','prior_geometry_config','scout_collection'):
        if m[k] is not None:verify(m[k])
    c=read(m['run_config']);collection=read(m['collection']);prep=read(m['preparation']);prior=read(m['prior_geometry_result'])
    require(collection['config']==m['run_config'] and collection['preparation']==m['preparation']==c['preparation'],'Actual collection/preparation identity differs')
    require(collection['stage']==m['stage'] and m['candidates']==list(CANDIDATES),'Wrong finite review scope')
    require([r['source_case_id'] for r in collection['rows']]==m['source_ids']==c['all_source_ids'] and len(m['source_ids'])==10,'Lost source denominator')
    spec=importlib.util.spec_from_file_location('_saved_Asn_geometry_helper',verify(m['geometry_helper']));h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    out=Path(m['output']);out.mkdir(exist_ok=False);start=time.monotonic();cpu=time.process_time();rows=[]
    ex={'status':'failed','config':mp,'job_id':job,'host':socket.gethostname(),'scontrol':raw,'new_energy_calls':0}
    try:
        for r in collection['rows']:
            try:
                p,=[v for v in prep['rows'] if v['source_case_id']==r['source_case_id']];old,=[v for v in prior['rows'] if v['source_case_id']==r['source_case_id']]
                rows.append(analyze(r,p,old,h))
            except Exception as e:rows.append({'source_case_id':r['source_case_id'],'status':'unavailable','reason':str(e),'traceback':traceback.format_exc()})
        checks={}
        if m['mode']=='scout_check':
            selected=c['selected_source_ids'];require(len(selected)==1 and selected[0].endswith('__conditioned_Ca__seed-1_sample-0'),'Wrong archived scout source')
            for r,rawrow in zip(rows,collection['rows']):
                available=r.get('candidates',{})
                if r['source_case_id'] in selected:
                    require(rawrow['status']==r['status']=='complete' and all(available[n]['status']=='available' for n in CANDIDATES),'Real Ca0 proposal replay failed')
                else:
                    require(rawrow['status']=='not_executed' and r['status']=='partial','Real unexecuted source status lost')
                    require(all(available[n]['status']=='available' for n in CANDIDATES[:2]),'Unexecuted source baseline lost')
                    require(all(available[n]['status']=='unavailable' and 'metrics' not in available[n] for n in CANDIDATES[2:]),'Missing proposal substituted')
            checks={'real_Ca0_four_geometries':True,'nine_real_unexecuted_sources_preserved':True,'named_identity_and_paired_replay':True,'prior4626_baseline_replay':True}
        result={'config':mp,'declared_sources':10,'candidates':list(CANDIDATES),'rows':rows,'new_energy_calls':0,'real_artifact_checks':checks,
            'all_candidate_geometries_available':all(r['status']=='complete' for r in rows),
            'interpretation':'Named pair distances establish only how saved contacts changed. No new clash threshold, energy attribution, equilibrium or basin validation.'}
        write(out/'RESULT.json',result)
        with (out/'TABLE.csv').open('x',newline='') as f:
            writer=csv.writer(f);writer.writerow(['source_case_id','candidate','status','boundary_flag','OD1_metal_A','ND2_metal_A',*pairs_header()])
            for r in rows:
                for name in CANDIDATES:
                    v=r.get('candidates',{}).get(name,{'status':'unavailable'});q=v.get('metrics',{});d=q.get('metal_donor_distances_A',{});pairs=q.get('named_nonbonded_pairs',{})
                    writer.writerow([r['source_case_id'],name,v['status'],v.get('boundary_flag'),d.get('OD1'),d.get('ND2'),*[pairs.get(k,{}).get('distance_A') for k in pairs_header()]])
        ex.update(status='complete' if checks or all(r['status']=='complete' for r in rows) else 'partial',
            all_candidate_geometries_available=result['all_candidate_geometries_available'],real_artifact_checks=checks,
            result=pin(out/'RESULT.json'),table=pin(out/'TABLE.csv'))
    except BaseException:
        ex['traceback']=traceback.format_exc();ex['completed_analysis_rows']=rows;raise
    finally:
        ex.update(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);write(out/'EXECUTION.json',ex)

def pairs_header():return ('HD21__Gly292_H','HD22__Gly292_H','HD21__Asn290_H','HD22__Asn290_H')

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True);q=s.add_parser('freeze');q.add_argument('--collection',required=True);q.add_argument('--output',required=True);q.add_argument('--scout-check',action='store_true')
    q=s.add_parser('execute');q.add_argument('--config',required=True);q.add_argument('--config-sha',required=True);a=p.parse_args()
    if a.action=='freeze':freeze(a.collection,a.output,a.scout_check)
    else:execute(a.config,a.config_sha)
