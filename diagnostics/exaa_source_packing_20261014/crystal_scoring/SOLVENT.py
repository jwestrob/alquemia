#!/usr/bin/env python3
"""Reuse crystal MACE pools; unchanged qualified frozen-charge solvent endpoints."""
import argparse
import copy
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
import traceback

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
HERE = Path(__file__).resolve().parent
MATRIX = HERE.parent/'MATRIX.py'
NATIVE_MODEL = ROOT/'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/manifest.json'
CANONICAL = NATIVE_MODEL.parent/'CANONICAL_COLLECTION.json'
REFERENCE_SHA = 'ee00fa03aba64c23e6a737fb2b3b4a9f40dc139e58bdb9225e8ed10c9af2b58d'
ORDER = ('origin','adaptive_Ca','adaptive_La')
METALS = ('Ca','La')


def pin(path):
    p=Path(path).resolve(); h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576), b''): h.update(b)
    return {'path':str(p),'sha256':h.hexdigest()}


def verify(p):
    assert pin(p['path'])==p, ('changed input',p)
    return Path(p['path'])


def read(p):
    return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())


def write(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(obj,f,indent=2,allow_nan=False);f.write('\n')


def load(name,p):
    spec=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m);return m


def prepare(native_input,output):
    """Bind actual completed native pools; no scientific imports or new fit."""
    ip=Path(native_input).resolve(); data=read(ip)
    assert data['denominator']==2 and len(data['rows'])==2
    assert [r['chain'] for r in data['rows']]==['A','B']
    native_config=read(data['config']); worker=read(data['worker_manifest'])
    collection=read(ip.parent/'FINAL_COLLECTION.json')
    assert collection['terminal'] and collection['matrix_input']==pin(ip)
    assert collection['worker_manifest']==data['worker_manifest']==native_config['worker_manifest']
    nm=read(NATIVE_MODEL); build=read(NATIVE_MODEL.parent/'canonical_v1/BUILD.json')
    canonical=read(CANONICAL); assert pin(CANONICAL)['sha256']==REFERENCE_SHA
    assert canonical['status']=='complete' and all(v['status']=='available' for v in canonical['calibrations'].values())
    source=read(ROOT/'workspaces/orca_parameter_audit_20261012/context_gfn2_v1/manifest.json')
    cells=[]
    for row in data['rows']:
        # Incomplete native pools retain their row but do not launch a new partial search.
        if row['status']!='prepared':continue
        admission=read(row['admission']);assert admission['local_state_identical_to_Q9']
        assert read(row['state_signature'])==read(row['reference_state_signature'])
        for z in METALS:
            assert set(row['matrix'][z])==set(ORDER)
            for q in ORDER:
                cell=row['matrix'][z][q];assert cell is not None
                cells.append({'task_id':row['case_id']+'__'+z+'__'+q, 'case_id':row['case_id'],
                    'metal':z,'candidate':q, **{k:cell[k] for k in ('xyz','charge','multiplicity','atom_count','valence_electrons','MACE','MACE_eV')},
                    **{k:row[k] for k in ('maps','state_signature','source_preparation')}})
    assert len(cells) in (0,6,12) and len({t['task_id'] for t in cells})==len(cells)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    c={'protocol_id':'Nikasha_ExaA_1FLG_conditional_current_frozen_distribution_v1',
       'implementation':pin(__file__),'plan':pin(HERE/'PLAN.md'),'wrapper':pin(HERE/'SOLVENT.sbatch'),
       'matrix':pin(MATRIX),'input':pin(ip),'native_config':data['config'],
       'native_collection':pin(ip.parent/'FINAL_COLLECTION.json'),'worker_manifest':data['worker_manifest'],
       'template':native_config['template'],'MACE_model':worker['model'],'rows':data['rows'],'cells':cells,
       'factorial':pin(ROOT/'scripts/context_factorial_gfn2.py'),
       'component':pin(ROOT/'diagnostics/orca_parameter_audit_20261012/expanded_transfer_inventory/REEVALUATE.py'),
       'runtime':source['runtime'],'settings':source['settings'],'native_manifest':pin(NATIVE_MODEL),
       'native_build':pin(NATIVE_MODEL.parent/'canonical_v1/BUILD.json'),'native_binary':build['binary'],
       'native_adapter':nm['implementation']['NATIVE.py'], 'calibration':pin(CANONICAL),
       'scorer':canonical['scorer'],'frozen_calibrations':canonical['calibrations'],
       'new_MACE_calls':0,'new_calibration':False,'production_changed':False}
    write(out/'CONFIG.json',c)
    print(json.dumps({'config':pin(out/'CONFIG.json'),'complete_native_pools':len(cells)//6,
                      'planned_SCC_calls':4*len(cells),'new_MACE_calls':0}))


def mace_receipts(c):
    mr={}
    for cell in c['cells']:
        mm=read(cell['MACE']);request=read(mm['request']);verify(cell['xyz'])
        assert mm['status']=='complete' and mm['model']==c['MACE_model']
        assert mm['energy_eV']==cell['MACE_eV']
        assert request['xyz']==cell['xyz'] and request['charge']==cell['charge'] and request['multiplicity']==cell['multiplicity']
        state=mm['state']
        assert state['atoms']==cell['atom_count'] and state['charge']==cell['charge']
        assert state['spin_multiplicity']==1 and state['head']=='omol'
        verify(mm['forces'])
        mr[cell['task_id']]={'receipt':cell['MACE'],'status':'complete','energy_eV':mm['energy_eV']}
    return mr


def execute(path):
    assert os.environ.get('SLURM_JOB_ID'),'actual Slurm required before scientific imports'
    c=read(path);assert c['implementation']==pin(__file__)
    for k in ('plan','wrapper','input','native_config','native_collection','calibration','scorer'):verify(c[k])
    assert c['calibration']['sha256']==REFERENCE_SHA
    assert read(c['calibration'])['calibrations']==c['frozen_calibrations']
    matrix=load('crystal_existing_matrix',c['matrix']);place=matrix.placement()
    base=Path(path).resolve().parent/'execution'
    with (base.parent/'SOLVENT.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);base.mkdir(exist_ok=False)
        start=time.monotonic();receipt={'status':'started','config':pin(path),'placement':place,'new_MACE_calls':0}
        write(base/'START.json',receipt)
        try:
            sys.path.insert(0,str(ROOT/'scripts'))
            f=load('crystal_factorial',c['factorial']);r=load('crystal_component',c['component'])
            n=load('crystal_frozen_native',c['native_adapter']);b=load('crystal_existing_worker',read(c['worker_manifest'])['implementation'])
            view=copy.copy(c);view['rows']=[row for row in c['rows'] if row['status']=='prepared']
            write(base/'PREFLIGHT.json',matrix.preflight(view,f,b))
            mr=mace_receipts(c);write(base/'REUSED_MACE.json',{'rows':mr,'new_MACE_calls':0})
            sr=matrix.solvent(c,base,f,r,n,mr)
            result=matrix.summarize(c,mr,sr);scorer=load('crystal_existing_frozen_scorer',c['scorer'])
            for row in result['rows']:
                prepared=next(p for p in c['rows'] if p['case_id']==row['case_id'])
                energies={z:{q:({'MACE_eV':v['MACE_eV'],'solvent_hartree':v['G_qvac_hartree']} if v else None)
                             for q,v in row['matrix'][z].items()} for z in METALS}
                pool=scorer.choose(energies)
                expected=prepared['source']['expected_class'] if prepared.get('source') else 'Ca'
                assert expected=='Ca'
                row.update(pool=pool, chain=prepared['chain'], expected_class=expected,
                    conditional_decisions={policy:scorer.decision(pool[policy]['R_model_kcal_mol'] if pool[policy] else None,
                        c['frozen_calibrations'][policy]['bands'],expected) for policy in scorer.POLICIES},
                    classification_reason='conditional_current_Q9_state_calibration_transfer_crystal_redox_unverified',
                    conditional_source_state=True, source_redox_experimentally_qualified=False,
                    admission=prepared['admission'], native_source_result=prepared.get('native_source_result'),
                    native_pool_complete=prepared['status']=='prepared',
                    incomplete_native_policy='No SCC dispatch for incomplete native pool; source remains in denominator.')
            result.update(config=pin(path),frozen_calibration=c['calibration'],frozen_calibrations=c['frozen_calibrations'],
                          new_MACE_calls=0,prior_native_search_receipt=c['native_collection'],
                          new_optimization_scope='solvent_stage_only; actual prior native accommodation is retained in its receipt',
                          scope='two_consumed_crystal_chains_conditional_source_state_not_new_validation')
            write(base/'FINAL_COLLECTION.json',result)
            receipt.update(status='complete',result=pin(base/'FINAL_COLLECTION.json'))
        except BaseException as exc:
            receipt.update(status='failed',reason=type(exc).__name__+': '+str(exc),traceback=traceback.format_exc());raise
        finally:
            elapsed=time.monotonic()-start
            receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=7*elapsed)
            write(base/'EXECUTION.json',receipt)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    for name in ('prepare','prepare-run'):
        a=s.add_parser(name);a.add_argument('--native',required=True);a.add_argument('--output',required=True)
    a=s.add_parser('execute');a.add_argument('--config',required=True)
    a=p.parse_args()
    if a.command=='execute':execute(a.config)
    else:
        if a.command=='prepare-run':assert os.environ.get('SLURM_JOB_ID'),'actual Slurm required'
        prepare(a.native,a.output)
        if a.command=='prepare-run':execute(Path(a.output)/'CONFIG.json')
