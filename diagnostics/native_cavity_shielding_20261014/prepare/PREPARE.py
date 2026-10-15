#!/usr/bin/env python3
"""Metadata-only freeze; all cavity coordinate/charge preparation requires Slurm."""
import argparse
import copy
import hashlib
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

PROTOCOL = 'native_ALPB_capped_source_union_four_origins_preparation_v1'
INVENTORY_SHA = '757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343'
THREADS = ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','OPENMM_CPU_THREADS')


def pin(path):
    p=Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def check(p):
    if pin(p['path'])!=p:raise ValueError('Changed input: '+p['path'])
    return Path(p['path'])


def read(p):
    return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())


def write(path,value):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def freeze(root,output):
    root=Path(root).resolve();out=Path(output).resolve()
    ip=pin(root/'workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json')
    if ip['sha256']!=INVENTORY_SHA:raise ValueError('Wrong source inventory')
    ap=pin(root/'diagnostics/exaa_metal_position_20261014/solvent_input_audit/EVIDENCE.json')
    inv=read(ip);audit=read(ap)
    out.mkdir(parents=True,exist_ok=False);(out/'inputs').mkdir();(out/'implementation').mkdir()
    snapshots={}

    def snapshot(p):
        source=check(p);key=p['sha256']
        if key not in snapshots:
            dest=out/'inputs'/(key+source.suffix);shutil.copyfile(source,dest)
            snapshots[key]={'original':p,'snapshot':pin(dest)}
        return snapshots[key]['snapshot']

    requests=[]
    for sample in (0,3):
        sid=f'q9z4j7-pqq-la_model__conditioned_2Ca_v1__seed-1_sample-{sample}'
        src,=[r for r in inv['rows'] if r['case_id']==sid]
        prep=read(src['source_preparation']);sp=prep['source']
        source_inputs={k:snapshot(src[k]) for k in ('source_preparation','state_signature','prepared_admission','spectator_admission','source_geometry_audit')}
        source_inputs['maps']={z:snapshot(src['maps'][z]) for z in ('Ca','La')}
        source_inputs['protonated_source']=snapshot(sp)
        source_inputs['protonation']=snapshot(pin(Path(sp['path']).with_name('protonated_protonation_manifest.json')))
        for metal in ('Ca','La'):
            cell=src['matrix'][metal]['origin']
            saved,=[r for r in audit['rows'] if (r['source_case_id'],r['metal'],r['candidate'])==(sid,metal,'origin')]
            if saved['xyz']!=cell['xyz'] or saved['status']!='available' or not saved['component_v2_acceptance']['accepted']:
                raise ValueError('Original component admission/identity differs')
            req={'id':f'sample{sample}_{metal}','sample':sample,'role':src['role'],'source_case_id':sid,
                 'metal':metal,'record':src,'cell':cell,'saved_solvent':saved,'inputs':copy.deepcopy(source_inputs)}
            req['inputs'].update(xyz=snapshot(cell['xyz']),charges=snapshot(saved['native_replay']['vacuum']['charges']),
                native_result=snapshot(saved['native_result']),native_execution=snapshot(saved['native_replay']['vacuum']['execution']),
                SCC_point=snapshot(saved['SCC_endpoints']['vacuum']['point']),SCC_execution=snapshot(saved['SCC_endpoints']['vacuum']['execution']),
                SCC_stdout=snapshot(saved['SCC_endpoints']['vacuum']['stdout']),SCC_charges=snapshot(saved['SCC_endpoints']['vacuum']['charges']),
                SCC_xyz=snapshot(saved['SCC_endpoints']['vacuum']['executed_xyz']))
            requests.append(req)
    impl={}
    for name in ('PREPARE.py','RUN.sbatch'):
        dest=out/'implementation'/name;shutil.copyfile(Path(__file__).with_name(name),dest);impl[name]=pin(dest)
    m={'protocol':PROTOCOL,'inventory':ip,'solvent_audit':ap,'requests':requests,'input_snapshots':snapshots,
       'implementation':impl,'python':pin('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python'),
       'feasibility':pin(root/'diagnostics/native_cavity_shielding_20261014/feasibility/FEASIBILITY.md'),
       'settings':{'source_coordinate_tolerance_A':1e-12,'core_atom_count':192,'cap_count':14,
                   'translation_A':[11,-7,3],'rotation':'proper +90 degrees about z: (x,y,z)->(-y,x,z)',
                   'permutation':'reverse all union indices with charge/mask/identity records',
                   'source_charge':'exact zero for all appended source atoms; not a physical protonation/charge assignment',
                   'deduplication':'exact source identity, never coordinate proximity',
                   'radius_changes':False,'charge_renormalization':False,'new_hydrogens':False},
       'scout_id':'sample0_Ca','declared_endpoints':4,'output':str(out/'execution'),
       'energy_calls':0,'SCC_calls':0,'MACE_calls':0,'DFT_calls':0,'scientific_preparation':'unrun'}
    write(out/'manifest.json',m);print(json.dumps(pin(out/'manifest.json')))


def allocation():
    job=os.environ.get('SLURM_JOB_ID')
    if not job or not job.isdigit():raise RuntimeError('Scientific preparation requires actual Slurm')
    text=subprocess.check_output(['scontrol','show','job','-o',job],text=True)
    d=dict(s.split('=',1) for s in text.split() if '=' in s)
    if (d.get('JobState'),d.get('NumNodes'),d.get('NumCPUs'))!=('RUNNING','1','1'):
        raise RuntimeError('Expected one CPU on one Slurm node')
    if int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=1 or any(os.environ.get(k)!='1' for k in THREADS):
        raise RuntimeError('Expected one numerical thread')
    return {'job_id':job,'host':socket.gethostname(),'scontrol':text,'affinity':sorted(os.sched_getaffinity(0)),'allocated_CPUs':1}


def key(meta):
    return (str(meta['chain']),int(meta['resnum']),str(meta.get('insertion_code',meta.get('icode',''))).strip(),
            str(meta['resname']),str(meta['atom']),str(meta.get('altloc','')).strip())


def pdb_atoms(path):
    atoms=[];seen=set();models=0
    for line in Path(path).read_text().splitlines():
        if line.startswith('MODEL '):models+=1
        if not line.startswith(('ATOM  ','HETATM')):continue
        meta={'chain':line[21:22].strip(),'resnum':int(line[22:26]),'insertion_code':line[26:27].strip(),
              'resname':line[17:20].strip(),'atom':line[12:16].strip(),'altloc':line[16:17].strip(),
              'element':line[76:78].strip().capitalize()}
        if meta['element'] not in ('H','C','N','O','S','Ca') or meta['altloc']:
            raise ValueError('Unqualified element/alternate location in physical source')
        ident=key(meta)
        if ident in seen:raise ValueError('Duplicate source atom identity')
        seen.add(ident);atoms.append({'source':meta,'element':meta['element'],'xyz':[float(line[a:b]) for a,b in ((30,38),(38,46),(46,54))]})
    if models>1 or not atoms:raise ValueError('Missing or multiple source models')
    return atoms


def xyz_atoms(path):
    lines=Path(path).read_text().splitlines(keepends=True);n=int(lines[0]);atoms=[]
    if len(lines)!=n+2:raise ValueError('XYZ record count differs')
    for line in lines[2:]:
        fields=line.split()
        if len(fields)!=4:raise ValueError('Malformed XYZ')
        atoms.append({'element':fields[0],'xyz':[float(v) for v in fields[1:]],'line':line})
    return atoms


def close(a,b,tol):
    return max(abs(x-y) for x,y in zip(a,b))<=tol


def emit(directory,name,centers,original_xyz=None):
    target=directory/name;target.mkdir()
    xp=target/'cavity.xyz';qp=target/'charges';mp=target/'masks';jp=target/'mapping.json'
    if original_xyz is not None:shutil.copyfile(original_xyz,xp)
    else:
        with xp.open('x') as f:
            f.write(f'{len(centers)}\nExplicit cavity diagnostic; no electronic evaluation\n')
            for c in centers:f.write(c.get('line') or (c['element']+' '+' '.join(format(v,'.17g') for v in c['xyz'])+'\n'))
    with qp.open('x') as f:f.write(''.join(c['charge_text']+'\n' for c in centers))
    with mp.open('x') as f:
        f.write(str(len(centers))+'\n')
        f.write(''.join(f"{int(c['cap'])} {int(c['uncertain'])}\n" for c in centers))
    identities=[{k:v for k,v in c.items() if k not in ('line','charge_text','xyz')} for c in centers]
    write(jp,{'variant':name,'zero_based':True,'centers':[{'index':i,**x} for i,x in enumerate(identities)]})
    actual=xyz_atoms(xp)
    if [(x['element'],x['xyz']) for x in actual]!=[(x['element'],x['xyz']) for x in centers]:
        raise ValueError('Coordinate serialization changed')
    q=[float(c['charge_text']) for c in centers]
    if not all(math.isfinite(v) for c in centers for v in c['xyz']) or not all(math.isfinite(v) for v in q):
        raise ValueError('Nonfinite coordinates or charges')
    if any(float(c['charge_text'])!=0 for c in centers if not c['uncertain']):
        raise ValueError('Appended environmental charge is not exact zero')
    coords=[tuple(c['xyz']) for c in centers]
    if len(coords)!=len(set(coords)):raise ValueError('Exactly coincident cavity centers; no silent deletion')
    return {'status':'prepared','xyz':pin(xp),'charges':pin(qp),'masks':pin(mp),'mapping':pin(jp),
            'atom_count':len(centers),'charge_sum_e':math.fsum(q),'cap_indices':[i for i,c in enumerate(centers) if c['cap']],
            'uncertain_charge_indices':[i for i,c in enumerate(centers) if c['uncertain']],
            'original_local_indices':[c.get('local_index') for c in centers],'coordinate_serialization':'exact float roundtrip'}


def prepare_row(req,directory,tol):
    inp=req['inputs'];src=req['record'];saved=req['saved_solvent'];metal=req['metal'];cell=req['cell']
    cp=read(inp['source_preparation']);state=read(inp['state_signature']);sp=read(inp['spectator_admission'])
    admission=read(inp['prepared_admission']);proton=read(inp['protonation']);audit=read(inp['source_geometry_audit'])
    if state!=src['state'] or state['endpoint_charges']!={'Ca':-3,'La':-2} or state['multiplicities']!={'Ca':1,'La':1} or state['water_inventory']:
        raise ValueError('Changed dry core state')
    if admission['status']!='prepared' or (admission['atom_count'],admission['cap_count'])!=(192,14) or not sp['passed']:
        raise ValueError('Original source admission unavailable')
    if read(inp['maps']['Ca'])!=read(inp['maps']['La']):raise ValueError('Original paired maps differ')
    if cp['source']!=sp['full_protonated_source'] or proton['output']!=cp['source']:
        raise ValueError('Protonated source identity differs')
    if proton['detected_missing_residues'] or proton['repaired_missing_atom_count'] or proton['repaired_missing_terminal_atom_count']:
        raise ValueError('Source completeness unsupported')
    native=read(inp['native_result'])
    if native['status']!='available' or not all(x['passed'] and x['informative'] for x in saved['native_diagonal_checks']) or not all(x['informative'] for x in saved['vacuum_precision']):
        raise ValueError('Original native solvent qualification missing')
    if not saved['component_v2_acceptance']['accepted']:raise ValueError('Original component-v2 failed')
    atoms=xyz_atoms(check(inp['xyz']));physical=pdb_atoms(check(inp['protonated_source']));lookup={key(x['source']):x for x in physical}
    qtext=check(inp['charges']).read_text().split()
    if len(atoms)!=192 or len(qtext)!=192 or atoms[0]['element']!=metal:raise ValueError('Original atom/charge inventory differs')
    if check(inp['charges']).read_bytes()!=check(inp['SCC_charges']).read_bytes() or check(inp['xyz']).read_bytes()!=check(inp['SCC_xyz']).read_bytes():
        raise ValueError('Actual SCC inputs/populations differ from native replay')
    if (cell['charge'],cell['multiplicity'],cell['GFN2_valence_electrons'])!=(state['endpoint_charges'][metal],1,588):raise ValueError('Endpoint charge/spin differs')
    if abs(math.fsum(float(q) for q in qtext)-cell['charge'])>1e-5:raise ValueError('Raw saved charge sum fails existing tolerance')
    active=key(sp['active_selector']);spectator=key(sp['spectator_selector'])
    if active==spectator or active not in lookup or spectator not in lookup or lookup[active]['element']!='Ca' or lookup[spectator]['element']!='Ca':
        raise ValueError('Missing/ambiguous active and spectator Ca')
    if not close(atoms[0]['xyz'],lookup[active]['xyz'],tol):raise ValueError('Active source metal does not match origin')
    real_used={active};mapping={x['qm_index']:x for x in cp['mapping']['source_to_qm']}
    if set(mapping)!=set(range(1,192)):raise ValueError('Incomplete local source map')
    ar,=[r for r in audit['rows'] if r['new_source_id']==req['source_case_id']]
    hydrogens={x['name']:x for x in ar['PQQ']['validated_microstate_metadata']['hydrogens']}
    if set(hydrogens)!={'H_N1','H_C3','H_C8'}:raise ValueError('Changed physical PQQ hydrogen inventory')
    centers=[];caps=[];carver=[]
    for i,(a,q) in enumerate(zip(atoms,qtext)):
        c={**a,'local_index':i,'union_index':i,'charge_text':q,'uncertain':True,'cap':False}
        if i==0:c.update(kind='exchanged_source_metal',source=lookup[active]['source'])
        else:
            m=mapping[i]
            if not close(a['xyz'],m['xyz_A'],tol):raise ValueError('Local map coordinates differ')
            if m['kind']=='sigma_link_H':
                if a['element']!='H':raise ValueError('Non-H synthetic cap')
                if key(m['retained']) not in lookup or key(m['omitted']) not in lookup:raise ValueError('Cap source endpoints missing')
                c.update(kind='synthetic_cap',cap=True,cap_mapping=m);caps.append(i)
            elif 'source' in m:
                k=key(m['source'])
                if k in real_used or k not in lookup:raise ValueError('Duplicate/missing physical core identity')
                if lookup[k]['element']!=a['element'] or not close(lookup[k]['xyz'],a['xyz'],tol):raise ValueError('Mapped physical core geometry differs')
                real_used.add(k);c.update(kind='physical_source',source=lookup[k]['source'])
            elif m['kind']=='opaque_cofactor' and a['element']=='H' and m['name'] in hydrogens:
                h=hydrogens[m['name']]
                if not close(a['xyz'],h['xyz_A'],tol):raise ValueError('Saved carver PQQ hydrogen differs')
                c.update(kind='physical_carver_PQQ_H',cofactor_name=m['name'],cofactor_parent=h['parent_atom']);carver.append(m['name'])
            else:raise ValueError('Unmapped local physical atom')
        centers.append(c)
    if len(caps)!=14 or sorted(carver)!=sorted(hydrogens):raise ValueError('Cap/PQQ hydrogen accounting differs')
    if spectator in real_used:raise ValueError('Spectator already occurs in local core')
    for a in physical:
        if key(a['source']) in real_used:continue
        centers.append({**a,'union_index':len(centers),'local_index':None,'charge_text':'0','uncertain':False,'cap':False,'kind':'source_excluded_volume_only'})
    if len([c for c in centers if 'source' in c])!=len(physical):raise ValueError('Full physical source not represented exactly once')
    # Each cap remains explicit. Its close omitted source partner is recorded, not filtered away.
    cap_pairs=[{'local_cap_index':i,'retained':mapping[i]['retained'],'omitted':mapping[i]['omitted'],
                'cap_to_omitted_A':math.dist(atoms[i]['xyz'],lookup[key(mapping[i]['omitted'])]['xyz'])} for i in caps]
    variants={'local':emit(directory,'local',centers[:192],check(inp['xyz'])),
              'union':emit(directory,'union',centers),
              'cap_removed':emit(directory,'cap_removed',[c for c in centers if not c['cap']])}
    if req['id']=='sample0_Ca':
        translated=copy.deepcopy(centers);rotated=copy.deepcopy(centers)
        for c in translated:c['xyz']=[v+t for v,t in zip(c['xyz'],(11,-7,3))];c.pop('line',None)
        for c in rotated:
            x,y,z=c['xyz'];c['xyz']=[-y,x,z];c.pop('line',None)
        variants.update(translation=emit(directory,'translation',translated),rotation=emit(directory,'rotation',rotated),
                        permutation=emit(directory,'permutation',list(reversed(centers))))
    if variants['local']['charge_sum_e']!=variants['union']['charge_sum_e']:raise ValueError('Union changes total core charge')
    realq=math.fsum(float(c['charge_text']) for c in centers if not c['cap'])
    if realq!=variants['cap_removed']['charge_sum_e']:raise ValueError('Cap-removed control changes a real charge')
    return {'status':'prepared','id':req['id'],'sample':req['sample'],'role':req['role'],'source_case_id':req['source_case_id'],
            'metal':metal,'core_charge':cell['charge'],'multiplicity':1,'state_signature':src['state_signature'],
            'source_preparation':src['source_preparation'],'protonated_source':cp['source'],'archived_native_result':saved['native_result'],
            'saved_solvent':saved,'variants':variants,'core_atom_count':192,'physical_source_atom_count':len(physical),
            'cap_count':14,'carver_PQQ_H':sorted(carver),'active_selector':sp['active_selector'],'spectator_selector':sp['spectator_selector'],
            'cap_source_separations':cap_pairs,'zero_environment_charge_is_diagnostic_only':True,
            'physical_source_preserved_exactly':True,'all_real_core_charges_preserved':True,'radius_changes':False}


def execute(manifest,sha):
    place=allocation();mp={'path':str(Path(manifest).resolve()),'sha256':sha};m=read(mp)
    if m['protocol']!=PROTOCOL or pin(__file__)!=m['implementation']['PREPARE.py']:raise ValueError('Wrong frozen preparation')
    out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);t0=time.monotonic();cpu=time.process_time()
    rows=[{'id':r['id'],'sample':r['sample'],'source_case_id':r['source_case_id'],'metal':r['metal'],'status':'unavailable','reason':'not prepared'} for r in m['requests']]
    receipt={'status':'running','manifest':mp,'placement':place,'energy_calls':0};write(out/'START.json',receipt)
    try:
        if pin(sys.executable)!=m['python']:raise ValueError('Wrong pinned interpreter')
        for p in (m['inventory'],m['solvent_audit'],m['feasibility'],*m['implementation'].values()):check(p)
        for p in m['input_snapshots'].values():
            check(p['original']);check(p['snapshot'])
            if p['original']['sha256']!=p['snapshot']['sha256']:raise ValueError('Snapshot differs')
        for i,r in enumerate(m['requests']):
            d=out/r['id'];d.mkdir()
            try:rows[i]=prepare_row(r,d,m['settings']['source_coordinate_tolerance_A'])
            except Exception as exc:rows[i].update(reason=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc())
            write(d/'PREPARATION.json',rows[i])
        # Both endpoints must represent the same physical center order and geometry.
        for sample in (0,3):
            paired=[r for r in rows if r['sample']==sample]
            if all(r['status']=='prepared' for r in paired):
                for name in ('local','union','cap_removed'):
                    aa,bb=[xyz_atoms(check(r['variants'][name]['xyz'])) for r in paired]
                    if aa[0]['xyz']!=bb[0]['xyz'] or [(x['element'],x['xyz']) for x in aa[1:]]!=[(x['element'],x['xyz']) for x in bb[1:]]:
                        for r in paired:r.update(status='unavailable',reason='Paired physical coordinates differ')
        receipt['status']='complete' if all(r['status']=='prepared' for r in rows) else 'failed_preparation'
    except BaseException as exc:
        receipt.update(status='failed',error=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc());raise
    finally:
        receipt.update(wall_seconds=time.monotonic()-t0,process_CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(out/'EXECUTION.json',receipt)
        write(out/'PREPARATION.json',{'protocol':PROTOCOL,'manifest':mp,'terminal':True,
            'status':'prepared' if all(r['status']=='prepared' for r in rows) else 'unavailable','rows':rows,
            'declared_endpoints':4,'prepared_endpoints':sum(r['status']=='prepared' for r in rows),'scout_id':m['scout_id'],
            'execution':pin(out/'EXECUTION.json'),'energy_calls':0,'SCC_calls':0,'MACE_calls':0,'DFT_calls':0})
    if receipt['status']!='complete':raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='command',required=True)
    f=s.add_parser('freeze');f.add_argument('--root',required=True);f.add_argument('--output',required=True)
    e=s.add_parser('execute');e.add_argument('--manifest',required=True);e.add_argument('--manifest-sha',required=True)
    a=p.parse_args()
    if a.command=='freeze':freeze(a.root,a.output)
    else:execute(a.manifest,a.manifest_sha)
