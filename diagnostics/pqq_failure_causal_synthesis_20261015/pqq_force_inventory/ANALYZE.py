#!/usr/bin/env python3
"""Project eight saved PQQ force arrays; scientific array loading is Slurm-only."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/pqq_force_inventory'
PROTOCOL='saved_PQQ_terminal_O_and_rigid_force_projection_v1'
LABELS=('C5_O1_terminal','C6_O2_terminal','translation_e1','translation_e2','translation_e3','rotation_e1','rotation_e2','rotation_e3')

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
    if pin(p['path'])!=p:raise ValueError('Changed pin: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,x):
    with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def require(ok,reason):
    if not ok:raise ValueError(reason)
def freeze(output):
    prior_path=ROOT/'workspaces/exaa_metal_position_20261014/force_review_v1/manifest.json';prior=read(prior_path)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('ANALYZE.py','RUN.sbatch'):shutil.copyfile(HERE/name,out/name)
    config={'protocol':PROTOCOL,'implementation':pin(out/'ANALYZE.py'),'wrapper':pin(out/'RUN.sbatch'),
      'inventory':pin(HERE/'INVENTORY.json'),'plan':pin(HERE/'PLAN.md'),'python':prior['python'],
      'environment_precedent':pin(prior_path),'declared_cells':8,'frame_norm_tolerance':1e-10,
      'rigid_projection_tolerance_eV_A':1e-9,'new_molecular_calls':0,'new_optimizations':0}
    verify(config['python']);write(out/'CONFIG.json',config);print(json.dumps(pin(out/'CONFIG.json')))

def placement():
    require(bool(os.environ.get('SLURM_JOB_ID')),'Actual Slurm required before array loading')
    text=subprocess.check_output(['scontrol','show','job','-o',os.environ['SLURM_JOB_ID']],text=True)
    import re
    fields=dict(re.findall(r'(\S+?)=(\S+)',text))
    require(fields['JobState']=='RUNNING' and int(fields['NumCPUs'])==1 and int(fields['NumNodes'])==1,'One actual shared CPU required')
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):require(os.environ.get(name)=='1','One numerical thread required')
    return {'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'scheduler':text,'CPUs':1}

def projection(np,xyz,force,named,elements,config):
    ids=list(named.values());heavy=[i for i in ids if elements[i]!='H'];require(len(ids)==27 and len(heavy)==24,'Intact PQQ count differs')
    center=np.mean(xyz[heavy],axis=0)
    e1=xyz[named['O1']]-xyz[named['C5']];e1=e1/np.linalg.norm(e1)
    e2=xyz[named['C6']]-xyz[named['C5']];e2=e2-e1*np.dot(e1,e2);e2=e2/np.linalg.norm(e2)
    e3=np.cross(e1,e2);frame=np.stack([e1,e2,e3])
    residual=float(np.max(np.abs(frame@frame.T-np.eye(3))))
    require(np.isfinite(frame).all() and residual<config['frame_norm_tolerance'],'Chemical frame failed')
    directions=[];lengths={}
    for carbon,oxygen in (('C5','O1'),('C6','O2')):
        u=xyz[named[oxygen]]-xyz[named[carbon]];lengths[carbon+'_'+oxygen]=float(np.linalg.norm(u));u=u/np.linalg.norm(u)
        direction=np.zeros_like(xyz);direction[named[oxygen]]=u;directions.append(direction)
    for axis in frame:
        direction=np.zeros_like(xyz);direction[ids]=axis;directions.append(direction)
    for axis in frame:
        direction=np.zeros_like(xyz);direction[ids]=np.cross(axis,xyz[ids]-center);directions.append(direction)
    projections={};norms={};normalization=[]
    for name,direction in zip(LABELS,directions):
        norm=float(np.sqrt(np.sum(direction[heavy]**2)));require(np.isfinite(norm) and norm>0,'Degenerate direction')
        direction=direction/norm;err=abs(float(np.sum(direction[heavy]**2))-1)
        require(err<config['frame_norm_tolerance'],'Heavy-atom direction normalization failed')
        normalization.append(err);norms[name]=norm;projections[name]=float(np.sum(force*direction))
    total=np.sum(force[ids],axis=0);torque=np.sum(np.cross(xyz[ids]-center,force[ids]),axis=0)
    return {'normalized_projections_eV_A':projections,'raw_direction_heavy_norms':norms,
      'raw_total_force_global_eV_A':total.tolist(),'raw_total_force_chemical_frame_eV_A':(frame@total).tolist(),
      'raw_torque_global_eV_rad':torque.tolist(),'raw_torque_chemical_frame_eV_rad':(frame@torque).tolist(),
      'raw_total_force_norm_eV_A':float(np.linalg.norm(total)),'raw_torque_norm_eV_rad':float(np.linalg.norm(torque)),
      'quinone_bond_lengths_A':lengths,'heavy_centroid_A':center.tolist(),'chemical_frame':frame.tolist(),
      'frame_orthonormality_residual':residual,'maximum_normalization_residual':max(normalization)}

def compare(rows):
    by={(r['sample'],r['candidate'],r['metal']):r for r in rows if r['status']=='available'}
    def values(s,q,z):return by[s,q,z]['projection']['normalized_projections_eV_A']
    def combine(a,b,sign=1,scale=1):return {k:(a[k]+sign*b[k])*scale for k in LABELS}
    paired={};motion={};source={}
    for sample in (0,3):
        paired[str(sample)]={}
        for candidate in ('origin','adaptive_La'):
            if all((sample,candidate,z) in by for z in ('Ca','La')):
                a,b=values(sample,candidate,'Ca'),values(sample,candidate,'La')
                paired[str(sample)][candidate]={'Ca_minus_La':combine(a,b,-1),'common_mean_load':combine(a,b,1,.5)}
            else:paired[str(sample)][candidate]=None
        motion[str(sample)]={}
        for z in ('Ca','La'):
            motion[str(sample)][z]=combine(values(sample,'adaptive_La',z),values(sample,'origin',z),-1) if all((sample,q,z) in by for q in ('origin','adaptive_La')) else None
        a,b=paired[str(sample)]['adaptive_La'],paired[str(sample)]['origin']
        motion[str(sample)]['paired_response']={k:combine(a[k],b[k],-1) for k in ('Ca_minus_La','common_mean_load')} if a and b else None
    for candidate in ('origin','adaptive_La'):
        a,b=paired['3'][candidate],paired['0'][candidate]
        source[candidate]={k:combine(a[k],b[k],-1) for k in ('Ca_minus_La','common_mean_load')} if a and b else None
    return {'same_geometry_metal_comparisons_eV_A':paired,'proposal_minus_origin_eV_A':motion,'sample3_minus_sample0_eV_A':source}

def execute(path):
    place=placement();c=read(path);require(c['implementation']==pin(__file__),'Wrong analysis snapshot')
    for key in ('wrapper','plan','python','environment_precedent'):verify(c[key])
    require(pin(sys.executable)==c['python'],'Actual Python differs from pinned saved-force environment')
    require(c['protocol']==PROTOCOL and c['frame_norm_tolerance']==1e-10 and c['rigid_projection_tolerance_eV_A']==1e-9,'Frozen directions/gates differ')
    import numpy as np
    inv=read(c['inventory']);require(inv['declared_endpoints']==8 and inv['DFT_forces_available']==0,'Exact saved-force scope differs')
    out=Path(path).parent/'execution';out.mkdir(exist_ok=False);start=time.monotonic();cpu=time.process_time()
    rows=[{'task_id':r['task_id'],'sample':r['sample'],'candidate':r['candidate'],'metal':r['metal'],'status':'unavailable','reason':'not analyzed'} for r in inv['rows']]
    contexts={r['sample']:r for r in inv['contexts']};coords={};receipt={'status':'failed','config':pin(path),'placement':place}
    try:
        schema=read(inv['cofactor_schema_authority']);require(schema['name_map']==inv['name_map_generic_to_CCD'] and schema['schema_bonds']==inv['cofactor_graph'],'Cofactor graph changed')
        for context in contexts.values():
            cp=read(context['source_preparation']);mapping=read(context['maps']['Ca']);read(context['state_signature']);read(context['original_core_preparation'])
            for p in context['maps'].values():verify(p)
            named=context['PQQ_indices_zero_based'];pqq=set(named.values())
            require(not any(pqq&set(m['moving_indices']) for m in mapping['context']['modes']),'PQQ overlaps a donor mode')
            actual={r['name']:r['qm_index'] for r in cp['mapping']['source_to_qm'] if r.get('fragment')==context['PQQ_fragment']}
            require(actual==named and context['microstate_id']=='pqq_ox_3minus_v1','PQQ state/mapping differs')
        for i,t in enumerate(inv['rows']):
            context=contexts[t['sample']];named=context['PQQ_indices_zero_based'];m=read(t['MACE']['receipt']);request=read(t['MACE']['request'])
            verify(t['source_xyz']);verify(request['xyz'])
            require(m['forces']==t['MACE']['forces'] and m['request']==t['MACE']['request'] and m['state']==t['MACE']['state'] and m['model']==t['MACE']['model'],'Saved MACE receipt identity differs')
            require(m['force_definition']=='negative Cartesian total vacuum OMOL energy gradient, eV/angstrom','Force convention differs')
            require(request['xyz']['sha256']==t['xyz']['sha256'] and request['charge']==t['charge'] and request['multiplicity']==t['multiplicity'],'Exact coordinate/state request differs')
            lines=verify(t['xyz']).read_text().splitlines();atoms=[line.split() for line in lines[2:]];elements=[a[0] for a in atoms];xyz=np.array([[float(v) for v in a[1:]] for a in atoms])
            force=np.load(verify(t['MACE']['forces']),allow_pickle=False)
            require(int(lines[0])==192 and xyz.shape==force.shape==(192,3) and np.isfinite(xyz).all() and np.isfinite(force).all(),'Finite192x3 coordinates/forces required')
            require(elements[0]==t['metal'],'Wrong metal identity')
            require(all(elements[a['qm_index']]==a['element'] for a in context['PQQ_atom_identity']),'Cofactor elements differ')
            key=(t['sample'],t['candidate'],t['metal']);coords[key]=(elements,xyz)
            projected=projection(np,xyz,force,named,elements,c)
            rotation=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
            transformed=projection(np,xyz@rotation.T+np.array([11.,-7.,3.]),force@rotation.T,named,elements,c)
            residual={label:transformed['normalized_projections_eV_A'][label]-projected['normalized_projections_eV_A'][label] for label in LABELS}
            require(max(abs(v) for v in residual.values())<c['rigid_projection_tolerance_eV_A'],'Joint rigid-transform projection invariance failed')
            rows[i].update(status='available',reason=None,force_receipt=t['MACE']['receipt'],forces=t['MACE']['forces'],xyz=t['xyz'],projection=projected,rigid_invariance_residuals_eV_A=residual,DFT_force_status='unavailable')
        for sample in (0,3):
            ids=list(contexts[sample]['PQQ_indices_zero_based'].values())
            for candidate in ('origin','adaptive_La'):
                ca,la=[coords[sample,candidate,z] for z in ('Ca','La')]
                require(ca[0][1:]==la[0][1:] and np.array_equal(ca[1],la[1]),'Ca/La physical coordinates differ')
            for z in ('Ca','La'):require(np.array_equal(coords[sample,'origin',z][1][ids],coords[sample,'adaptive_La',z][1][ids]),'PQQ moved during donor proposal')
        receipt['status']='complete'
    except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc());raise
    finally:
        valid=receipt['status']=='complete'
        write(out/'RESULT.json',{'protocol':PROTOCOL,'config':pin(path),'inventory':c['inventory'],'terminal':True,'status':receipt['status'],
          'rows':rows,'available_projection_cells':sum(r['status']=='available' for r in rows),'paired_coordinate_and_fixed_PQQ_checks_passed':valid,
          'comparisons':compare(rows) if valid else None,'DFT_force_available_cells':0,'new_molecular_calls':0,
          'interpretation':'Proposal-Hamiltonian residual forces only; common load is not strain energy; no composite gradient, energy correction or class label'})
        receipt.update(elapsed_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu,new_molecular_calls=0)
        write(out/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    q=sub.add_parser('freeze');q.add_argument('--output',type=Path,required=True)
    q=sub.add_parser('execute');q.add_argument('--config',type=Path,required=True);q.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.action=='freeze':freeze(a.output)
    else:execute(verify({'path':str(a.config.resolve()),'sha256':a.config_sha}))
