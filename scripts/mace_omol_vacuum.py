"""Matched vacuum/CPCM diagnostic on eight frozen real core centers."""
from __future__ import annotations
import argparse
import datetime as dt
import json
from pathlib import Path
import re
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact, BOHR_TO_A, HA_TO_KCAL, cache_key, energy, paired, read_json, record, verify, write_new, xyz
from affordable_response import read_engrad
from mace_file_checks import cached_file_checks

PROTOCOL='native_r2scan3c_matched_vacuum_response_diagnostic_v1'
SOURCE_SHA='4b3420ef72fc742b25cc407b1992073e3a55da2d28f25f15c9055781bb7db814'
METHOD='r2SCAN-3c NoAutostart DefGrid3 TightSCF EnGrad'
CASES=('GGR_extended','GGR_connected','ALPHA_1F6S','ALPHA_6IP9')
TOL={'partition_diagnostic_kcal_scale':2.,'algebra_kcal_scale':1e-7}
NUMBERS={'H':1,'C':6,'N':7,'O':8,'Ca':20,'La':57}
UNAVAILABLE={'aqueous_score':None,'calibrated_class':None,'predictive_improvement_claimed':False,
             'response_status':'response_model_not_validated','relaxation_correction_kcal_mol':None,
             'baseline_changed':False,'conditional_whole_inference_eligible':False}


def scientific_input(charge):
    return f'! {METHOD}\n* xyzfile {charge} 1 core.xyz\n'


@cached_file_checks
def sources(path):
    from ggr_sensitivity import endpoint_output_validation, executed
    from run_orca_task_manifest import load_manifest_tasks, _completed_attempt_is_valid
    if record(path)['sha256']!=SOURCE_SHA:raise InvalidArtifact('exact frozen eight-center source required')
    m=read_json(path);sources={};costs=[]
    if {t['task_id'] for t in m['tasks']}!={f'{c}_center_{m}' for c in CASES for m in ('Ca','La')}:
        raise InvalidArtifact('source center inventory differs')
    for t in m['tasks']:
        for pin in (t['xyz'],t['core_mapping'],t['core_jacobians']):verify(pin)
        ref=t['DFT_gradient']['artifacts'];old=verify(ref['input']);text=old.read_text()
        if text.replace('CPCM(Water) ','')!=scientific_input(t['charge']):
            raise InvalidArtifact('reference method is not the exact CPCM counterpart')
        if t['charge']!=(-1 if t['metal']=='Ca' else 0) or t['spin_multiplicity']!=1:
            raise InvalidArtifact('frozen physical electronic state differs')
        for pin in ref.values():verify(pin)
        if t['xyz']['sha256']!=ref['xyz']['sha256'] or t['DFT_source']['output']!=ref['output']:
            raise InvalidArtifact('source coordinates/output differ')
        receipt=read_json(verify(t['DFT_source']['receipt']));mp=verify(receipt['manifest']);dm,nt=load_manifest_tasks(mp)
        actual=next(a for a in nt if a['task_id']==receipt['task_id'])
        if not _completed_attempt_is_valid(verify(t['DFT_source']['receipt']),verify(ref['output']),
                manifest_sha256=record(mp)['sha256'],task=actual,
                runner_identity=dm['execution_policy']['task_runner'],runtime_renderer_identity=dm['execution_policy']['runtime_renderer']):
            raise InvalidArtifact('CPCM reference receipt is not valid')
        if receipt['orca_executable']!=dm['orca']:raise InvalidArtifact('source ORCA differs')
        for p in receipt['artifacts'].values():verify(p)
        endpoint_output_validation({'input':ref['input'],'xyz':t['xyz'],'charge':t['charge'],
                                    'metal':t['metal'],'multiplicity':1},verify(ref['output']),require_gradient=True)
        raw=read_engrad(verify(ref['engrad']));value=energy(verify(ref['output']))
        g=raw['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A
        if value!=t['DFT_source']['energy_hartree'] or g.tolist()!=t['DFT_gradient']['gradient_kcal_mol_per_A']:
            raise InvalidArtifact('source DFT energy/gradient does not replay')
        jac=np.load(verify(t['core_jacobians']))
        if jac.shape!=(2,len(xyz(verify(t['xyz']))),3) or not np.isfinite(jac).all():raise InvalidArtifact('physical projection invalid')
        sources[t['task_id']]=t
        wall=(dt.datetime.fromisoformat(receipt['finished_at_utc'])-dt.datetime.fromisoformat(receipt['started_at_utc'])).total_seconds()
        costs.append({'task_id':t['task_id'],'receipt':t['DFT_source']['receipt'],'wall_seconds':wall,
                      'nprocs':receipt['parallelism']['nprocs'],'rank_seconds':wall*receipt['parallelism']['nprocs']})
    return m,sources,costs


def prepare(source,agreement,output):
    from mace_global_benchmark import snapshot
    from ggr_sensitivity import ORCA
    start=time.monotonic();cpu=time.process_time();sm,old,costs=sources(source)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    extra=['mace_omol_vacuum.py','ggr_sensitivity.py','ggr_preparation.py','carve_generic.py','coordination_policy.py',
           'pqq_microstates.py','site_mechanics.py','run_orca_task_manifest.py','render_orca_runtime_input.py','result_protocol.py']
    pins=snapshot(out,(*sm['implementation'],*extra,*(p.name for p in Path(__file__).parent.glob('affordable_*.py'))))
    tasks=[]
    for tid,t in old.items():
        d=out/tid;d.mkdir();xp=d/'core.xyz';ip=d/'endpoint.inp'
        shutil.copyfile(verify(t['xyz']),xp);ip.write_text(scientific_input(t['charge']))
        task={'task_id':tid,'case_id':t['case_id'],'metal':t['metal'],'charge':t['charge'],'multiplicity':1,
              'input':record(ip),'xyz':record(xp),'output_path':str(d/'endpoint.out'),'engrad_path':str(d/'endpoint.engrad'),
              'task_type':'analytic_gradient','core_mapping':t['core_mapping'],'core_jacobians':t['core_jacobians'],
              'CPCM_source':t['DFT_source'],'CPCM_gradient':t['DFT_gradient'],'source_xyz':t['xyz'],
              'evidence':t['evidence'],'evidence_use':'consumed_method_development'}
        task['cache_key']=cache_key({'task':task,'protocol':PROTOCOL,'method':METHOD,'implementation':pins})
        tasks.append(task)
    m={'protocol_id':PROTOCOL,'agreement':record(agreement),'source_manifest':record(source),'implementation':pins,
       'execution_policy':{'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},
       'orca':record(ORCA),'tasks':tasks,'energy_scope':'isolated_vacuum_endpoint','tolerances':TOL,
       'method':METHOD,'new_DFT_calls':8,'new_MACE_calls':0,'cost_tracking':{'compute_budget':None,'wall_time_limit':None,
       'matched_CPCM_center_receipts':costs},**UNAVAILABLE}
    write_new(out/'manifest.json',m)
    result=validate(out/'manifest.json');result.update(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu)
    write_new(out/'preparation.json',result);return result


@cached_file_checks
def validate(manifest):
    from affordable_workflow import dry_run
    m=read_json(manifest);_,old,_=sources(verify(m['source_manifest']))
    if m['protocol_id']!=PROTOCOL or m['method']!=METHOD or m['tolerances']!=TOL or len(m['tasks'])!=8:
        raise InvalidArtifact('vacuum diagnostic method/inventory differs')
    for pin in m['implementation'].values():verify(pin)
    if {t['task_id'] for t in m['tasks']}!=set(old):raise InvalidArtifact('task inventory changed')
    for t in m['tasks']:
        s=old[t['task_id']];p={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':p,'protocol':PROTOCOL,'method':METHOD,'implementation':m['implementation']}):
            raise InvalidArtifact('scientific cache changed')
        if verify(t['input']).read_text()!=scientific_input(s['charge']) or t['xyz']['sha256']!=s['xyz']['sha256']:
            raise InvalidArtifact('only CPCM removal is allowed')
        for k in ('charge','metal','case_id','core_mapping','core_jacobians','evidence'):
            if t[k]!=s[k]:raise InvalidArtifact('physical source state differs: '+k)
        if t['CPCM_source']!=s['DFT_source'] or t['CPCM_gradient']!=s['DFT_gradient'] or t['multiplicity']!=1:
            raise InvalidArtifact('reference/paired state differs')
    for c in CASES:
        t={v['metal']:v for v in m['tasks'] if v['case_id']==c}
        paired(verify(t['La']['xyz']),verify(t['Ca']['xyz']),0,-1)
    return dry_run(manifest)


def parse_endpoint(task,output,engrad):
    if verify(task['input']).read_text()!=scientific_input(task['charge']):raise InvalidArtifact('exact vacuum EnGrad input required')
    text=Path(output).read_text();e=energy(output)
    if not re.search(r'Program Version\s+6\.1\.1\b',text):raise InvalidArtifact('ORCA version differs')
    if re.search(r'^\s*(?:CPCM SOLVATION MODEL|SMD SOLVATION(?: MODEL)?|COSMO SOLVATION(?: MODEL)?)\s*$',text,re.M|re.I):
        raise InvalidArtifact('vacuum endpoint contains solvent model')
    if re.search(r'numerical (?:gradient|differentiation)',text,re.I):raise InvalidArtifact('numerical gradients unsupported')
    atoms=xyz(verify(task['xyz']));numbers=[NUMBERS[a[0]] for a in atoms]
    electrons=sum(numbers)-task['charge']-(46 if task['metal']=='La' else 0)
    ecp=re.findall(r'Type\s+(\w+)\s+ECP\s+(\S+)\s+\(replacing\s+(\d+)\s+core electrons',text)
    if ecp!=([('La','Def2-ECP','46')] if task['metal']=='La' else []) or electrons%2:raise InvalidArtifact('native ECP/parity differs')
    for pattern,wanted in [(r'Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)',task['charge']),
                           (r'Multiplicity\s+Mult\s+\.{2,}\s+(\d+)',1),
                           (r'Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)',electrons)]:
        if re.findall(pattern,text)!=[str(wanted)]:raise InvalidArtifact('actual electronic state differs')
    patterns={'analytic_scf':r'ORCA SCF GRADIENT CALCULATION','native_dispersion':r'DISPERSION GRADIENT',
              'native_gcp':r'gCP correction\s+\.{2,}\s+done','total_cartesian':r'CARTESIAN GRADIENT',
              'native_D4':r'DFTD4','native_gcp_energy':r'gCP correction\s+[-+0-9.]'}
    if task['metal']=='La':patterns['native_ecp']=r'ECP gradient\s+\(SHARK\)\s+\.{2,}\s+done'
    components={k:bool(re.search(p,text,re.I)) for k,p in patterns.items()}
    if not all(components.values()):raise InvalidArtifact('native gradient/energy component absent')
    raw=read_engrad(engrad)
    if raw['atom_count']!=len(atoms) or abs(raw['energy_Ha']-e)>1e-8 or raw['atomic_numbers'].tolist()!=numbers:
        raise InvalidArtifact('analytic gradient energy/order differs')
    if not np.allclose(raw['coordinates_bohr']*BOHR_TO_A,[a[1:] for a in atoms],atol=1e-6,rtol=0):
        raise InvalidArtifact('gradient coordinates differ from frozen input')
    return {'energy_hartree':e,'gradient_kcal_mol_per_A':(raw['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A).tolist(),
            'quantity':'gradient_not_force','energy_scope':'isolated_vacuum_endpoint','explicit_electrons':electrons,
            'ecp_core_electrons':46 if ecp else 0,'native_components':components,
            'engrad':record(engrad),'environment_gradient_included':False}


def collect(manifest):
    from ggr_sensitivity import executed
    validate(manifest);m,rows=executed(manifest)
    for t in m['tasks']:
        r=rows[t['task_id']]
        if r['status']!='complete':continue
        try:r.update(parse_endpoint(t,verify(r['output']),t['engrad_path']))
        except (ValueError,OSError) as exc:r.update(status='invalid',reason=str(exc),energy_hartree=None)
    return {'manifest':record(manifest),'protocol_id':PROTOCOL,'rows':rows,
            'status':'complete' if all(v['status']=='complete' for v in rows.values()) else 'incomplete',**UNAVAILABLE}


def report(manifest,masked,neutral,output):
    from mace_hybrid import EV_TO_KCAL, accepted_attempt
    start=time.monotonic();cpu=time.process_time();c=collect(manifest);m=read_json(manifest)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);reports={}
    for label,path in [('raw_zero',masked),('shared_neutral',neutral)]:
        r=read_json(path);mp=verify(r['manifest']);mm=read_json(mp)
        for pin in mm['implementation'].values():verify(pin)
        for t in m['tasks']:
            st=next(v for v in mm['tasks'] if v['task_id']==t['task_id']);row=r['rows'][t['task_id']]
            if st['xyz']['sha256']!=t['xyz']['sha256'] or st['charge']!=t['charge']:raise InvalidArtifact('learned comparison input differs')
            if not any(accepted_attempt(a,st,mp)==row for a in (mp.parent/'execution'/t['task_id']).glob('attempt_*')):
                raise InvalidArtifact('learned comparison lacks actual receipt')
        reports[label]=r
    cases={};checks=[];partition=None
    if c['status']=='complete':
        for name in CASES:
            ends={};ts={t['metal']:t for t in m['tasks'] if t['case_id']==name}
            for metal,t in ts.items():
                v=c['rows'][t['task_id']];jac=np.load(verify(t['core_jacobians']));g=np.array(v['gradient_kcal_mol_per_A']);cg=np.array(t['CPCM_gradient']['gradient_kcal_mol_per_A'])
                project=lambda x:np.einsum('aij,ij->a',jac,x).tolist()
                ev=v['energy_hartree'];ec=t['CPCM_source']['energy_hartree']
                models={}
                for label,r in reports.items():
                    row=r['rows'][t['task_id']];gradient=row['gradient'];gmodel=np.load(verify(gradient))*EV_TO_KCAL
                    models[label]={'energy_eV':row['energy_eV'],'gradient_projection':project(gmodel),
                        'vacuum_minus_learned_projection':project(g-gmodel),'CPCM_minus_learned_projection':project(cg-gmodel)}
                ends[metal]={'vacuum_energy_hartree':ev,'CPCM_energy_hartree':ec,'CPCM_minus_vacuum_kcal_mol':(ec-ev)*HA_TO_KCAL,
                    'vacuum_gradient_projection':project(g),'CPCM_gradient_projection':project(cg),
                    'CPCM_minus_vacuum_gradient_projection':project(cg-g),'models':models}
            ca,la=ends['Ca'],ends['La'];rr={kind:(ca[kind+'_energy_hartree']-la[kind+'_energy_hartree'])*HA_TO_KCAL for kind in ('vacuum','CPCM')}
            rr.update({k:(ca['models'][k]['energy_eV']-la['models'][k]['energy_eV'])*EV_TO_KCAL for k in reports})
            cases[name]={'endpoints':ends,'R':rr,'projection_units':['kcal_scale_per_A','kcal_scale_per_radian'],
                'vacuum_grad_R_projection':(np.array(ca['vacuum_gradient_projection'])-la['vacuum_gradient_projection']).tolist()}
        a,b=cases['GGR_connected']['R'],cases['GGR_extended']['R'];delta={k:a[k]-b[k] for k in a}
        solvent=delta['CPCM']-delta['vacuum'];partition={'convention':'connected_minus_extended','shifts':delta,'CPCM_minus_vacuum':solvent,'models':{}}
        for model in reports:
            residual=delta['vacuum']-delta[model];mixed=delta['CPCM']-delta[model];closure=mixed-(solvent+residual)
            partition['models'][model]={'CPCM_mixed_residual':mixed,'vacuum_residual':residual,
                'vacuum_residual_within_2kcal_diagnostic':abs(residual)<=TOL['partition_diagnostic_kcal_scale'],'closure_error':closure}
            checks.append({'name':model+'_partition_component_closure','error':closure,'pass':abs(closure)<=TOL['algebra_kcal_scale']})
    result={**c,'cases':cases,'partition':partition,'checks':checks,'tolerances':TOL,
        'learned_reports':{'raw_zero':record(masked),'shared_neutral':record(neutral)},'implementation':record(__file__),
        'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    write_new(out/'result.json',result);return {'status':c['status'],'partition':partition}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    for k in ('source','agreement','output'):a.add_argument('--'+k,required=True)
    for op in ('dry-run','execute','collect'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
    a=s.add_parser('report')
    for k in ('manifest','masked','neutral','output'):a.add_argument('--'+k,required=True)
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.source,a.agreement,a.output)
    elif a.op=='dry-run':r=validate(a.manifest)
    elif a.op=='execute':
        from affordable_workflow import execute
        validate(a.manifest);r=execute(a.manifest)
    elif a.op=='collect':r=collect(a.manifest)
    else:r=report(a.manifest,a.masked,a.neutral,a.output)
    if a.op in ('dry-run','execute','collect') and a.output:write_new(a.output,r)
    print(json.dumps(r if a.op!='collect' else {'status':r['status'],'rows':len(r['rows'])},indent=2))
