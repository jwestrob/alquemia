#!/usr/bin/env python3
"""Join the eight dispatched native endpoints; no electronic/model calls."""
from pathlib import Path
import argparse, datetime as dt, hashlib, importlib.util, json, math, os, re, resource, socket, subprocess, time, traceback
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'workspaces/pede_solvent_response_20261013'
EXPECTED={f'sample{s}_{z}_{v}' for s in (1,3) for z in ('Ca','La') for v in ('water','vacuum')}
NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
def read(p): return json.loads(Path(p).read_text())
def pin(p):
 p=Path(p).resolve();h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return {'path':str(p),'sha256':h.hexdigest()}
def check(p):
 if pin(p['path'])!=p:raise ValueError('Changed input '+p['path'])
 return Path(p['path'])
def write(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def load_module(p,name):
 spec=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def freeze(workspace):
 w=Path(workspace).resolve();w.mkdir(parents=True,exist_ok=True)
 dispatch=BASE/'REFERENCE_DISPATCH_v1.json';d=read(dispatch);submission=BASE/'reference_continuation_v2/SUBMISSION.json';sub=read(submission)
 recovery=BASE/'reference_collection_recovery_v2/execution'
 gf=ROOT/'workspaces/pede_failure_diagnosis_20261013/energetics_v1/execution/RESULT.json';g=read(gf)
 jobs={t:str(j['job_id']) for j in sub['new_jobs'] for t in j['tasks']};cells=[]
 for c in d['cells']:
  mp=check(c['manifest']);m=read(mp);t=next(x for x in m['tasks'] if x['task_id']==c['task_id'])
  assert t['input']['sha256']==c['scientific_input_sha256'] and t['xyz']['sha256']==c['xyz_sha256'] and t['output_path']==c['output']
  cells.append({**c,'expected_molecular_job':'3354' if c['reuse_existing_molecular_execution'] else jobs[c['task_id']],
    'collection':str(recovery/'COLLECTION.json') if c['reuse_existing_molecular_execution'] else str(mp.parent/'continuation_v2/COLLECTION.json'),
    'per_task_collection':None if c['reuse_existing_molecular_execution'] else str(mp.parent/'continuation_v2'/f"{c['task_id']}_COLLECTION.json"),
    'lane_execution':str(recovery/'EXECUTION.json') if c['reuse_existing_molecular_execution'] else str(mp.parent/'continuation_v2/EXECUTION.json')})
 assert {x['task_id'] for x in cells}==EXPECTED and len(cells)==8
 original=read(check(d['original_manifest']));inputs=[pin(dispatch),pin(submission),pin(gf),pin(recovery/'COLLECTION.json'),pin(recovery/'EXECUTION.json')]
 inputs += [c['manifest'] for c in cells]
 for s in original['sources']:
  inputs += [s[k] for k in ('prepared_admission','source_preparation','state_signature')]+list(s['endpoints'].values())
  row=next(x for x in g['rows'] if x['source_case_id']==s['source_case_id'])
  for z in ('Ca','La'):inputs+=list(row['component_input_pins'][z+'__origin'].values())
 for c in cells:
  m=read(check(c['manifest']));t=next(t for t in m['tasks'] if t['task_id']==c['task_id']);inputs += [t[k] for k in ('input','xyz','source_xyz')]
 inputs += [g['sources']['recovered_manifest'],g['runtime']]
 for p in inputs:check(p)
 cfg={'schema':'pede_response_comparison.v1','driver':pin(__file__),'python':pin('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python'),'plan':pin(ROOT/'diagnostics/pede_solvent_response_20261013/comparison/PLAN.md'),'dispatch':pin(dispatch),'submission':pin(submission),
  'original_manifest':d['original_manifest'],'protocol_id':d['protocol_id'],'gf_saved_result':pin(gf),'cells':cells,
  'original_parser':original['implementation']['REFERENCE.py']['snapshot'],
  'corrected_parser':sub['implementation']['COLLECT_V2.py'],
  'partial_test_collection':pin(recovery/'COLLECTION.json'),'partial_test_execution':pin(recovery/'EXECUTION.json'),
  'inputs':list({p['path']:p for p in inputs}.values()),'output':str(w/'execution'),'partial_test_output':str(w/'partial_test'),
  'units':{'energy':'hartree','display_conversion_HA_TO_KCAL':g['settings']['HA_TO_KCAL'],'source':'unchanged saved corrected GFN2 settings'},
  'policy':{'no_molecular_calls':True,'missing_required_cell':'corresponding contrast unavailable','collection_precedence':'final collection, else own per-task collection; reused scout always immutable recovery collection',
   'population_schemes':['Mulliken','Loewdin'],'population_precision_e':1e-6,'mapping_coordinate_tolerance_A':1e-6,
   'contrast':'D=(U_Ca(1)-U_La(1))-(U_Ca(3)-U_La(3)); U=final_SP_water-final_SP_vacuum',
   'labels':'none; no score bands/aquo offsets or classifier decisions'}}
 check(cfg['original_parser']);check(cfg['corrected_parser']);write(w/'CONFIG.json',cfg);print(json.dumps({'config':pin(w/'CONFIG.json'),'executed':False}))

def source_groups(source,parser):
 prep=read(check(source['source_preparation']));admission=read(check(source['prepared_admission']));state=read(check(source['state_signature']))
 assert admission['atom_count']==192 and admission['cap_count']==14 and state['water_inventory']==[]
 ca=parser.readxyz(check(source['endpoints']['Ca']));la=parser.readxyz(check(source['endpoints']['La']))
 assert ca[0][0]=='Ca' and la[0][0]=='La' and ca[0][1:]==la[0][1:] and ca[1:]==la[1:]
 mapping=prep['mapping']['source_to_qm'];assert sorted(x['qm_index'] for x in mapping)==list(range(1,192))
 groups={'metal':[0],'PQQ':[],'caps':[]};residues={}
 for x in mapping:
  i=x['qm_index'];kind=x['kind'];assert max(abs(float(a)-b) for a,b in zip(ca[i][1:],x['xyz_A']))<=1e-6
  if kind=='source':
   s=x['source'];assert ca[i][0]==s['element'];key=f"{s['chain']}:{s['canonical_resname']}{s['resnum']}{s['insertion_code']}";groups.setdefault('residue:'+key,[]).append(i);residues[key]={k:s[k] for k in ('chain','canonical_resname','resnum','insertion_code')}
  elif kind=='opaque_cofactor':assert ':PQQ' in x['fragment'];assert ca[i][0]==x['element'];groups['PQQ'].append(i)
  elif kind=='sigma_link_H':assert ca[i][0]=='H';groups['caps'].append(i)
  else:raise ValueError('Unsupported mapping kind '+kind)
 assert len(groups['PQQ'])==27 and len(groups['caps'])==14 and sorted(i for ids in groups.values() for i in ids)==list(range(192))
 return {'indices':groups,'residue_metadata':residues,'group_semantics':'source-residue atoms retained in local region only; synthetic caps separate; not full-protein/formal residue charges','symbols':{'Ca':[r[0] for r in ca],'La':[r[0] for r in la]},'source_preparation':source['source_preparation'],'state_signature':source['state_signature']}

def native_cell(c,cfg,parser,corrected,groups,partial=False):
 out={'task_id':c['task_id'],'dispatch':c,'status':'unavailable','energy_hartree':None}
 try:
  m=read(check(c['manifest']));t=next(t for t in m['tasks'] if t['task_id']==c['task_id']);assert m['protocol_id']==cfg['protocol_id']
  out.update({k:t[k] for k in ('sample','metal','medium','charge','multiplicity')})
  original=read(check(cfg['original_manifest']));ot=next(t0 for t0 in original['tasks'] if t0['task_id']==t['task_id'])
  for k in ('sample','metal','medium','charge','multiplicity'):assert t[k]==ot[k]
  for k in ('input','xyz','source_xyz'):check(t[k]);assert t[k]['sha256']==ot[k]['sha256']
  assert t['input']['sha256']==c['scientific_input_sha256'] and t['xyz']['sha256']==c['xyz_sha256']
  assert t['output_path']==c['output']
  if partial:cp=check(cfg['partial_test_collection'])
  else:
   cp=Path(c['collection'])
   if not cp.exists() and c['per_task_collection'] and Path(c['per_task_collection']).exists():cp=Path(c['per_task_collection'])
  if not cp.exists():out['reason']='collection_missing';return out
  col=read(cp);out['collection']=pin(cp);assert col['schema']=='pede_fixed_origin_dft_collection.v2'
  row=next(r for r in col['cells'] if r['task_id']==t['task_id']);out['collected_status']=row['status']
  if row['status']!='complete_diagnostic':out['reason']='endpoint_not_admitted';return out
  assert col['manifest']==c['manifest'] and all(row['checks'].values())
  for k in ('sample','metal','medium','charge','multiplicity'):assert row[k]==t[k]
  ep=check(row['execution']);execution=read(ep);raw=check(row['output']);assert str(raw)==c['output']
  assert execution['task_id']==c['task_id'] and execution['manifest']==c['manifest'] and execution['orca_executable']==m['orca']
  assert execution['allocation']['slurm_job_id']==c['expected_molecular_job']
  assert execution['returncode']==0 and execution['normal_termination'] and execution['scf_converged'] and execution['orca_version']=='6.1.1'
  for k,v in [('output',row['output']),('template_input',t['input']),('xyz',t['xyz'])]:assert execution['artifacts'][k]==v;check(v)
  for k in ('runtime_input','runtime_sidecar'):check(execution['artifacts'][k])
  text=raw.read_text(errors='replace');energies=re.findall(r'FINAL SINGLE POINT ENERGY\s+('+NUMBER+')',text);assert energies
  energy=parser.numeric(energies[-1]);assert math.isfinite(energy) and energy==row['energy_hartree']
  assert 'ORCA TERMINATED NORMALLY' in text and 'def2-mTZVPP' in text
  assert t['metal']!='La' or re.search(r'replacing\s+46\s+core',text)
  header='! r2SCAN-3c NoAutostart DefGrid3 TightSCF'+(' CPCM(Water)' if t['medium']=='water' else '')
  assert check(t['input']).read_text()==header+f"\n* xyzfile {t['charge']} 1 core.xyz\n"
  if row['final_solver']=='TRAH':assert corrected.gradient_criterion(text,row['final_trah_macro']['orbital_residual'])['passed']
  g=groups[t['sample']];symbols=g['symbols'][t['metal']];atomic={'H':1,'C':6,'N':7,'O':8,'S':16,'Ca':20,'La':57}
  electrons=sum(atomic[s] for s in symbols)-t['charge']-(46 if t['metal']=='La' else 0)
  assert electrons==row['printed_explicit_electrons'] and electrons%2==0 and t['multiplicity']==1
  pop={};group_charges={};closures={}
  for scheme,title in [('Mulliken','MULLIKEN ATOMIC CHARGES'),('Loewdin','LOEWDIN ATOMIC CHARGES')]:
   p=parser.population(text,title,192);assert p==row['populations'][scheme] and p is not None
   assert [x['element'] for x in p]==symbols and all(math.isfinite(x['charge']) for x in p)
   residual=math.fsum(x['charge'] for x in p)-t['charge'];bound=192*.5*cfg['policy']['population_precision_e']+1e-9
   assert abs(residual)<=bound
   pop[scheme]=p;group_charges[scheme]={k:math.fsum(p[i]['charge'] for i in ids) for k,ids in g['indices'].items()};closures[scheme]={'residual_e':residual,'printed_rounding_bound_e':bound}
  out.update(status='available',energy_hartree=energy,output=row['output'],native_execution=row['execution'],execution=execution,
    xyz=t['xyz'],state_signature=g['state_signature'],explicit_electrons=electrons,checks=row['checks'],
    final_solver=row['final_solver'],final_TRAH_criterion=row.get('final_TRAH_criterion'),scf_summary=row['scf_summary'],
    population_schemes=pop,group_charges_e=group_charges,population_sum=closures,
    CPCM_dielectric_hartree_diagnostic_only=row['CPCM_dielectric_hartree'],density_response_qualification=row['density_response_qualification'])
  lp=Path(c['lane_execution'])
  if lp.exists():out['lane_execution']={'pin':pin(lp),'record':read(lp)}
 except Exception as exc:out.update(status='unavailable',energy_hartree=None,reason=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc())
 return out

def scalar(terms,rows,factor):
 missing=[k for k in terms if rows[k]['status']!='available']
 if missing:return {'status':'unavailable','missing_cells':missing,'energy_hartree':None,'energy_kcal_mol':None}
 value=math.fsum(a*rows[k]['energy_hartree'] for k,a in terms.items())
 return {'status':'available','coefficients':terms,'energy_hartree':value,'energy_kcal_mol':value*factor}
def contrasts(rows,factor):
 u={f'sample{s}_{z}':scalar({f'sample{s}_{z}_water':1,f'sample{s}_{z}_vacuum':-1},rows,factor) for s in (1,3) for z in ('Ca','La')}
 paired={f'sample{s}':scalar({f'sample{s}_{z}_{v}':(1 if z=='Ca' else -1)*(1 if v=='water' else -1) for z in ('Ca','La') for v in ('water','vacuum')},rows,factor) for s in (1,3)}
 d=scalar({f'sample{s}_{z}_{v}':(1 if s==1 else -1)*(1 if z=='Ca' else -1)*(1 if v=='water' else -1) for s in (1,3) for z in ('Ca','La') for v in ('water','vacuum')},rows,factor)
 response={z:{v:scalar({f'sample1_{z}_{v}':1,f'sample3_{z}_{v}':-1},rows,factor) for v in ('water','vacuum')} for z in ('Ca','La')}
 for z in ('Ca','La'):response[z]['transfer']=scalar({f'sample{s}_{z}_{v}':(1 if s==1 else -1)*(1 if v=='water' else -1) for s in (1,3) for v in ('water','vacuum')},rows,factor)
 return {'U_water_minus_vacuum':u,'paired_U_Ca_minus_La':paired,'D_sample1_minus_sample3':d,'source_response_sample1_minus_sample3':response}

def charge_response(rows,groups):
 result={}
 for scheme in ('Mulliken','Loewdin'):
  pairs={}
  for s in (1,3):
   for z in ('Ca','La'):
    a,b=rows[f'sample{s}_{z}_water'],rows[f'sample{s}_{z}_vacuum'];key=f'sample{s}_{z}'
    if a['status']!='available' or b['status']!='available':pairs[key]={'status':'unavailable','reason':'water/vacuum pair incomplete'};continue
    assert set(a['group_charges_e'][scheme])==set(b['group_charges_e'][scheme])
    pairs[key]={'status':'available','water_minus_vacuum_e':{k:a['group_charges_e'][scheme][k]-b['group_charges_e'][scheme][k] for k in a['group_charges_e'][scheme]}}
  source={}
  for z in ('Ca','La'):
   source[z]={}
   for v in ('water','vacuum'):
    a,b=rows[f'sample1_{z}_{v}'],rows[f'sample3_{z}_{v}']
    if a['status']!='available' or b['status']!='available':source[z][v]={'status':'unavailable','reason':'source pair incomplete'};continue
    assert groups[1]['indices']==groups[3]['indices'] and set(a['group_charges_e'][scheme])==set(b['group_charges_e'][scheme])
    source[z][v]={'status':'available','sample1_minus_sample3_e':{k:a['group_charges_e'][scheme][k]-b['group_charges_e'][scheme][k] for k in a['group_charges_e'][scheme]}}
  result[scheme]={'water_minus_vacuum':pairs,'source_changes':source,'interpretation':'within-scheme population response only; no equivalence to other schemes or integer redox assignment'}
 return result

def gfn_rows(cfg,original):
 saved=read(check(cfg['gf_saved_result']));check(saved['runtime']);manifest=read(check(saved['sources']['recovered_manifest']));rows={}
 for source in original['sources']:
  c=next(x for x in saved['rows'] if x['source_case_id']==source['source_case_id']);assert c['state_signature']['sha256']==source['state_signature']['sha256']
  mc=next(x for x in manifest['cases'] if x['source_case_id']==source['source_case_id'])
  for z in ('Ca','La'):
   assert mc['matrix'][z]['origin']['xyz']==source['endpoints'][z]
   for v,field in [('water','GFN2_ALPB_Eh'),('vacuum','GFN2_vacuum_Eh')]:
    p=c['component_input_pins'][z+'__origin'][v];r=read(check(p));t=r['task'];b=r['backends']['xtb']
    assert t['xyz']==source['endpoints'][z] and t['metal']==z and t['medium']==v and t['charge']==(-3 if z=='Ca' else -2) and t['multiplicity']==1
    assert b['status']=='available' and b['energy_hartree']==c['cells'][z]['origin'][field]
    rows[f"sample{source['sample']}_{z}_{v}"]={'status':'available','energy_hartree':b['energy_hartree'],'point':p,'xyz':t['xyz'],'state_signature':source['state_signature']}
 return rows

def compare(cfg,partial=False):
 parser=load_module(cfg['original_parser'],'frozen_reference');corrected=load_module(cfg['corrected_parser'],'corrected_reference')
 original=read(check(cfg['original_manifest']));assert original['protocol_id']==cfg['protocol_id']
 groups={s['sample']:source_groups(s,parser) for s in original['sources']};assert groups[1]['state_signature']['sha256']==groups[3]['state_signature']['sha256']
 assert groups[1]['indices']==groups[3]['indices']
 cells={c['task_id']:native_cell(c,cfg,parser,corrected,groups,partial) for c in cfg['cells']};assert set(cells)==EXPECTED
 factor=cfg['units']['display_conversion_HA_TO_KCAL'];dft=contrasts(cells,factor);gr=gfn_rows(cfg,original);gfn=contrasts(gr,factor)
 discrepancy={'status':'unavailable','energy_hartree':None,'energy_kcal_mol':None}
 if dft['D_sample1_minus_sample3']['status']=='available':
  v=dft['D_sample1_minus_sample3']['energy_hartree']-gfn['D_sample1_minus_sample3']['energy_hartree'];discrepancy={'status':'available','energy_hartree':v,'energy_kcal_mol':v*factor}
 return {'status':'complete' if all(c['status']=='available' for c in cells.values()) else 'partial','protocol_id':cfg['protocol_id'],'declared_cells':8,'available_cells':sum(c['status']=='available' for c in cells.values()),
  'cells':cells,'DFT_CPCM':dft,'GFN2_ALPB':{'cells':gr,**gfn},'D_DFT_minus_GFN2':discrepancy,'group_maps':groups,'population_changes':charge_response(cells,groups),
  'units':cfg['units'],'new_molecular_calls':0,'interpretation':'fixed-origin response diagnostic; CPCM and ALPB differ in electronic and solvent models; no affinity, class decision, replacement solvent correction or numerical-density qualification claimed'}

def execute(config,partial=False):
 job=os.environ.get('SLURM_JOB_ID','');assert job.isdigit(),'Actual Slurm required for all scientific execution/tests'
 allocation=subprocess.check_output(['scontrol','show','job',job,'-o'],text=True);assert 'JobId='+job+' ' in allocation
 cfg=read(config);check(cfg['driver']);check(cfg['python']);check(cfg['plan'])
 for p in cfg['inputs']:check(p)
 out=Path(cfg['partial_test_output'] if partial else cfg['output']);out.mkdir(parents=True,exist_ok=False)
 t=time.monotonic();cpu=time.process_time();receipt={'config':pin(config),'job':job,'host':socket.gethostname(),'allocation':allocation,'new_molecular_calls':0,'mode':'actual_partial_parser_test' if partial else 'final_dispatch_join'}
 try:
  result=compare(cfg,partial)
  if partial:
   tests={'one_actual_accepted_cell':result['available_cells']==1,'seven_actual_unavailable':sum(c['status']=='unavailable' for c in result['cells'].values())==7,
    'no_invented_U':all(v['status']=='unavailable' and v['energy_hartree'] is None for v in result['DFT_CPCM']['U_water_minus_vacuum'].values()),
    'no_invented_D':result['DFT_CPCM']['D_sample1_minus_sample3']['energy_hartree'] is None,
    'raw_SP_preserved':result['cells']['sample1_La_water']['energy_hartree']==next(r for r in read(check(cfg['partial_test_collection']))['cells'] if r['task_id']=='sample1_La_water')['energy_hartree'],
    'GFN_saved_transfer_algebra':abs(result['GFN2_ALPB']['D_sample1_minus_sample3']['energy_kcal_mol']-21.256628301423902)<1e-8}
   write(out/'TESTS.json',{'actual_fixture':cfg['partial_test_collection'],'tests':tests,'all_passed':all(tests.values()),'synthetic_scientific_data':False})
   assert all(tests.values()),tests
  write(out/'RESULT.json',result);receipt.update(status='complete',result=pin(out/'RESULT.json'),scientific_matrix_status=result['status'])
 except Exception as exc:receipt.update(status='failed',error=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc());raise
 finally:
  receipt.update(wall_seconds=time.monotonic()-t,process_CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);write(out/'RECEIPT.json',receipt)
 print(json.dumps({'status':receipt['status'],'result':receipt['result'],'scientific_matrix_status':receipt['scientific_matrix_status']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);s.add_parser('freeze').add_argument('--workspace',required=True)
 for name in ('partial-test','collect'):s.add_parser(name).add_argument('--config',required=True)
 a=p.parse_args();freeze(a.workspace) if a.op=='freeze' else execute(a.config,a.op=='partial-test')
