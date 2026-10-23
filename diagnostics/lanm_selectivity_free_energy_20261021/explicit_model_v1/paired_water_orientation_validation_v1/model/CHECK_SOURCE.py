#!/usr/bin/env python3
"""Narrow saved input/guard controls; no molecular program or reserved outcome."""
import argparse
import ast
import copy
import importlib.util
import json
from pathlib import Path
import time

started=time.monotonic();cpu=time.process_time();here=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('water_orientation_source',here/'RUN.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
p=argparse.ArgumentParser();p.add_argument('--bound',required=True);p.add_argument('--out',required=True);args=p.parse_args()
bound=m.read(m.pin(args.bound));cfg=m.read(bound['config']);m.bound_inputs(cfg,bound);checks={}
def check(name,ok):
    checks[name]=bool(ok);assert ok,name
def rejects(name,fn):
    try:fn()
    except (AssertionError,ValueError):checks[name]=True
    else:raise AssertionError(name)
transfer=m.read(cfg['transfer_result']);old=m.read(transfer['bound']);body=m.module(cfg['body_helper'],'water_orientation_control_body')
typing=m.module(cfg['typing_helper'],'water_orientation_control_typing')
for gid,row in bound['rows'].items():
    g=row['geometry'];lines=m.verify(row['files']['xyz']).read_text().splitlines();parent=row['parent_training_id'];meta=row['input']
    check('atoms_'+gid,int(lines[0].split()[0])==g['atoms']==25)
    check('coordinates_'+gid,[[float(x) for x in line.split()[2:5]] for line in lines[1:]]==g['positions_A'])
    check('types_'+gid,[int(line.split()[5]) for line in lines[1:]]==meta['types'])
    check('bonds_'+gid,[[int(x)-1 for x in line.split()[6:]] for line in lines[1:]]==meta['bonds_zero_based'])
    check('water_only_types_'+gid,set(meta['types'])=={36,37,500 if g['metal']=='La' else 501})
    check('selected_key_'+gid,row['files']['key']['sha256']==old['inputs'][parent]['files']['key']['sha256'])
    check('unchanged_parameters_'+gid,row['files']['parameters']['sha256']==old['inputs'][parent]['files']['parameters']['sha256'])
    check('metadata_files_'+gid,all(meta[k]==v for k,v in row['files'].items()))
    check('body_graph_'+gid,body.bodies(row)==row['bodies'])
    check('no_FD_'+gid,meta['derivative_test'] is False)
    check('state_'+gid,g['charge']==3 and g['physical_multiplicity']==(1 if g['metal']=='La' else 6))
    peer=gid[:-2]+('Dy' if gid.endswith('La') else 'La')
    check('paired_coordinates_'+gid,g['positions_A']==bound['rows'][peer]['geometry']['positions_A'])
    check('anchor_'+gid,bound['anchor_model'][parent]==transfer['rows'][parent]['model'])
check('four_exact_cells',set(bound['rows'])==set(m.CELLS) and len(m.CELLS)==4)
check('finite_eight_programs',cfg['maximum_native_programs']==2*len(m.CELLS)==8)
check('four_workers',cfg['maximum_workers']==4)
check('zero_new_FD',cfg['new_FD_energies']==0)
check('all_eight_raw_replayed',bound['all_eight_selected_water_raw_replayed'])
check('no_reserved_outcomes',bound['reserved_outcomes_opened'] is False)
check('zero_bind_native',bound['new_native_calls']==0)
check('donor_failures_retained',bound['donor_transfer_failures_preserved'] is True)
check('full_validation_unavailable',bound['full16_validation_complete'] is False)
for gid in m.ANCHORS:
    record=bound['aquo_cartesian_checks'][gid]
    check('existing_cartesian_'+gid,record==transfer['rows'][gid]['cartesian_check'] and record['fd_passed'])
bad=copy.deepcopy(bound);bad['parameters'][4]*=1.01
rejects('changed_water_candidate_rejected',lambda:m.bound_inputs(cfg,bad))
bad=copy.deepcopy(bound);bad['rows'].pop(next(iter(bad['rows'])))
rejects('missing_slot_rejected',lambda:m.bound_inputs(cfg,bad))
source=ast.parse((here/'RUN.py').read_text());run=next(n for n in source.body if isinstance(n,ast.FunctionDef) and n.name=='run')
runtime=next(n.value.args[0] for n in run.body if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and getattr(n.value.func,'id',None)=='require' and 'SLURM_JOB_ID' in ast.unparse(n.value))
gate=compile(ast.Expression(runtime),'<actual runtime condition>','eval')
class FakeOS:pass
fake=FakeOS();base={k:'1' for k in m.THREADS};base.update(SLURM_JOB_ID='1',SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='4');fake.environ=base
check('runtime_positive',eval(gate,{'os':fake,'THREADS':m.THREADS}))
for key,value in (('SLURM_JOB_ID',''),('SLURM_NTASKS','2'),('SLURM_CPUS_PER_TASK','8'),('OMP_NUM_THREADS','2')):
    fake.environ=dict(base);fake.environ[key]=value;check('runtime_rejects_'+key,not eval(gate,{'os':fake,'THREADS':m.THREADS}))
result={'status':'passed','implementation':m.pin(here/'RUN.py'),'control_source':m.pin(__file__),'configuration':bound['config'],'bound':m.pin(args.bound),
    'checks':checks,'check_count':len(checks),'native_calls':0,'new_FD_energies':0,'reserved_outcomes_opened':False,'donor_fit_result_opened':False,
    'local_threads':1,'elapsed_seconds':time.monotonic()-started,'CPU_seconds':time.process_time()-cpu}
out=Path(args.out);out.mkdir(parents=True,exist_ok=False);m.save(out/'RESULT.json',result)
print(json.dumps({'result':m.pin(out/'RESULT.json'),'checks':len(checks),'elapsed_seconds':result['elapsed_seconds']}))
