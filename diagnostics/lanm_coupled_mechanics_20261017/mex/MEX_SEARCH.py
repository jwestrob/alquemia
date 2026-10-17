#!/usr/bin/env python3
"""Mex all-physical-coordinate search adapter; root binds and submits on Slurm.

The immutable MEX_CARTESIAN evaluator remains separately qualified. This module
only supplies physical guards, shared-optimizer API and checkpoint integration.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import socket
import time
import traceback
import xml.etree.ElementTree as ET

HERE=Path(__file__).resolve().parent
PROTOCOL='mex_coupled_cartesian_search_v1'
SCOPE={'site_membership_distance_A':3.5,'water_membership_distance_A':3.5,
       'metals':['La','Dy'],'all_physical_atoms_movable':2049,
       'search_limits':{'accepted_steps':80,'evaluations':200}}


def require(ok,message):
    if not ok:raise ValueError(message)


def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path'])['sha256']==record['sha256'],'Changed pin: '+record['path'])
    return Path(record['path'])


def read(record):return json.loads(verify(record).read_text())


def write(path,value):
    with Path(path).open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def build_guard(ev,optimizer,out):
    """Read the already adopted ff19SB/TIP3P graph and equilibrium bond lengths."""
    np=ev.np;native=read(ev.old.export['atoms']);native_by={a['id']:a for a in native}
    require(len(native)==2045 and [a['index'] for a in native]==list(range(2045)),
            'Native MM particle ordering differs')
    atoms=[];protein=[];water=[]
    for i,a in enumerate(ev.old.cells['La']['physical']):
        item=dict(native_by.get(a['id'],{}),**a);atoms.append(item)
        if item.get('kind')=='protein_source':protein.append(i)
        if item.get('kind')=='retained_crystal_water' and item['element']=='O':water.append(i)
    require(len(atoms)==2049 and protein and len(water)==167,'Physical protein/water inventory differs')
    root=ET.fromstring(verify(ev.old.export['system']).read_text())
    require(len(root.find('Particles'))==2045,'Native XML particle denominator differs')
    forces=[f for f in root.find('Forces') if f.get('type')=='HarmonicBondForce']
    require(len(forces)==1,'Unique native covalent bond force required')
    native_ids=[a['id'] for a in native];bonds=[];equilibrium=[]
    for b in forces[0].find('Bonds'):
        a,c=int(b.get('p1')),int(b.get('p2'))
        bonds.append([ev.index[native_ids[a]],ev.index[native_ids[c]]]);equilibrium.append(float(b.get('d'))*10.)
    original=read(ev.old.export['bonds'])
    require({tuple(sorted((ev.ids[a],ev.ids[b]))) for a,b in bonds}=={tuple(sorted(p)) for p in original},
            'Native force bond graph differs from full physical topology')
    source=ev.old.cells['La']['x'];acceptors=[i for i in protein if atoms[i]['element'] in ('N','O')]
    target=ev.old.cells['La']['target_id'];spectators=ev.old.source['spectator_occupancy']
    specs=[]
    for sid,label in [(target,'EF3_target')]+[(a['id'],a['site']+'_Nd_spectator') for a in spectators]:
        index=ev.index[sid]
        donors=[i for i in acceptors if np.linalg.norm(source[i]-source[index])<=SCOPE['site_membership_distance_A']]
        waters=[i for i in water if np.linalg.norm(source[i]-source[index])<=SCOPE['water_membership_distance_A']]
        require(donors,'No source-derived protein N/O donors for '+sid)
        specs.append({'index':index,'label':label,'source_id':sid,'site_donor_indices':donors,
                      'coordinating_water_indices':waters,'site_donor_source_ids':[ev.ids[i] for i in donors],
                      'coordinating_water_source_ids':[ev.ids[i] for i in waters]})
    geometry=ev.loaded[3].geometry_module
    guard=optimizer.PhysicalGuard(atoms,bonds,equilibrium,source,geometry,protein,specs)
    checks={m:guard.check(ev.starts[m]) for m in ('La','Dy')}
    require(all(c['pass'] and c['state_compatible'] for c in checks.values()),'Existing terminal start violates broad physical/state gate')
    write(out/'GUARD.json',{'scope':SCOPE,'source_inputs':ev.old.cells['La']['mapping']['inputs'],
          'native_system':ev.old.export['system'],'native_bonds':ev.old.export['bonds'],
          'physical_source_ids':ev.ids,'protein_indices':protein,'water_oxygen_indices':water,
          'bond_pairs':bonds,'equilibrium_A':equilibrium,'state_specs':specs,
          'start_checks':checks,'optimizer_settings':optimizer.SETTINGS,
          'geometry_helper':pin(geometry.__file__),'new_molecular_calls':0})
    return guard,geometry


class SearchAdapter:
    def __init__(self,ev,metal):self.ev=ev;self.metal=metal
    @property
    def new_evaluations(self):return self.ev.new_evaluations
    @property
    def evaluation_requests(self):return self.ev.requests
    @staticmethod
    def pin(path):return pin(path)
    def evaluate(self,positions):
        value=self.ev.evaluate(positions,self.metal)
        value['positions_A']=self.ev.np.asarray(positions).tolist()
        value['provenance']=dict(value['provenance'],cartesian_request=self.ev.requests[-1])
        return value
    def difference(self,new,old):return self.ev.difference(new,old)
    def read_point(self,provenance):
        require(provenance['identity_sha256']==self.ev.identity_hash,'Checkpoint model differs')
        record=provenance['cartesian_request'];value=self.ev.restore_cartesian(record)
        request=read(record);require(request['state']['metal']==self.metal,'Checkpoint target differs')
        value['positions_A']=request['all_physical_coordinates_A']
        value['provenance']=provenance
        return value


def gate(ev,cfg):
    result=read(cfg['qualification_result'])
    require(result['status']=='passed' and result['stage']=='qualification'
            and result['all_declared_passed'] and result['identity_sha256']==ev.identity_hash,
            'Complete same-model Cartesian derivative qualification required')
    verify(result['implementation']);verify(ev.cfg['implementation'])
    require(result['implementation']['sha256']==ev.cfg['implementation']['sha256'] and read(result['config'])==ev.cfg,
            'Qualified implementation/config differ')
    require(len(result['declared_checks'])==10 and all(c['status']=='passed' for c in result['declared_checks']),
            'Qualification denominator differs')
    require(result['counters']['new_composite_evaluations_attempted']==34,'Expected finite qualification not executed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--stage',choices=('preflight','run'),required=True)
    parser.add_argument('--metal',choices=('La','Dy'),default='La')
    args=parser.parse_args();require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    cfg=json.loads(args.config.read_text())
    require(cfg['protocol']==PROTOCOL and cfg['scope']==SCOPE,'Search scope differs')
    verify(cfg['implementation'])
    require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Search source differs')
    verify(cfg['plan']);verify(cfg['specification'])
    tick=time.monotonic();out=args.output;out.mkdir(parents=True,exist_ok=False)
    result={'status':'running','implementation':pin(__file__),'config':pin(args.config),'stage':args.stage,
            'metal':args.metal,'slurm_job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname()}
    ev=None
    try:
        optimizer=load('mex_shared_cartesian_optimizer',verify(cfg['optimizer']))
        ecfg=read(cfg['evaluator_config']);verify(cfg['evaluator']);verify(ecfg['implementation'])
        require(ecfg['implementation']['sha256']==cfg['evaluator']['sha256'],'Evaluator source/config differ')
        # Preserve the qualified evaluator's exact identity while allowing this
        # driver's own byte-identical execution snapshot at another path.
        evaluator=load('mex_full_cartesian_evaluator',verify(ecfg['implementation']))
        ev=evaluator.MexEvaluator(ecfg,out/'evaluations')
        guard,geometry=build_guard(ev,optimizer,out);result['guard']=pin(out/'GUARD.json')
        if args.stage=='preflight':result['status']='passed'
        else:
            gate(ev,cfg);adapter=SearchAdapter(ev,args.metal)
            start=adapter.evaluate(ev.starts[args.metal])
            resume=read(cfg['checkpoint']) if cfg.get('checkpoint') is not None else None
            if resume is not None:
                # Retain exact rejected-but-complete evaluations across interruption;
                # importing their immutable receipts adds no molecular calls.
                for provenance in resume.get('evaluated_points',[]):adapter.read_point(provenance)
            search=optimizer.search(adapter,guard,start,out/'search',SCOPE['search_limits'],geometry.lbfgs_direction,resume=resume)
            write(out/'SEARCH.json',search);result.update(status=search['status'],search=pin(out/'SEARCH.json'))
    except Exception as exc:result.update(status='failed',failure=repr(exc),traceback=traceback.format_exc())
    finally:
        if ev is not None:result.update(counters=ev.counters,requests=ev.requests,new_points=ev.new_points,
                                       native_points=ev.old.points,native_requests=ev.old.requests)
        result['elapsed_seconds']=time.monotonic()-tick
        write(out/'RESULT.json',result)
    raise SystemExit(0 if result['status'] in ('passed','complete') else 1)


if __name__=='__main__':main()
