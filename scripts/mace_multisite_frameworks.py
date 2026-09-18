"""Reuse actual native multisite parameter receipts through a verified identity map."""
from __future__ import annotations
import argparse
import copy
from pathlib import Path
import shutil
import time
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from mace_density_inputs import audit
from mace_density_frameworks import key,validate
from mace_multisite_density import audit_alias
from mace_tinker_framework_solver import check_parameters,parse_parameters

PROTOCOL='source_multisite_AMOEBA2018_GK_verified_framework_import_v1'


def prepare(preparation,capability,plan,output):
    start=time.monotonic();cpu=time.process_time();audit(preparation)
    p=read_json(preparation);cr=read_json(capability);cm=read_json(verify(cr['manifest']))
    if not cr['complete'] or cr['protocol']!='source_multisite_background_Ca_AMOEBA2018_native_GK_capability_v2':
        raise InvalidArtifact('actual complete multisite native capability required')
    if set(p['cases'])!=set(cm['cases']):raise InvalidArtifact('framework/source case inventory differs')
    parent=read_json(verify(cm['native_source']));root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    impl=root/'implementation';impl.mkdir();paths={k:verify(v) for k,v in p['implementation'].items()}
    for n in ('mace_multisite_frameworks.py','mace_density_frameworks.py','mace_multisite_amoeba.py','mace_amoeba_capability.py','mace_tinker_framework_solver.py','mace_tinker_capability.py'):
        paths[n]=Path(__file__).with_name(n)
    for n,s in paths.items():shutil.copyfile(s,impl/n)
    m=dict(protocol=PROTOCOL,preparation=record(preparation),reference_solver=cm['native_source'],parent=record(capability),
        software=parent['software'],executable=parent['executable'],native_parameters=parent['native_parameters'],plan=record(plan),
        implementation={f.name:record(f) for f in impl.glob('*.py')},tasks=[],new_energy_calls=0,new_response_solves=0,
        new_parameter_reads=0,reused_native_parameter_reads=len(p['cases']),numerical_score=None,baseline_changed=False)
    for case,pin in p['cases'].items():
        c=read_json(verify(pin));physical=read_json(verify(c['normalized_global_preparation']));audit_alias(physical)
        r=next(r['audit'] for r in cr['rows'] if r['case_id']==case);old=read_json(verify(r['mapping']));nr=r['native_GK']
        if physical['identity_parent']!=r['source_preparation'] or cm['cases'][case]!=r['source_preparation']:
            raise InvalidArtifact('framework and source bridge derive from different preparations')
        before=read_json(verify(physical['identity_parent']))
        if old['physical_atoms']!=before['physical_atoms']:raise InvalidArtifact('actual native mapping has different physical system')
        if nr['source_parameters']!=m['native_parameters'] or nr['software']!=m['software']:
            raise InvalidArtifact('native parameters/backend differ')
        receipt=read_json(verify(nr['receipt']));text=verify(receipt['log']).read_text()
        if receipt['returncode'] or 'ALQUEMIA_PREFLIGHT_COMPLETE' not in text:raise InvalidArtifact('native parameter execution unavailable')
        params=parse_parameters(text)
        if params!=read_json(verify(nr['parameters'])):raise InvalidArtifact('native output and saved parameters differ')
        # The selected metal is absent from this actual parameterized framework;
        # only its ledger identity is aliased. Every native atom is unchanged.
        mapping=dict(physical_atoms=physical['physical_atoms'],system_atom_ids=old['system_atom_ids'],
            source_mapping=r['mapping'],source_preparation=c['normalized_global_preparation'],
            identity_parent=physical['identity_parent'],selected_metal_identity_alias=physical['selected_metal_identity_alias'],
            support_state=c['state'],imported_capability=record(capability))
        t={**nr['task'],'task_id':case+'_all_standard','case_id':case,'variant':'all_standard',
            'parameters':nr['parameters'],'preflight_receipt':nr['receipt'],'reused_from':record(capability)}
        check_parameters(params,t,mapping)
        d=root/'tasks'/case;d.mkdir(parents=True);write_new(d/'mapping.json',mapping);t['mapping']=record(d/'mapping.json')
        t['cache_key']=key(t,m);m['tasks'].append(t)
    m['preparation_receipt']=dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu)
    write_new(root/'manifest.json',m);return validate(root/'manifest.json')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('preparation','capability','plan','output'):p.add_argument('--'+n,required=True)
    a=p.parse_args();r=prepare(a.preparation,a.capability,a.plan,a.output)
    print({'status':'pass','reused_parameter_reads':r['reused_native_parameter_reads'],'new_scientific_calls':0})
