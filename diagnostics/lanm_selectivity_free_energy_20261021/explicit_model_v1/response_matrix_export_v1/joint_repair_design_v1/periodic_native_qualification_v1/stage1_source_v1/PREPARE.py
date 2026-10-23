"""Literal Stage1 input/probe preparation only; no native, factor or operator call."""
import ast
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import time
import numpy as np


def pin(path):
    path=Path(path).resolve()
    with path.open('rb') as stream:
        h=hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(path=str(path),sha256=h)


def verify(item):
    if pin(item['path'])!=item:raise ValueError('Pin changed: '+item['path'])
    return Path(item['path'])


def load(item,name):
    source=verify(item)
    spec=importlib.util.spec_from_file_location(name,source)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod


def probes(source_item,n):
    # Import no factorization API: compile only the identical frozen pure helper.
    tree=ast.parse(verify(source_item).read_text())
    node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='ten_probes')
    namespace={'np':np}
    exec(compile(ast.Module(body=[node],type_ignores=[]),source_item['path'],'exec'),namespace)
    return namespace['ten_probes'](n,(0,1,2))


def prepare(authorities_path,out):
    tick,cpu=time.monotonic(),time.process_time()
    authorities_path=Path(authorities_path).resolve();auth=json.loads(authorities_path.read_text())
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    static=json.loads(verify(auth['static_input']).read_text())
    composition=json.loads(verify(static['composition_manifest']).read_text())
    if composition['candidate']['selected']!=static['candidate'] or not composition['selected_parameter_key'].startswith('e41fcd3d'):
        raise ValueError('Exact old candidate identity')
    native_build=json.loads(verify(auth['parent_build']).read_text())
    if native_build['manifest']!=auth['parent_source'] or native_build['status']!='complete':
        raise ValueError('Exact admitted build/source')
    library=native_build['library'];verify(library)
    row=next(x for x in composition['files'] if x['id']=='La_aqueous_BIO')['prepared']
    xyz=verify(row['xyz']).read_text();key=verify(row['key']).read_text();parameters=verify(row['parameters']).read_text()
    if row['xyz']!=static['geometry']:raise ValueError('Static geometry identity')
    shared=authorities_path.parent
    pme=shared.parent/'periodic_stability_definition_v1/pme_export_v1'
    derive=load(pin(pme/'full_box_initialization_v1/DERIVE.py'),'parity_input_derivation')
    sites_helper=load(pin(pme/'full_box_initialization_v1/EXECUTE_v2.py'),'parity_site_parser')
    old_key,provenance=derive.derive(key)
    # Preserve both exact original records; only add the already-qualified PME settings.
    actual_key=key+derive.SETTINGS
    removed=''.join(item['literal_line'] for item in provenance['removed'])
    if len(provenance['removed'])!=2 or not removed:raise ValueError('Literal response records')
    sites,cell=sites_helper.input_sites(xyz,parameters,actual_key)
    if cell!=[49.324]*3+[90.]*3:raise ValueError('Cell changed')
    for name,text in [('system.xyz',xyz),('system.key',actual_key),('amoebabio09.prm',parameters)]:
        (out/name).write_text(text)
    coefficients={}
    for item in provenance['removed']:
        fields=item['literal_line'].split();coefficients[int(fields[1])]=list(map(float,fields[3:]))
    if coefficients[900]!=static['metal_stiffness']['coefficients_a_b_c_d']:
        raise ValueError('La static coefficient match')
    with (out/'expected_state.dat').open('x') as stream:
        stream.write('PERIODIC_PARITY_EXPECTED_V1\n')
        stream.write(f"{len(sites)} {static['metal_stiffness']['normalized_beta']:.17g}\n")
        for metal in (900,901):stream.write(' '.join(format(x,'.17g') for x in coefficients[metal])+'\n')
        for site in sites:
            stream.write(f"{site['atom']} {site['type']} "+' '.join(format(x,'.17g') for x in site['xyz']+[site['alpha']])+'\n')
        stream.write('END\n')
    ps=probes(auth['probe_source'],3*len(sites))
    with (out/'probes.dat').open('x') as stream:
        stream.write(f'PERIODIC_PARITY_PROBES_V1\n{len(sites)} 10\n')
        for j,(label,v) in enumerate(ps,1):
            stream.write(f'{j} {label}\n')
            for i,values in enumerate(v.reshape(-1,3),1):
                stream.write(str(i)+' '+' '.join(format(x,'.17g') for x in values)+'\n')
        stream.write('END\n')
    source=Path(__file__).parent
    record=dict(status='prepared_stage1_parity_inputs_no_native',implementation=pin(__file__),
                authorities=pin(authorities_path),static_input=auth['static_input'],
                static_root_admission=auth['static_root_admission'],candidate=static['candidate'],
                composition=static['composition_manifest'],parent_build=auth['parent_build'],library=library,
                raw_inputs={k:pin(out/k) for k in ('system.xyz','system.key','amoebabio09.prm','expected_state.dat','probes.dat')},
                original_inputs=row,settings=pin(pme/'full_box_initialization_v1/DERIVE.py'),
                input_site_helper=pin(pme/'full_box_initialization_v1/EXECUTE_v2.py'),
                snapshot_reader=pin(pme/'SNAPSHOT.py'),coefficients=coefficients,probe_source=auth['probe_source'],
                labels=[x[0] for x in ps],probe_dimension=3*len(sites),
                exact_composed_key_retained=True,appended_settings_sha256=hashlib.sha256(derive.SETTINGS.encode()).hexdigest(),
                operator_only_initialization_key_identity=provenance,
                numerical_contract=pin(shared/'NUMERICAL_CONTRACT.md'),
                nominal_contract=auth['nominal_contract'],root_intent=auth['root_intent'],
                source={name:pin(source/name) for name in ('PARITY.f90','PARITY_SETUP.f','SETUP_PROVENANCE.json')},
                native_calls=0,builds=0,operator_applications=0,factorizations=0,
                counts=dict(native_programs=1,complete_operator=30,separate_reciprocal=30,total_reciprocal=60,
                            direct_initializations=1,induced_solves=0,energies=0,gradients=0),
                cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
    (out/'MANIFEST.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


if __name__=='__main__':
    import sys
    result=prepare(*sys.argv[1:]);print(result['status'])
