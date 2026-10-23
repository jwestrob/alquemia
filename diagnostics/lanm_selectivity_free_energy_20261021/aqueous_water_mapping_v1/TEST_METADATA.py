#!/usr/bin/env python3
"""Actual prepared-input and collision/override rejection controls; no energies."""
import argparse
import collections
import copy
import importlib.util
from pathlib import Path
import time

t=time.monotonic();cpu=time.process_time();here=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('aqueous_mapping',here/'PREPARE.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--out',required=True);a=p.parse_args()
mf=m.read(m.pin(a.manifest));cfg=m.read(mf['configuration']);inventory=m.module(cfg['inventory_helper'],'water_map_inventory');bio=inventory.parse(mf['BIO_parameters_preserved_byte_for_byte']);checks={}
def check(name,ok):checks[name]=bool(ok);assert ok,name
def reject(name,fn):
    try:fn()
    except (ValueError,AssertionError):checks[name]=True
    else:raise AssertionError(name)
check('exact_source',mf['implementation']==m.pin(here/'PREPARE.py'))
check('exact_candidate',mf['values']==m.VALUES)
check('finite_inventory',len(mf['cells'])==2 and len(mf['aqueous_endpoints'])==2)
for row in mf['cells']+mf['aqueous_endpoints']:
    for k,q in row['files'].items():m.verify(q);check('flat_file_'+row['id']+'_'+k,q==row[k])
    check('BIO_bytes_'+row['id'],row['parameters']['sha256']==mf['BIO_parameters_preserved_byte_for_byte']['sha256'])
    key=m.verify(row['key']).read_text();resolved=m.validate_key(key,mf['intrinsic_source_definitions'],bio)
    check('local_basename_'+row['id'],'parameters '+row['parameter_filename']+'\n' in key)
    check('mass_convention_'+row['id'],resolved['atom'][900]['mass']==resolved['atom'][901]['mass']==138.91)
    check('monopole_'+row['id'],resolved['metal_multipoles'][900]==resolved['metal_multipoles'][901]==[3.]+[0.]*9)
    check('no_metal_group_'+row['id'],not resolved['polarize'][900]['groups'] and not resolved['polarize'][901]['groups'])
    for typ in m.NONMETALS:
        check('BIO_alpha_'+row['id']+'_'+str(typ),resolved['polarize'][typ]['alpha_A3']==float(bio['polarize'][typ][0]['tokens'][2]))
for row in mf['cells']:
    lines=m.verify(row['xyz']).read_text().splitlines()[1:];old=m.verify(row['reference_xyz']).read_text().splitlines()[1:]
    check('coordinates_'+row['id'],[x.split()[2:5] for x in lines]==[x.split()[2:5] for x in old])
    check('bonds_'+row['id'],[x.split()[6:] for x in lines]==[x.split()[6:] for x in old])
    check('types_'+row['id'],[int(x.split()[5]) for x in lines]==row['types'])
    check('identity_map_'+row['id'],all(x['source_index1']==x['mapped_index1']==i+1 for i,x in enumerate(row['atom_map'])))
    check('saved_analytic_only_'+row['id'],row['reference_native']['derivative_test'] is False and row['si_model']['native']['finite_difference'] is None)
    for k in ('analyze_output','gradient_output'):m.verify(row['reference_native'][k]);check('SI_native_'+row['id']+'_'+k,row['reference_native'][k]==row['si_model']['native'][k])
for row in mf['aqueous_endpoints']:
    check('box_exact_'+row['id'],row['xyz']['sha256']==row['source']['sha256'])
    check('box_not_released_'+row['id'],row['prepared_only'] and not row['numerical_execution_authorized'])
    counts=collections.Counter(row['types']);check('neutral_inventory_'+row['id'],counts[405]==7 and counts[414]==10 and counts[402]==3980 and counts[403]==7960 and row['formal_charge']==0)
intrinsic=mf['intrinsic_source_definitions'];key=m.verify(mf['cells'][0]['key']).read_text()
for label,changed in [('water_H_instead_of_O',key.replace('polpair 402 900','polpair 403 900')),('wrong_class_instead_of_type',key.replace('vdwpr 73 900','vdwpr 402 900')),('donor_override',key+'vdwpr 31 900 3.6 0.2\n'),('salt_override',key+'polpair 414 901 0.1\n'),('oxygen_polarizability_override',key+'polarize 402 1.931 0.39 403\n'),('water_bonded_change',key+'ureybrad 74 73 74 -7.60 1.5337\n')]:
    reject(label,lambda changed=changed:m.validate_key(changed,intrinsic,bio))
occupied=copy.deepcopy(bio);occupied['atom'][900]=copy.deepcopy(bio['atom'][402]);reject('occupied_type',lambda:m.check_namespace(occupied))
occupied=copy.deepcopy(bio);occupied['atom'][402]['class']=901;reject('occupied_class',lambda:m.check_namespace(occupied))
row=mf['cells'][0];wrong=m.verify(row['reference_xyz']).read_text().replace(' 500',' 900');reject('already_mapped_or_foreign_source',lambda:m.map_aquo(wrong,'La'))
check('Urey_difference_preserved',mf['water_model_difference']['SI']['normalized_tokens'][-1]=='1.5337' and mf['water_model_difference']['BIO']['normalized_tokens'][-1]=='1.5537' and not mf['water_model_difference']['physical_difference_erased'])
check('no_donor_or_periodic_release',not mf['donor_parameters_bound'] and not mf['periodic_numerical_policy_bound'])
check('zero_native',mf['native_calls']==mf['energy_calls']==mf['force_calls']==mf['dynamics_steps']==0)
out=Path(a.out);out.mkdir(parents=True,exist_ok=False);result={'status':'metadata_controls_passed','implementation':mf['implementation'],'control_source':m.pin(__file__),'manifest':m.pin(a.manifest),'checks':checks,'check_count':len(checks),'native_calls':0,'reserved_outcomes_opened':False,'local_threads':1,'elapsed_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu};m.save(out/'RESULT.json',result);print(m.pin(out/'RESULT.json'))
