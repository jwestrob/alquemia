"""Liquid-water reference tied to real ORCA gas energy and NIST thermochemistry.

The harmonic ZPE approximation is explicit. A liquid-water reference alone does
not supply bound-state entropy, bound vibrational changes or missing cavity work.
"""
from __future__ import annotations
import argparse
import math
from pathlib import Path
import re
import shutil

from affordable_common import HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new, xyz
from hydration_square import endpoint

R_J_MOL_K = 8.31446261815324
CM_TO_J_MOL = 6.62607015e-34 * 299792458.0 * 100.0 * 6.02214076e23
PROTOCOL = 'native_r2scan3c_water_gas_to_liquid_298K_harmonic_reference_v1'


def thermochemical_terms(data):
    t = data['temperature_K']
    if t != 298.15 or data['standard_pressure_bar'] != 1.0:
        raise InvalidArtifact('tabulated H/S support only 298.15 K and 1 bar')
    a,b,c = data['vapor_pressure']['antoine_ABC']
    lo,hi = data['vapor_pressure']['temperature_range_K']
    if not lo <= t <= hi:
        raise InvalidArtifact('temperature outside Antoine fit')
    p = 10**(a-b/(t+c))
    zpe = .5*sum(data['harmonic_frequencies_cm_inverse'])*CM_TO_J_MOL/4184
    dh = data['gas_H298_minus_H0_kJ_mol']/4.184
    minus_ts = -t*data['gas_entropy_J_mol_K']/4184
    vapor = R_J_MOL_K*t*math.log(p/data['standard_pressure_bar'])/4184
    # Pure-liquid chemical potential at 1 bar, incompressible Poynting correction.
    pressure = data['liquid_molar_volume_cm3_mol']*1e-6*(data['liquid_pressure_bar']-p)*1e5/4184
    return {'harmonic_gas_ZPE_kcal_mol':zpe,'gas_thermal_enthalpy_increment_kcal_mol':dh,
            'minus_T_gas_entropy_kcal_mol':minus_ts,'vapor_pressure_bar':p,
            'gas_pressure_term_kcal_mol':vapor,'liquid_pressure_term_kcal_mol':pressure,
            'mu_liquid_minus_Egas_kcal_mol':zpe+dh+minus_ts+vapor+pressure,
            'mu_liquid_minus_Egas_minus_gas_ZPE_kcal_mol':dh+minus_ts+vapor+pressure,
            'approximation':'experimental harmonic ZPE plus tabulated gas H/S; ideal vapor, incompressible liquid',
            'uncertainty_kcal_mol':None}


def prepare(config, thermochemistry, output):
    cfg=read_json(config);source=read_json(verify(cfg['source_config']))
    data=read_json(thermochemistry);thermochemical_terms(data)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    coords=verify(source['water_reference_geometry'])
    if sorted(a[0] for a in xyz(coords))!=['H','H','O']:
        raise InvalidArtifact('not the pinned neutral H2O fixture')
    tasks=[]
    for name,keywords in [('gas','! r2SCAN-3c NoAutostart DefGrid3 TightSCF Opt'),
                          ('cpcm_comparator','! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF')]:
        d=out/name;d.mkdir();xp=d/'water.xyz';xp.write_bytes(coords.read_bytes())
        ip=d/'endpoint.inp';ip.write_text(keywords+'\n'+('%geom MaxIter 80 end\n' if name=='gas' else '')+'* xyzfile 0 1 water.xyz\n')
        tasks.append({'task_id':name,'input':record(ip),'xyz':record(xp),
                      'charge':0,'multiplicity':1,'output_path':str(d/'endpoint.out')})
    shutil.copyfile(__file__,out/'hydration_reference.py')
    shutil.copyfile(Path(__file__).with_name('hydration_square.py'),out/'hydration_square.py')
    m={'protocol_id':PROTOCOL,'configuration':record(config),'thermochemistry':record(thermochemistry),
       'agreement':cfg['agreement'],'orca':source['orca'],'execution_policy':source['execution_policy'],
       'tasks':tasks,'source_geometry':source['water_reference_geometry'],
       'execution_resources':{'mpi_ranks':1,'concurrent_tasks':1},
       'implementation':record(out/'hydration_reference.py'),'baseline_reference_changed':False}
    write_new(out/'manifest.json',m)
    return {'manifest':record(out/'manifest.json'),'tasks':len(tasks)}


def collect(manifest,output):
    m=read_json(manifest);terms=thermochemical_terms(read_json(verify(m['thermochemistry'])))
    results={}
    for t in m['tasks']:
        op=Path(t['output_path']);rp=Path(str(op)+'.execution.json')
        e=endpoint(record(op),record(rp),t['xyz'],t['input'])
        text=op.read_text()
        if t['task_id']=='gas' and 'THE OPTIMIZATION HAS CONVERGED' not in text:
            raise InvalidArtifact('gas water geometry not converged')
        if re.search(r'numerical gradient|numerical differentiation',text,re.I):
            raise InvalidArtifact('unexpected numerical derivative')
        e['scf_energy_evaluation_count']=len(re.findall('FINAL SINGLE POINT ENERGY',text))
        if t['task_id']=='gas':e['optimized_geometry']=record(op.parent/'endpoint.runtime.xyz')
        results[t['task_id']]=e
    e_gas=results['gas']['energy_hartree']
    result={'status':'approximate_bulk_reference_available','protocol_id':PROTOCOL,
            'manifest':record(manifest),'endpoints':results,'terms':terms,
            'mu_water_hartree':e_gas+terms['mu_liquid_minus_Egas_kcal_mol']/HA_TO_KCAL,
            'cpcm_comparator_minus_gas_minimum_kcal_mol':(results['cpcm_comparator']['energy_hartree']-e_gas)*HA_TO_KCAL,
            'bound_water_free_energy_corrections':None,'occupancy_probabilities':None,
            'limitations':['harmonic experimental gas ZPE is an approximation',
                           'bound-water basin entropy and vibrational shifts still unavailable',
                           'site CPCM lacks non-electrostatic cavity/dispersion terms',
                           'CPCM comparator uses original fixed water geometry; gas minimum optimized'],
            'production_reference_changed':False}
    write_new(output,result)
    return result


if __name__=='__main__':
    import json
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='operation',required=True)
    q=s.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--thermochemistry',required=True);q.add_argument('--output',required=True)
    q=s.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    args=vars(p.parse_args());op=args.pop('operation');print(json.dumps(globals()[op](**args),indent=2))
