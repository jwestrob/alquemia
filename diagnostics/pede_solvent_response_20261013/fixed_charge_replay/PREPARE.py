#!/usr/bin/env python3
"""Metadata-only freeze: copies no coordinates or charges and evaluates nothing."""
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[3]
OWN = Path(__file__).resolve().parent
WORK = ROOT/'workspaces/pede_solvent_response_20261013/fixed_charge_replay_v1'


def read(p):
    return json.loads(Path(p).read_text())


def pin(p):
    p = Path(p).resolve()
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}


def main():
    if (WORK/'manifest.json').exists():
        raise SystemExit('refusing to replace frozen manifest')
    invpath = ROOT/'diagnostics/pede_failure_diagnosis_20261013/model_review/ORIGIN_CHARGE_ARTIFACTS.json'
    resultpath = ROOT/'workspaces/pede_solvent_response_20261013/solvent_audit_v1/execution/RESULT.json'
    result = read(resultpath)
    inventory = read(invpath)
    energetics = read(result['inputs']['energetics']['path'])
    paritypath = ROOT/'workspaces/orca_parameter_audit_20261012/context_gfn2_v1/manifest.json'
    parity = read(paritypath)
    implementation = WORK/'implementation'
    implementation.mkdir()
    pins = {}
    for name in ('REPLAY.py','REPLAY.f90','CMakeLists.txt','RUN.sbatch','PREPARE.py'):
        shutil.copyfile(OWN/name,implementation/name)
        pins[name] = pin(implementation/name)
    cases = []
    for sample in (1,3):
        sid = f'q88jh5-pqq-la_model__conditioned_2Ca_v1__seed-1_sample-{sample}'
        archived = next(r for r in result['rows'] if r['source_case_id']==sid)
        for metal in ('Ca','La'):
            old = archived['metals'][metal]
            task = old['water']['task']
            item = {'case_id':f'sample{sample}_{metal}','source_case_id':sid,
                    **{k:task[k] for k in ('metal','charge','multiplicity','valence_electrons','atom_count','xyz')},
                    'archived_components_Eh':{k:old['terms_Eh'][k] for k in ('Gelec','Ghb','Gsasa','Gshift')},
                    'intrinsic_remainder_Eh':old['terms_Eh']['intrinsic_remainder'],
                    'U_Eh':old['U_Eh'],'endpoints':{}}
            for medium in ('water','vacuum'):
                row = next(r for r in inventory['rows'] if (r['source_case_id'],r['metal'],r['medium'])==(sid,metal,medium))
                b = row['backends']['xtb']
                item['state_signature'] = row['state_signature']
                item['endpoints'][medium] = {'point':row['result'],'charges':b['population_file'],
                    'stdout':b['stdout'],'execution':b['execution'],
                    'executed_xyz':pin(Path(b['stdout']['path']).parent/'core.xyz')}
            cases.append(item)
    sources = {}
    for name,version,commit,remote in (
        ('xtb','6.7.1','26b28010e805f7d1aeeef39813feb473e69cc4be','https://github.com/grimme-lab/xtb'),
        ('mctc-lib','0.4.2','51139755818af265803e600059cf3c44fe14995e','https://github.com/grimme-lab/mctc-lib'),
        ('test-drive','0.4.0','9c3401e30dbd2da1add77aaa252a4c6928fe39a1','https://github.com/fortran-lang/test-drive')):
        sources[name] = {'version':version,'commit':commit,'remote':remote,
                         'archive':pin(WORK/'upstream'/f'{name}-{version}.tar')}
    conversion = energetics['settings']['HA_TO_KCAL']
    m = {'protocol':'xtb671_native_fixed_charge_alpb_replay_v1',
         'status':'prepared_unexecuted','cases':cases,'component_evaluation_count':8,
         'execution_order':'four water-charge diagonals; only if all informative and pass, four vacuum-charge evaluations',
         'implementation':pins,'HA_TO_KCAL':conversion,
         'provenance':{'scope':pin(ROOT/'diagnostics/pede_solvent_response_20261013/FIXED_CHARGE_REPLAY_PLAN.md'),
                       'inventory':pin(invpath),'decomposition':pin(resultpath),
                       'gfN_energy_parity_precedent':pin(paritypath),
                       'runtime':pin(ROOT/'workspaces/orca_parameter_audit_20261012/context_gfn2_v1/RUNTIME.json')},
         'tolerances':{'charge_rounding_e':5e-9,'native_component_Eh':1e-9,
                       'native_allowance_source':'existing archived solvent-audit component-closure tolerance',
                       'maximum_informative_bound_kcal_mol':parity['settings']['energy_parity_Eh']*conversion,
                       'informative_ceiling_source':'existing context GFN2 backend energy parity1e-6Eh',
                       'Gelec_bound':'d*sum(abs(B_e*q)) + d*d/2*sum(abs(B_e)); B_e=native bornMat−diag(2*hbw)',
                       'Ghb_bound':'2*d*sum(abs(hbw*q)) + d*d*sum(abs(hbw))',
                       'Gsasa_Gshift_bound':'zero charge-rounding contribution; native allowance only',
                       'd':'half the stored 1e-8e rounding unit; no invented extra precision'},
         'native_model':{'version':'6.7.1','method_level':2,'solvent':'water','alpb':True,
                         'kernel':'native gbKernel%p16','state':'native solutionState%gsolv',
                         'temperature_K':298.15,'nAng':230,'ionStrength':0,'ionRad':0,
                         'parameter_source':'internal GFN2-xTB/ALPB factory; no manual radii or dielectric substitution',
                         'qsh':'zero-length unused argument; exact pinned getEnergyParts never references qsh'},
         'required_setup_lines':['internal GFN2-xTB/ALPB','gsolv [1 M gas/solution]',
                                 'Temperature                        2.9815E+02 K',
                                 'Interaction kernel             P16','Grid points                               230 per atom',
                                 'Ion screening                  false','H-bond correction              true'],
         'build':{'sources':sources,'tools':{k:pin('/usr/bin/'+v) for k,v in
                       (('cmake','cmake'),('gfortran','gfortran'),('gcc','gcc'),('make','make'))},
                  'libraries':{'blas':pin('/home/jwestrob/miniconda3/envs/lanm_bench/lib/libopenblasp-r0.3.32.so')},
                  'library_directory':'/home/jwestrob/miniconda3/envs/lanm_bench/lib',
                  'mctc_version_basis':'same0.4.2 dependency as installed xTB6.7.1 package',
                  'optional_backends':'tblite/CPCM-X disabled; native ALPB routine unchanged',
                  'installs':False,'network_at_build':False,'maximum_compile_parallelism':12,
                  'numerical_threads':1}}
    with (WORK/'manifest.json').open('x') as f:
        json.dump(m,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps(pin(WORK/'manifest.json')))


if __name__=='__main__':
    main()
