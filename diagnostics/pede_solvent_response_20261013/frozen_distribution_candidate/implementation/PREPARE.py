#!/usr/bin/env python3
"""Metadata-only freeze; no scientific imports, compilation, or evaluation."""
import hashlib
import json
from pathlib import Path
import shutil

OWN = Path(__file__).resolve().parent
ROOT = OWN.parents[3]
WORK = ROOT/'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1'
PRIOR = ROOT/'workspaces/pede_solvent_response_20261013/fixed_charge_replay_v3'
INVENTORY = ROOT/'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_inventory_v1'


def read(p):
    return json.loads(Path(p).read_text())


def pin(p):
    p = Path(p).resolve()
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return {'path':str(p),'sha256':h.hexdigest()}


def main():
    if (WORK/'manifest.json').exists():
        raise SystemExit('refusing to replace frozen manifest')
    old = read(PRIOR/'manifest.json')
    build = read(PRIOR/'execution_v1/BUILD.json')
    inventory = read(INVENTORY/'INVENTORY.json')
    models = [source['model'] for source in inventory['sources'].values()]
    assert all(model == models[0] for model in models)
    WORK.mkdir(parents=True,exist_ok=True)
    snapshot = WORK/'implementation'
    snapshot.mkdir()
    implementation = {}
    for name in ('NATIVE.py','REPLAY.f90','PREPARE.py','RUN.sbatch'):
        shutil.copyfile(OWN/name,snapshot/name)
        implementation[name] = pin(snapshot/name)
    base = PRIOR/'execution_v1'
    include = [base/'sources/xtb/include',base/'build/xtb/include',base/'sources/mctc-lib/include',base/'build/mctc-lib/include']
    files = sorted({p.resolve() for d in include for p in d.rglob('*') if p.is_file()})
    compiler = old['build']['tools']['gfortran']
    xtb = pin(base/'build/xtb/libxtb.a')
    mctc = pin(base/'build/mctc-lib/libmctc-lib.a')
    blas = old['build']['libraries']['blas']
    timer = pin('/usr/bin/time')
    m = {'protocol':'native_frozen_distribution_components_v1','status':'prepared_unexecuted',
         'inventory':pin(INVENTORY/'INVENTORY.json'),'reuse':pin(INVENTORY/'NATIVE_REUSE.json'),
         'implementation':implementation,'task_counts':inventory['summary']['unique_native_inputs_by_stage'],
         'candidate_energy':'E_OMOL,vac + G_nativeALPB(q_GFN2,vac)',
         'MACE_model':models[0],
         'native_model':old['native_model'],'tolerances':old['tolerances'],
         'required_setup_lines':old['required_setup_lines'],'HA_TO_KCAL':old['HA_TO_KCAL'],
         'provenance':{'plan':pin(OWN.parent/'PLAN.md'),'run_plan':pin(OWN.parent/'RUN_PLAN.md'),
                       'qualified_manifest':pin(PRIOR/'manifest.json'),
                       'qualified_collection':pin(base/'FINAL_COLLECTION.json'),
                       'qualified_build':pin(base/'BUILD.json'),
                       'qualified_wrapper':old['implementation']['REPLAY.f90']},
         'build':{'compiler':compiler,'xtb_static':xtb,'mctc_static':mctc,'blas':blas,'timer':timer,
                  'include_directories':[str(p) for p in include],
                  'library_directory':old['build']['library_directory'],
                  'linked_libraries':build['linked_libraries'],
                  'pins':[compiler,xtb,mctc,blas,timer,pin('/usr/bin/ldd'),*build['linked_libraries'],*[pin(p) for p in files]],
                  'engine_rebuilt':False,'installs':False,'network':False},
         'execution':{'canonical_native_inputs':150,'canonical_component_evaluations':300,
                      'transfer_native_inputs':1392,'transfer_reused_inputs':4,
                      'transfer_new_component_evaluations':2776,
                      'workers':'actual SLURM_NTASKS, range1..12','numerical_threads_per_worker':1,
                      'gate':'each water diagonal passes before its vacuum replay; canonical calibration passes before transfer',
                      'vacuum_precision':'same analytic rounding formula/native allowance and unchanged informative ceiling; declared before candidate evaluation',
                      'retry':'none; fresh output only; partial evidence retained',
                      'SCC_calls':0,'DFT_calls':0},
         'score_collection_owner':'root; no selection or calibration in this adapter',
         'production_changed':False}
    with (WORK/'manifest.json').open('x') as f:
        json.dump(m,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps(pin(WORK/'manifest.json')))


if __name__ == '__main__':
    main()
