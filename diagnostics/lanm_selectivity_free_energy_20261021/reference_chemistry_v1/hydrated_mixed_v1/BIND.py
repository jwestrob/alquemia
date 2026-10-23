#!/usr/bin/env python3
"""Bind six native and eighteen MACE cells without importing model runtimes."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('reference',HERE.parent/'REFERENCE.py')
ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
MIXED=ref.ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/hydrated_mixed_v1/prepared_v2/PANEL.json'
MIXED_SHA='1191e3a00e3c02aa861faf83aae22eb250c94064e7ac780bf2fbebb00eb6c21a'
ORIGINAL=ref.ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/prepared_v2/PANEL.json'
MACE=ref.ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/mace_chemistry_v1/CONFIG.json'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();tick=time.monotonic()
    ref.check(ref.pin(MIXED)['sha256']==MIXED_SHA,'Mixed physical panel changed')
    physical=json.loads(MIXED.read_text());old=ref.load_panel(ORIGINAL)
    ref.check(old['adapter']['sha256']=='eec0e60d982bda4e1a8685c6120f25c6ab74cc57e114cd3d8f370b4baccba129','Native source changed')
    a.out.mkdir(parents=True,exist_ok=False);(a.out/'xyz').mkdir()
    geometries=[]
    for row in physical['geometries']:
        g=copy.deepcopy(row);src=ref.verify(g['xyz']);dst=a.out/'xyz'/src.name;shutil.copyfile(src,dst)
        g['source_xyz']=g['xyz'];g['xyz']=ref.pin(dst);geometries.append(g)
    central=[g for g in geometries if g['distortion']=='central'];cells=[]
    for g in central:
        core=46 if g['metal']=='La' else 55
        ne=sum(ref.Z[s] for s in g['symbols'])-g['charge']-core
        ref.check(ne%2==0,'Native effective-singlet electron parity')
        cells.append({'cell_id':'ecp_pbe0__'+g['geometry_id'],'recipe':'ecp_pbe0','geometry_id':g['geometry_id'],
            'charge':g['charge'],'multiplicity':1,'physical_multiplicity':g['physical_multiplicity'],
            'metal':g['metal'],'ecp_core':core,'explicit_electrons':ne,'geometry':g,'pilot':True})
    native={k:copy.deepcopy(old[k]) for k in ('protocol','adapter','published_coordinate_source','archived_recipe_source','orca','assets')}
    native.update(created_utc=ref.now(),specification=ref.pin(HERE/'RUN_SPEC.md'),parent_panel=ref.pin(ORIGINAL),
        source_inventory=ref.pin(MIXED),geometries=central,free_fragments=physical['free_ligands']+[physical['free_water']],
        cells=cells,pilot_cell_ids=[c['cell_id'] for c in cells],reused_aquo8=physical['reused_controls'],
        interpretation='Six new hydrated central ECP cells; original aquo8 results reused separately; no explicit-f expansion')
    ref.check(len(cells)==6 and len(geometries)==18,'Finite inventory differs')
    native_pin=ref.write(a.out/'NATIVE_PANEL.json',native)
    ref.write(a.out/'NATIVE_VALIDATION.json',ref.validate(native_pin['path']))
    cfg=json.loads(MACE.read_text())
    for key in ('implementation','native_helper','native_config'):ref.verify(cfg[key])
    cfg.update(specification=ref.pin(HERE/'RUN_SPEC.md'),source_inventory=ref.pin(MIXED),maximum_calls=18,
        source_config=ref.pin(MACE),tasks=[])
    for g in geometries:
        symbols,pos=ref.xyz_read(ref.verify(g['xyz']))
        ref.check(symbols==g['symbols'] and pos==g['positions_A'],'MACE physical fixture mismatch')
        ref.check(g['physical_multiplicity']==(1 if g['metal']=='La' else 6),'MACE spin mismatch')
        ref.check((sum(ref.Z[s] for s in symbols)-g['charge']-g['physical_multiplicity']+1)%2==0,'MACE parity')
        task={k:g[k] for k in ('case_id','distortion','metal','metal_index','charge','physical_multiplicity','xyz','source_xyz')}
        task['id']=g['geometry_id'];cfg['tasks'].append(task)
    mace_pin=ref.write(a.out/'MACE_CONFIG.json',cfg)
    receipt={'status':'ready_for_root_review','created_utc':ref.now(),'native_panel':native_pin,'mace_config':mace_pin,
        'new_native_cells':6,'new_mace_calls':18,'reused_controls_per_method':2,'source':ref.pin(__file__),
        'command':sys.argv,'local_numeric_threads':1,'elapsed_seconds':time.monotonic()-tick,'molecular_calls':0,
        'no_native_or_model_runtime_import':True,'aquo8_reuse':'Prior initial-panel results, never resubmit',
        'native_command_template':f"/usr/bin/python3 {old['adapter']['path']} run --panel {native_pin['path']} --cell CELL --out EXECUTION",
        'mace_command_template':f"PINNED_PYTHON {cfg['implementation']['path']} --config {mace_pin['path']} --output EXECUTION"}
    ref.write(a.out/'READY.json',receipt);print(json.dumps(receipt))

if __name__=='__main__':main()
