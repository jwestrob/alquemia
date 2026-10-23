#!/usr/bin/env python3
"""Delete only water atoms; bind existing executors without molecular calls."""
import collections
import copy
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
OUT = BASE/'matched_water_removal_v1/prepared_v1'
HYDRATED = BASE/'reference_chemistry_v1/hydrated_mixed_v1/bound_v1'
SOURCE = HYDRATED/'NATIVE_PANEL.json'
MACE_CONFIG = HYDRATED/'MACE_CONFIG.json'
ADMISSION = BASE/'hydrated_comparison_v1/partial_v1/RESULT.json'
Z = {'H':1, 'C':6, 'N':7, 'O':8, 'La':57, 'Dy':66}


def pin(path):
    p=Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    if pin(record['path']) != record:
        raise ValueError('Changed artifact: '+record['path'])
    return Path(record['path'])


def write(path,data):
    with Path(path).open('x') as handle:
        json.dump(data,handle,indent=2,allow_nan=False);handle.write('\n')
    return pin(path)


def check(ok,message):
    if not ok:
        raise ValueError(message)


def main():
    start=time.monotonic()
    OUT.mkdir(parents=True,exist_ok=False);(OUT/'xyz').mkdir();(OUT/'preview_inputs').mkdir()
    parent=json.loads(SOURCE.read_text());mace=json.loads(MACE_CONFIG.read_text())
    admitted=json.loads(ADMISSION.read_text())
    check(admitted['available']=={'MACE':20,'AMOEBA_transfer':20,'ecp_pbe0':8},'Source endpoints unavailable')
    for p in admitted['inputs']:
        verify(p)
    for p in (parent['adapter'],mace['implementation'],mace['native_helper'],mace['native_config']):
        verify(p)
    check(parent['adapter']['sha256']=='eec0e60d982bda4e1a8685c6120f25c6ab74cc57e114cd3d8f370b4baccba129',
          'Reference executor changed')
    check(mace['implementation']['sha256']=='6f2b6bd6220c0bf251c329f1ef5445c3bbce280d5bb6da8b027f5ae9cb16f48a',
          'MACE executor changed')
    new_geometries=[];cells=[];tasks=[];checks=[]
    for donor in ('acetate_mono','acetate_bi','acetamide'):
        for metal in ('La','Dy'):
            old_id=f'hydrated_{donor}_central_{metal}'
            matches=[g for g in parent['geometries'] if g['geometry_id']==old_id]
            check(len(matches)==1,'Missing/duplicate parent geometry')
            old=matches[0]
            lines=verify(old['xyz']).read_text().splitlines()
            atom_lines=[line for line in lines[2:] if line.strip()]
            check(len(atom_lines)==old['atoms']==int(lines[0]),'Parent XYZ count')
            keep=[i for i in range(old['atoms']) if i==old['metal_index'] or i in old['ligand_indices']]
            deleted=[i for i in range(old['atoms']) if i not in keep]
            nw=old['water_count']
            check(collections.Counter(old['symbols'][i] for i in deleted)=={'O':nw,'H':2*nw},'Deleted inventory is not water only')
            check(nw==(6 if donor=='acetate_bi' else 7),'Wrong removed-water count')
            new_id=f'matched_dry_{donor}_central_{metal}'
            target=OUT/'xyz'/f'{new_id}.xyz'
            target.write_text(str(len(keep))+'\n'+new_id+'; exact central water deletion\n'+
                              '\n'.join(atom_lines[i] for i in keep)+'\n')
            check(target.read_text().splitlines()[2:]==[atom_lines[i] for i in keep], 'Surviving XYZ line changed')
            index={old_i:new_i for new_i,old_i in enumerate(keep)}
            g={'geometry_id':new_id,'case_id':f'matched_dry_{donor}','donor_family':donor,
               'distortion':'central','metal':metal,'metal_index':index[old['metal_index']],
               'charge':old['charge'],'physical_multiplicity':old['physical_multiplicity'],
               'atoms':len(keep),'symbols':[old['symbols'][i] for i in keep],
               'positions_A':[old['positions_A'][i] for i in keep],'xyz':pin(target),
               'source_geometry_id':old_id,'source_xyz':old['xyz'],
               'retained_source_indices':keep,'deleted_water_source_indices':deleted,
               'water_count':0,'removed_water_count':nw,
               'ligand_indices':[index[i] for i in old['ligand_indices']],
               'donor_indices':[index[i] for i in old['donor_indices']],
               'radial_direction':old['radial_direction'],'radial_shift_A':0.0,
               'source_tinker_types':[old['source_tinker_types'][i] for i in keep],
               'source_bonds_zero_based':[[index[j] for j in old['source_bonds_zero_based'][i] if j in index] for i in keep],
               'construction':'Only delete water atoms; preserve surviving XYZ text and order exactly.'}
            check(g['charge']==(3 if donor=='acetamide' else 2),'Charge change')
            check(g['physical_multiplicity']==(1 if metal=='La' else 6),'Physical state change')
            check(g['atoms']==(10 if donor=='acetamide' else 8),'Dry atom count')
            core=46 if metal=='La' else 55
            electrons=sum(Z[s] for s in g['symbols'])-g['charge']-core
            check(electrons==40,'Wrong explicit electron count')
            check((sum(Z[s] for s in g['symbols'])-g['charge']-g['physical_multiplicity']+1)%2==0,
                  'Physical electron/spin parity')
            new_geometries.append(g)
            cells.append({'cell_id':'ecp_pbe0__'+new_id,'recipe':'ecp_pbe0','geometry_id':new_id,
                          'charge':g['charge'],'multiplicity':1,'physical_multiplicity':g['physical_multiplicity'],
                          'metal':metal,'ecp_core':core,'explicit_electrons':electrons,'geometry':g,'pilot':True})
            tasks.append({k:g[k] for k in ('case_id','distortion','metal','metal_index','charge','physical_multiplicity','xyz','source_xyz')})
            tasks[-1]['id']=new_id
            checks.append({'geometry':new_id,'surviving_source_indices':keep,'removed_waters':nw,
                           'surviving_line_equality':True,'atoms':g['atoms'],'charge':g['charge'],
                           'physical_multiplicity':g['physical_multiplicity'],'effective_multiplicity':1,'explicit_electrons':electrons})
    check(len(new_geometries)==len(cells)==len(tasks)==6,'Finite denominator')
    for donor in ('acetate_mono','acetate_bi','acetamide'):
        pair=[g for g in new_geometries if g['donor_family']==donor]
        check(pair[0]['positions_A']==pair[1]['positions_A'] and pair[0]['symbols'][1:]==pair[1]['symbols'][1:],
              'La/Dy paired geometry differs')
    reuse=[p for p in admitted['inputs'] if p['path'].endswith('/RESULT.json')]
    inventory={'protocol':'matched_water_removal_v1','status':'prepared_no_molecular_calls',
               'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),
               'parent_native_panel':pin(SOURCE),'parent_mace_config':pin(MACE_CONFIG),
               'admitted_hydrated_comparison':pin(ADMISSION),'reused_result_pins':reuse,
               'geometries':new_geometries,'checks':checks,'native_cell_count':6,'MACE_evaluation_count':6,
               'new_molecular_calls':0,'local_threads':1,
               'comparison':'Cwet-Cdry=(E_Dy_wet-E_Dy_matched_dry)-(E_La_wet-E_La_matched_dry); same aquo8 cancels.',
               'interpretation':'Metal-dependent addition of fixed waters, not unique pair interaction or protein affinity.'}
    inventory_pin=write(OUT/'PANEL.json',inventory)
    native_panel=copy.deepcopy(parent)
    native_panel.update(geometries=new_geometries,cells=cells,free_fragments=[],
                        specification=pin(HERE/'SPEC.md'),source_inventory=inventory_pin,
                        parent_panel=pin(SOURCE),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                        pilot_cell_ids=[c['cell_id'] for c in cells])
    native_pin=write(OUT/'NATIVE_PANEL.json',native_panel)
    mace.update(specification=pin(HERE/'SPEC.md'),source_inventory=inventory_pin,
                maximum_calls=6,tasks=tasks)
    mace_pin=write(OUT/'MACE_CONFIG.json',mace)
    module_spec=importlib.util.spec_from_file_location('frozen_reference',verify(parent['adapter']))
    reference=importlib.util.module_from_spec(module_spec);module_spec.loader.exec_module(reference)
    validation=reference.validate(str(OUT/'NATIVE_PANEL.json'))
    for c in cells:
        (OUT/'preview_inputs'/f'{c["cell_id"]}.inp').write_text(reference.deck(native_panel,c,4,3000))
    preview_note={'ranks':4,'maxcore_MiB':3000,'purpose':'Static method preview only. Runtime MaxCore is derived by frozen executor from actual Slurm memory grant.'}
    write(OUT/'preview_inputs/PREVIEW_ONLY.json',preview_note)
    write(OUT/'REFERENCE_VALIDATION.json',validation)
    native_helper=importlib.util.spec_from_file_location('frozen_mace_reader',verify(mace['implementation']))
    mace_runner=importlib.util.module_from_spec(native_helper);native_helper.loader.exec_module(mace_runner)
    for task,g in zip(tasks,new_geometries):
        symbols,coords=mace_runner.xyz(task['xyz'])
        check(symbols==g['symbols'] and coords==g['positions_A'],'Frozen MACE coordinate reader mismatch')
    report={'status':'prepared_and_saved_validation_passed','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'inventory':inventory_pin,'native_panel':native_pin,'mace_config':mace_pin,
            'native_validation':pin(OUT/'REFERENCE_VALIDATION.json'),'preparation':pin(__file__),
            'source_rows':checks,'native_calls':0,'MACE_calls':0,'submissions':0,'threads':1,
            'elapsed_seconds':time.monotonic()-start,'root_review_required_before_submission':True}
    write(OUT/'RESULT.json',report)
    print(json.dumps(report))


if __name__=='__main__':
    main()
