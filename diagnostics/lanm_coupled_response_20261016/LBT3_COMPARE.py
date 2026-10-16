#!/usr/bin/env python3
"""Saved-only LBT3 representation/source comparison and ordered Mex exchanges."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time

EV_TO_KCAL = 23.06054783061903
SOURCES = ('7CCO','7CCN')
REPS = ('whole','matched_core')
METALS = ('La','Dy')
SITES = ('EF1','EF2','EF3')
Z = {'H':1,'C':6,'N':7,'O':8,'La':57,'Dy':66}
VERIFIED = set()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path=Path(path).resolve(); h=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1048576),b''):
            h.update(block)
    return {'path':str(path),'sha256':h.hexdigest()}


def verify(p):
    require(isinstance(p,dict) and isinstance(p.get('path'),str) and
            isinstance(p.get('sha256'),str),'invalid artifact pin schema')
    key=(p['path'],p['sha256'])
    if key not in VERIFIED:
        # Legacy pins can retain symlink spellings and extra metadata. Validate
        # the target's bytes; do not rewrite the recorded provenance dictionary.
        path=Path(p['path']).resolve(strict=True)
        require(path.is_file() and pin(path)['sha256']==p['sha256'],
                'missing/changed artifact: '+p['path'])
        VERIFIED.add(key)
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def write(path,value):
    with Path(path).open('x') as handle:
        json.dump(value,handle,indent=2,sort_keys=True,allow_nan=False); handle.write('\n')


def xyz(p):
    text=verify(p).read_text().splitlines()
    rows=[line.split() for line in text[2:] if line.strip()]
    require(len(rows)==int(text[0]) and all(len(a)==4 for a in rows),'invalid XYZ inventory')
    result=[(a[0],*[float(v) for v in a[1:]]) for a in rows]
    require(all(a[0] in Z and all(math.isfinite(v) for v in a[1:]) for a in result),'invalid element/coordinate')
    return result


def keyed(rows,key):
    result={r[key]:r for r in rows}
    require(len(result)==len(rows),'duplicate '+key)
    return result


def state_check(row, coords, q, mult, metal, mi, electrons):
    require(coords[mi][0]==metal and sum(a[0] in METALS for a in coords)==1,'target identity/count differs')
    require(sum(Z[a[0]] for a in coords)-q==electrons and (electrons-mult+1)%2==0,'electron/spin parity differs')
    batch=row['input_state_check']
    expected={'atoms':len(coords),'charge':q,'spin_multiplicity':mult,'head':'omol',
              'selected_metal_index':mi,'graph_count':1,'charge_prediction_claimed':False}
    require(all(batch.get(k)==v for k,v in expected.items()),'actual native batch state differs')
    require(batch.get('metal_neighbor_edge_count',0)>0,'no represented target neighbor edges')


def paired(a,b,mi):
    return len(a)==len(b) and all(x[1:]==y[1:] and (i==mi or x[0]==y[0]) for i,(x,y) in enumerate(zip(a,b)))


def number(value):
    return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)


def difference(positive,negative):
    if any(v is None for v in positive+negative):
        return {'status':'unavailable','eV':None,'kcal_mol':None,'reason':'required cell unavailable'}
    value=math.fsum(positive+[-v for v in negative])
    return {'status':'available','eV':value,'kcal_mol':value*EV_TO_KCAL}


def no_coordinates(value):
    if isinstance(value,dict):
        return {k:no_coordinates(v) for k,v in value.items() if k!='xyz_A'}
    if isinstance(value,list):
        return [no_coordinates(v) for v in value]
    return value


def preparations(prep,c):
    require(prep['status']=='complete' and len(prep['rows'])==2,'two completed source preparations required')
    rows=keyed(prep['rows'],'source_id'); require(set(rows)==set(SOURCES),'prepared sources differ')
    signatures={}; sources={}; flat=[]
    for source in SOURCES:
        row=rows[source]; require(row['status']=='prepared','source preparation unavailable')
        p=read(row['preparation']); sources[source]=p
        require(p['tasks']==row['tasks'] and p['source_id']==source,'prepared task lineage differs')
        flat.extend(row['tasks'])
        require(p['sequence']=='FIDTNNDGWIEGDELLA' and p['charge']==-2 and p['protein_charge']==-5,
                'whole-peptide chemical state differs')
        require(p['waters']==p['spectators']==p['caps']==[] and p['deposited_H_coordinates_preserved']
                and p['deposited_heavy_coordinates_preserved'],'prepared inventory or source preservation differs')
        for source_pin in p['source_pins'].values(): verify(source_pin)
        verify(p['native_peptide_system']); verify(p['forcefield']); verify(p['hydrogen_audit'])
        physical=read(p['physical_atoms']); native=read(p['native_peptide_atoms']); bonds=read(p['covalent_bonds'])
        require(len(physical)==256 and len(native)==255 and physical[-1]['element']=='M','whole atom inventory differs')
        signatures[(source,'whole')]={'atoms':[(a['id'],a['element'],a.get('resname'),a.get('kind')) for a in physical],
            'native_charges':[(a['id'],a['charge_e']) for a in native],
            'bonds':sorted((tuple(sorted((b['atom_a_id'],b['atom_b_id']))),b['forcefield_equilibrium_A']) for b in bonds),
            'terminal_chemistry':p['terminal_chemistry'],'protein_charge':p['protein_charge'],'charge':p['charge'],
            'forcefield':p['forcefield'],'pH':p['hydrogen_preparation']['pH']}
        mapping=read(p['endpoints']['matched_core']['La']['source_cap_map'])
        require(mapping['shared_real_atoms_exactly_equal'] and len(mapping['mapping'])==48,'matched source map differs')
        verify(mapping['legacy_builder'])
        signatures[(source,'matched_core')]={'map':no_coordinates(mapping['mapping']),
            'fragments':no_coordinates(mapping['fragments']),'charge':-1,'builder':mapping['legacy_builder']}
        for rep in REPS:
            for metal in METALS:
                t=p['endpoints'][rep][metal]; coords=xyz(t['xyz'])
                require(len(coords)==(256 if rep=='whole' else 48),'representation size differs')
                if rep=='whole':
                    wanted=[(metal if a['element']=='M' else a['element'],*a['xyz_A']) for a in physical]
                    require(coords==wanted,'whole coordinates differ from physical map')
                else:
                    for a in mapping['mapping']:
                        i=a['core_index']
                        if a['kind']=='source':
                            parent=physical[a['physical_index']]
                            require(a['physical_id']==parent['id'] and coords[i]==
                                    (metal if parent['element']=='M' else parent['element'],*parent['xyz_A']),
                                    'matched core real atom differs from whole peptide')
                        else:
                            require(a['kind']=='synthetic_cap' and coords[i]==('H',*a['xyz_A']),
                                    'matched core cap differs from explicit map')
    prepared=keyed(flat,'task_id')
    require(len(c['tasks'])==8 and set(prepared)=={t['task_id'] for t in c['tasks']},'eight-cell task denominator differs')
    for t in c['tasks']:
        require({k:v for k,v in t.items() if k!='cache_key'}==prepared[t['task_id']],'scored task differs from preparation')
    pools={rep:{'compatible':signatures[(SOURCES[0],rep)]==signatures[(SOURCES[1],rep)],
                'inventory_checks':'ordered identities/elements, chemical state, cap rules; whole native charges/bonds/termini',
                'source_preparations':[rows[s]['preparation'] for s in SOURCES]} for rep in REPS}
    return pools


def run(config_path,out):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'saved numerical analysis requires Slurm')
    start=time.monotonic(); cfg=json.loads(config_path.read_text()); verify(cfg['plan'])
    r=read(cfg['lbt3_result']); c=read(r['config']); prep=read(c['preparation'])
    verify(c['plan'])
    require(c['plan']['sha256']==cfg['plan']['sha256'] and r['plan']==c['plan'],'frozen plan differs')
    verify(c['implementation']); verify(c['python'])
    m=read(c['native_precedent']); old=read(c['native_reuse']); software=read(m['software'])
    require(old['manifest']==c['native_precedent'],'Mex result uses another manifest')
    require(c['model']==r.get('model')==m['model'],'native model/recipe differs')
    model=c['model']; verify(model['checkpoint'])
    require(model['family']=='MACE_OMOL' and model['dtype']=='float64' and model['head']=='omol'
            and model['solvent'] is None and model['execution_adapter'] is None and not model['extra_dispersion']
            and not model['periodic'] and model['spin_input_convention']=='multiplicity'
            and model['energy_component']=='MACE_OMOL_total_vacuum_energy','unsupported physical model recipe')
    require(r.get('model_parameters_unchanged') is True,'new model mutation check unavailable')
    require(r['versions']==c['expected_versions']==old['versions'],'actual software versions differ')
    require(software['checkpoint']==model['checkpoint'] and software['python']['sha256']==c['python']['sha256'],
            'archived software/checkpoint/Python differs')
    verify(software['python']); verify(software['requirements']); verify(software['backend_source_inventory'])
    for name in ('mace_omol.py','lanm_series_followup.py'):
        require(c['helper_pins'][name]==m['implementation'][name],'native helper lineage differs')
        verify(c['helper_pins'][name])
    require(r['new_solvent_calls']==r['new_QM_calls']==r['optimizer_steps']==0,'unexpected molecular scope')
    pools=preparations(prep,c)
    found=keyed(r['rows'],'task_id'); declared=keyed(c['tasks'],'task_id')
    require(set(found)<=set(declared),'undeclared new result row')
    expected={(s,p,z) for s in SOURCES for p in REPS for z in METALS}
    require({(t['source_id'],t['representation'],t['metal']) for t in c['tasks']}==expected,'new matrix differs')
    cells=[]; energy={}; coordinates={}
    for t in c['tasks']:
        key=(t['source_id'],t['representation'],t['metal']); row=found.get(t['task_id'])
        item={k:t[k] for k in ('task_id','source_id','representation','metal','xyz','charge','multiplicity','atom_count','metal_index')}
        item.update(status='unavailable',energy_eV=None)
        try:
            require(row is not None,'missing saved cell')
            require(row['status']=='complete' and row.get('call_attempted') is True,row.get('reason','cell not completed'))
            require(all(row.get(k)==v for k,v in t.items()),'result task metadata differs')
            require(number(row['energy_eV']),'nonfinite energy')
            a=xyz(t['xyz']); coordinates[key]=a
            require(len(a)==t['atom_count'] and t['charge']==(-2 if key[1]=='whole' else -1)
                    and t['multiplicity']==(1 if key[2]=='La' else 6),'new physical state differs')
            state_check(row,a,t['charge'],t['multiplicity'],t['metal'],t['metal_index'],t['all_electron_count'])
            item.update(status='available',energy_eV=row['energy_eV'],input_state_check=row['input_state_check'],
                        wall_seconds=row['wall_seconds'])
        except (ValueError,KeyError,TypeError,OSError) as exc:
            item['reason']=str(exc)
        energy[key]=item['energy_eV']; cells.append(item)
    for source in SOURCES:
        for rep in REPS:
            keys=[(source,rep,z) for z in METALS]
            if all(k in coordinates for k in keys):
                require(paired(*(coordinates[k] for k in keys),255 if rep=='whole' else 0),'new paired coordinates differ')
    old_found=keyed(old['rows'],'endpoint_id'); old_endpoints=keyed(m['endpoints'],'endpoint_id')
    mex=[]; mex_energy={}
    for site in SITES:
        for metal in METALS:
            tid=f'Mex_{site}__{metal}'; e=old_endpoints[tid]; row=old_found.get(tid)
            item={'endpoint_id':tid,'site':site,'metal':metal,'status':'unavailable','energy_eV':None,'xyz':e['xyz']}
            try:
                require(tid in old['selected_endpoints'] and row is not None,'missing saved Mex cell')
                require(row['status']=='complete' and row.get('model_call_attempted') is True,row.get('reason','Mex cell incomplete'))
                require(all(row.get(k)==e[k] for k in ('xyz','state','cache_key')),'Mex task/state lineage differs')
                verify(e['source_preparation']); a=xyz(e['xyz'])
                mult=1 if metal=='La' else 6
                require(e['case_id']==f'Mex_{site}' and e['metal']==metal and e['charge']==-1
                        and e['spin_multiplicity']==mult and e['state']['physical_multiplicity']==mult
                        and e['state']['charge']==-1 and e['state']['oxidation_state']==3 and len(a)==44,
                        'old Mex physical state/inventory differs')
                state_check(row,a,-1,mult,metal,0,e['state']['all_electron_count'])
                require(number(row['energy_eV']),'nonfinite Mex energy')
                item.update(status='available',energy_eV=row['energy_eV'],input_state_check=row['input_state_check'],
                            state=e['state'],source_preparation=e['source_preparation'],wall_seconds=row['wall_seconds'])
                coordinates[('Mex',site,metal)]=a
            except (ValueError,KeyError,TypeError,OSError) as exc:
                item['reason']=str(exc)
            mex_energy[(site,metal)]=item['energy_eV']; mex.append(item)
        keys=[('Mex',site,z) for z in METALS]
        if all(k in coordinates for k in keys): require(paired(*(coordinates[k] for k in keys),0),'Mex paired geometry differs')
    contrast=[]; shifts=[]; spread=[]; exchange=[]
    for source in SOURCES:
        for rep in REPS:
            contrast.append({'source_id':source,'representation':rep,
                **difference([energy[(source,rep,'Dy')]],[energy[(source,rep,'La')]])})
            vector=[difference([mex_energy[(site,'Dy')],energy[(source,rep,'La')]],
                               [mex_energy[(site,'La')],energy[(source,rep,'Dy')]]) for site in SITES]
            exchange.append({'source_id':source,'representation':rep,'site_order':list(SITES),'R':vector})
        shifts.append({'source_id':source,**difference(
            [energy[(source,'whole','Dy')],energy[(source,'matched_core','La')]],
            [energy[(source,'whole','La')],energy[(source,'matched_core','Dy')]])})
    for rep in REPS:
        spread.append({'representation':rep,'signed_order':'7CCN_minus_7CCO',**difference(
            [energy[('7CCN',rep,'Dy')],energy[('7CCO',rep,'La')]],
            [energy[('7CCN',rep,'La')],energy[('7CCO',rep,'Dy')]])})
        pool=pools[rep]
        if pool['compatible'] and all(energy[(s,rep,z)] is not None for s in SOURCES for z in METALS):
            minima={z:min(energy[(s,rep,z)] for s in SOURCES) for z in METALS}
            pool.update(status='available',minimum_energy_eV=minima,
                selected_sources={z:[s for s in SOURCES if energy[(s,rep,z)]==minima[z]] for z in METALS},
                contrast=difference([minima['Dy']],[minima['La']]),
                Mex_exchange=[difference([mex_energy[(site,'Dy')],minima['La']],
                                         [mex_energy[(site,'La')],minima['Dy']]) for site in SITES])
        else:
            pool.update(status='unavailable',reason='chemical inventories incompatible or a required cell unavailable')
        pool['interpretation']='mathematical two-source minimum for each metal; no population or affinity claim'
    receipts=cfg.get('execution_receipts',{})
    for p in receipts.values(): verify(p)
    result={'protocol_id':'lbt3_whole_matched_core_saved_comparison_v1','config':pin(config_path),'plan':cfg['plan'],
        'implementation':pin(__file__),'status':'complete' if all(v['status']=='available' for v in cells+mex) else 'partial',
        'new_molecular_calls':0,'model_imports':0,'eV_to_kcal_mol':EV_TO_KCAL,
        'sources':{'LBT3_result':cfg['lbt3_result'],'LBT3_config':r['config'],'preparation':c['preparation'],
                   'Mex_result':c['native_reuse'],'Mex_manifest':c['native_precedent'],'software':m['software']},
        'compatibility':{'model':model,'actual_versions':r['versions'],'native_state_checks':'per admitted cell',
            'Mex_representation':'three distinct consumed 44-atom capped Nd-source site cores; older MACE, not later embedded QM',
            'GFN_components_used':False},
        'LBT3_cells':cells,'Mex_cells':mex,'LBT3_Dy_minus_La':contrast,
        'whole_minus_matched_core_contrast':shifts,'source_contrast_spread':spread,
        'Mex_minus_LBT3_exchange':exchange,'common_two_source_pools':pools,
        'formulae':{'D':'E_Dy - E_La','representation_shift':'D_whole - D_matched_core',
                    'source_spread':'D_7CCN - D_7CCO','R_i':'D_Mex_i - D_LBT3'},
        'limits':['Conditional dry vacuum model; no affinity, threshold or independent site labels.',
                  'No site summation. Mex protein-level preference is a qualitative expectation only.',
                  'Whole/core difference combines omitted scaffold, D13/terminal charge and caps.',
                  'Free charged termini are assumed; main-text/supplement charge discrepancy remains.',
                  'Raw unlike-metal energies include elemental offsets; interpret balanced differences.',
                  'Both source geometries remain reported; no favorable-source fallback.'],
        'costs':{'new_LBT3_model_calls':r['new_model_calls_attempted'],'LBT3_application_elapsed_s':r['elapsed_seconds'],
            'LBT3_model_load_s':r['model_load_seconds'],'LBT3_allocation':r['allocation'],'LBT3_device':r['device'],
            'reused_Mex_cells':sum(v['status']=='available' for v in mex),'historical_result_all_calls':old['molecular_calls_attempted'],
            'historical_result_all_wall_s':old['wall_seconds'],'historical_cost_scope':'includes other original cells, not attributed to six Mex rows',
            'preparation_work':prep['preparation_work'],'execution_receipts':receipts,
            'comparison_elapsed_s':time.monotonic()-start,'comparison_job_id':os.environ['SLURM_JOB_ID']}}
    out.mkdir(parents=True,exist_ok=False); write(out/'RESULT.json',result)
    def fmt(v): return 'unavailable' if v is None else f'{v:.6f}'
    lines=['# LBT3 whole-scaffold saved comparison','',
        f"Status: {result['status']}. Eight declared LBT3 cells; six saved Mex cells. No new model calls.",'',
        'Raw energies are eV; differences below are kcal/mol. D = E(Dy) − E(La). Raw D is not an affinity.','',
        '| Source | Representation | La energy (eV) | Dy energy (eV) | D (kcal/mol) |',
        '|---|---|---:|---:|---:|']
    for x in contrast:
        s,p=x['source_id'],x['representation']
        lines.append(f"| {s} | {p} | {fmt(energy[(s,p,'La')])} | {fmt(energy[(s,p,'Dy')])} | {fmt(x['kcal_mol'])} |")
    lines += ['','| Balanced comparison | Δ (kcal/mol) |','|---|---:|']
    lines += [f"| {x['source_id']}: D whole − matched core | {fmt(x['kcal_mol'])} |" for x in shifts]
    lines += [f"| {x['representation']}: D 7CCN − 7CCO | {fmt(x['kcal_mol'])} |" for x in spread]
    lines += ['','Rᵢ = D(Mex EFᵢ) − D(LBT3). Positive is conditional relative Mex La direction; no EF labels or summation.','',
              '| LBT3 source | Representation | R₁ | R₂ | R₃ |','|---|---|---:|---:|---:|']
    for x in exchange:
        lines.append(f"| {x['source_id']} | {x['representation']} | "+' | '.join(fmt(v['kcal_mol']) for v in x['R'])+' |')
    lines += ['','| Mathematical pool | La-selected source | Dy-selected source | D (kcal/mol) |','|---|---|---|---:|']
    for rep,pool in pools.items():
        selected=pool.get('selected_sources',{})
        lines.append(f"| {rep} | {','.join(selected.get('La',[])) or 'unavailable'} | {','.join(selected.get('Dy',[])) or 'unavailable'} | {fmt(pool.get('contrast',{}).get('kcal_mol'))} |")
    lines += ['','Pools require matching chemical inventory and every cell; they are minima over two structures, not populations.',
        'Whole/core changes include scaffold, D13/terminal charge and caps. Dry vacuum, terminal-state and physical-spin limitations remain. '
        'This consumed peptide control does not establish LanM profiling accuracy.','',
        f"LBT3: {r['new_model_calls_attempted']} model calls, {r['elapsed_seconds']:.3f}s application, job {r['slurm_job_id']}; "
        f"six historical Mex cells reused. Preparation's missing-H minimizations are recorded separately; their internal force-call count is unknown.",
        f"Exact inputs, runtime receipts, availability reasons and all matrices: `{out/'RESULT.json'}`."]
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json'),'report':pin(out/'REPORT.md')}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); run(args.config.resolve(),args.output.resolve())
