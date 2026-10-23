#!/usr/bin/env python3
"""Self-contained saved-only matched hydration comparison, with explicit gaps."""
import argparse
from collections import Counter
import datetime
import hashlib
import json
import math
from pathlib import Path
import re
import time

EV=23.060547830619
HA=627.509474
BOHR=0.529177210903
Z={'H':1,'C':6,'N':7,'O':8,'La':57,'Dy':66}
METHODS=('MACE','AMOEBA_transfer','ecp_pbe0')
DONORS=('acetate_mono','acetate_bi','acetamide')
PINS={}
AMOEBA_INPUTS={}


def require(ok,message):
    if not ok: raise ValueError(message)


def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path'])==record,'Changed artifact: '+record['path'])
    PINS[record['path']]=record
    return Path(record['path'])


def load(record):
    return json.loads(verify(record).read_text())


def write(path,data):
    with Path(path).open('x') as f:
        json.dump(data,f,indent=2,allow_nan=False);f.write('\n')


def xyz(record):
    lines=verify(record).read_text().splitlines()
    rows=[s.split() for s in lines[2:] if s.strip()]
    require(len(rows)==int(lines[0]) and all(len(r)==4 for r in rows),'XYZ dimensions')
    return [r[0] for r in rows],[[float(x) for x in r[1:]] for r in rows]


def geometry(record,g):
    symbols,positions=xyz(record)
    require(symbols==g['symbols'] and positions==g['positions_A'],'Geometry differs: '+g['geometry_id'])


def flat(f): return [x for r in f for x in r]


def radial(force,g):
    return -math.fsum(force[i][a]*g['radial_direction'][a]
                      for i in g['ligand_indices'] for a in range(3))


def balanced(ids,signs,geometries):
    counts=Counter();charge=0
    for gid,sign in zip(ids,signs):
        g=geometries[gid];charge+=sign*g['charge']
        for symbol in g['symbols']: counts[symbol]+=sign
    require(charge==0 and all(n==0 for n in counts.values()),'Unbalanced exchange')


def cycle(data,method,ids,geometries):
    signs=[1,1,-1,-1]
    balanced(ids,signs,geometries)
    missing=[gid for gid in ids if gid not in data[method]]
    if missing: return {'status':'unavailable','missing_terms':missing,'value_kcal_mol':None}
    return {'status':'complete','missing_terms':[],
            'value_kcal_mol':math.fsum(s*data[method][gid]['energy_kcal_mol'] for s,gid in zip(signs,ids)),
            'terms':[{'geometry':gid,'coefficient':s,'energy_kcal_mol':data[method][gid]['energy_kcal_mol'],
                      'result':data[method][gid]['result']} for s,gid in zip(signs,ids)]}


def add(data,method,gid,energy,force,g,result_pin):
    require(len(force)==len(g['symbols']) and all(len(row)==3 for row in force),'Force dimensions')
    require(all(math.isfinite(v) for v in [energy]+flat(force)),'Nonfinite endpoint')
    value={'energy_kcal_mol':energy,'force_kcal_mol_A':force,'result':result_pin,
           'geometry':g['xyz'],'charge':g['charge'],'physical_multiplicity':g['physical_multiplicity']}
    if gid in data[method]:
        prior=data[method][gid]
        require(prior['energy_kcal_mol']==energy and prior['force_kcal_mol_A']==force,'Conflicting duplicate endpoint')
    else: data[method][gid]=value


def ingest_mace(record,geometries,data,expected_config=None,expected_ids=None):
    result=load(record);verify(result['implementation']);cfg=load(result['config'])
    if expected_config is not None: require(result['config']==expected_config,'MACE config binding differs')
    require(result['implementation']==cfg['implementation'],'MACE implementation mismatch')
    for key in ('specification','source_inventory','native_helper','native_config'): verify(cfg[key])
    require(result['energy_units']=='eV' and result['force_units']=='eV/Angstrom','MACE units differ')
    require(result['status'] in ('complete','incomplete','failed'),'Unknown MACE status')
    rows=result.get('rows',[]);ids=[r['id'] for r in rows]
    require(len(ids)==len(set(ids)),'Duplicate MACE rows')
    tasks={r['id']:r for r in cfg['tasks']}
    require(len(tasks)==cfg['maximum_calls'],'MACE planned count differs')
    if expected_ids is not None: require(set(tasks)==set(expected_ids),'MACE six-cell inventory differs')
    require(result.get('native_calls',0)<=len(tasks),'MACE exceeded finite calls')
    if result['status']=='complete':
        require(set(ids)==set(tasks) and result['native_calls']==len(tasks),'MACE complete denominator inconsistent')
    for row in rows:
        if row['id'] not in geometries: continue
        require(row['id'] in tasks and row['task']==tasks[row['id']],'MACE row/task binding')
        if row['status']!='complete': continue
        g=geometries[row['id']];geometry(row['task']['xyz'],g)
        require(row['task']['charge']==g['charge'] and
                row['task']['physical_multiplicity']==g['physical_multiplicity'],'MACE chemical state')
        add(data,'MACE',row['id'],row['energy_eV']*EV,[[x*EV for x in r] for r in row['forces_eV_A']],g,record)
    return {'status':result['status'],'row_statuses':{r['id']:r['status'] for r in rows},
            'config':result['config'],'model':result.get('model')}


def ingest_amoeba(record,geometries,data,expected_config=None,expected_ids=None):
    result=load(record);verify(result['implementation']);load(result['config'])
    if expected_config is not None: require(result['config']==expected_config,'AMOEBA config binding differs')
    rows=result.get('rows',[]);ids=[r['id'] for r in rows]
    require(len(ids)==len(set(ids)),'Duplicate AMOEBA rows')
    if expected_ids is not None:
        require(set(ids)<=set(expected_ids),'Unexpected AMOEBA row')
        if result['status']=='complete':
            require(set(ids)==set(expected_ids) and result['complete_rows']==6,'AMOEBA complete denominator')
    for row in rows:
        gid=row['id']
        if gid not in geometries or row['status']!='complete': continue
        g=geometries[gid];inp=row['input'];native=row['native']
        require(inp['charge']==g['charge'],'AMOEBA charge')
        geometry(inp['source_geometry']['xyz'],g)
        for key in ('key','parameters'): verify(inp[key])
        if gid.startswith('matched_dry_'):
            parent=AMOEBA_INPUTS[g['source_geometry_id']]
            require(inp['hydrated_parent_id']==g['source_geometry_id'],'Wrong AMOEBA hydrated parent')
            for key in ('key','parameters'):
                require(inp['parent_'+key]==parent[key], 'AMOEBA parent pin differs')
                verify(inp['parent_'+key])
                require(inp[key]['sha256']==parent[key]['sha256'],'AMOEBA water removal changes parameters/key')
            keep=g['retained_source_indices'];mapping={j:i for i,j in enumerate(keep)}
            require(inp['types']==[parent['types'][i] for i in keep], 'AMOEBA water removal changes types')
            require(inp['bonds_zero_based']==[[mapping[j] for j in parent['bonds_zero_based'][i] if j in mapping] for i in keep],
                    'AMOEBA water removal changes bonds')
        lines=verify(inp['xyz']).read_text().splitlines();atoms=[s.split() for s in lines[1:] if s.strip()]
        require(int(lines[0].split()[0])==len(atoms)==len(g['symbols']),'Tinker atom count')
        for i,a in enumerate(atoms):
            require(int(a[0])==i+1 and a[1]==g['symbols'][i] and
                    list(map(float,a[2:5]))==g['positions_A'][i],'Actual Tinker geometry')
            require(int(a[5])==inp['types'][i] and [int(v)-1 for v in a[6:]]==inp['bonds_zero_based'][i],
                    'Actual Tinker topology/types')
        grad_text=verify(native['gradient_output']).read_text()
        analyze=verify(native['analyze_output']).read_text()
        require(not re.search(r'(?i)not converged|fatal|\bNaN\b|\bInfinity\b',grad_text+analyze),'Native AMOEBA failure')
        energies=[]
        for text in (grad_text,analyze):
            found=re.findall(r'Total Potential Energy\s*:\s*([-+0-9.Ee]+)',text)
            require(len(found)==1,'Tinker energy count');energies.append(float(found[0]))
        require(energies[0]==energies[1]==native['total_energy_kcal_mol'],'Tinker energy replay')
        forces=[]
        for line in grad_text.splitlines():
            a=line.split()
            if len(a)==6 and a[0]=='Anlyt' and a[1].isdigit():
                require(int(a[1])==len(forces)+1,'Tinker gradient ordering')
                forces.append([-float(x) for x in a[2:5]])
        require(forces==native['force_kcal_mol_A'],'Tinker force replay')
        add(data,'AMOEBA_transfer',gid,energies[0],forces,g,record)
        AMOEBA_INPUTS[gid]=inp
    return {'status':result['status'],'row_statuses':{r['id']:r['status'] for r in rows},'config':result['config']}


def ingest_reference(record,panels,geometries,data,expected_cell=None):
    result=load(record)
    require(result['panel']['path'] in panels,'Reference panel outside admitted inventory')
    panel_record,panel=panels[result['panel']['path']]
    require(result['panel']==panel_record,'Reference panel changed')
    cells={c['cell_id']:c for c in panel['cells']}
    require(result['cell_id'] in cells,'Unknown reference cell')
    if expected_cell is not None: require(result['cell_id']==expected_cell,'Wrong new reference cell')
    cell=cells[result['cell_id']]
    if cell['recipe']!='ecp_pbe0' or cell['geometry_id'] not in geometries:
        return {'status':'outside_scoped_reference'}
    if result['status']!='complete': return {'status':result['status']}
    require(result['adapter']==panel['adapter'],'Reference adapter binding');verify(result['adapter'])
    require(result['checks'] and all(result['checks'].values()),'Contradictory reference admission')
    g=geometries[cell['geometry_id']]
    require(result['geometry']==cell['geometry']['xyz'],'Reference cell geometry pin')
    geometry(result['geometry'],g)
    for key in ('input','output'): verify(result[key])
    require(cell['charge']==g['charge'] and cell['physical_multiplicity']==g['physical_multiplicity'] and
            cell['multiplicity']==1 and cell['ecp_core']==(46 if g['metal']=='La' else 55),'Reference state differs')
    require(cell['explicit_electrons']==sum(Z[s] for s in g['symbols'])-g['charge']-cell['ecp_core'],
            'Reference explicit electron count')
    lines=[s.strip() for s in verify(result['engrad']).read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
    n=int(lines[0]);energy=float(lines[1]);require(n==len(g['symbols']) and len(lines)==2+4*n,'Engrad dimensions')
    grad=[float(s) for s in lines[2:2+3*n]]
    force=[[-grad[3*i+a]*HA/BOHR for a in range(3)] for i in range(n)]
    require(energy==result['energy_hartree'] and force==result['forces_kcal_mol_A'],'Native reference units/sign replay')
    for i,line in enumerate(lines[2+3*n:]):
        a=line.split();require(int(a[0])==Z[g['symbols'][i]],'Engrad element identity')
        require(max(abs(float(v)*BOHR-x) for v,x in zip(a[1:],g['positions_A'][i]))<1e-7,'Engrad coordinate precision')
    add(data,'ecp_pbe0',cell['geometry_id'],energy*HA,force,g,record)
    return {'status':'complete','cell':cell['cell_id'],'physical_state_qualified':result['physical_state_qualified'],
            'method_accuracy_qualified':result['method_accuracy_qualified'],'printed_scf_criteria':result['printed_scf_criteria']}


def evaluate(data,geometries):
    contrasts=[];hydration=[];placement=[];radials=[]
    for method in METHODS:
        for donor in DONORS:
            contexts={}
            for context,prefix in [('wet',f'hydrated_{donor}_central'),('matched_dry',f'matched_dry_{donor}_central'),
                                   ('original_dry',f'{donor}_central')]:
                ids=[prefix+'_Dy','aquo8_central_La',prefix+'_La','aquo8_central_Dy']
                value=cycle(data,method,ids,geometries)
                value.update(method=method,donor=donor,context=context);contrasts.append(value);contexts[context]=value
            for label,a,b,store in [('wet_minus_matched_dry','wet','matched_dry',hydration),
                                    ('matched_minus_original_dry','matched_dry','original_dry',placement)]:
                ar,br=contexts[a],contexts[b]
                entry={'method':method,'donor':donor,'contrast':label,'status':'unavailable','value_kcal_mol':None,
                       'missing_terms':sorted(set(ar['missing_terms']+br['missing_terms']))}
                if ar['status']==br['status']=='complete':
                    value=ar['value_kcal_mol']-br['value_kcal_mol']
                    a_ids=[ar['terms'][0]['geometry'],ar['terms'][2]['geometry']]
                    b_ids=[br['terms'][0]['geometry'],br['terms'][2]['geometry']]
                    direct_ids=[a_ids[0],b_ids[1],a_ids[1],b_ids[0]]
                    direct=cycle(data,method,direct_ids,geometries)
                    require(abs(value-direct['value_kcal_mol'])<2e-8,'Common-aquo8 cancellation failed')
                    entry.update(status='complete',value_kcal_mol=value,direct_cycle=direct,
                                 cancellation_residual_kcal_mol=value-direct['value_kcal_mol'])
                store.append(entry)
            ids=[f'hydrated_{donor}_central_{m}' for m in ('La','Dy')]+[
                 f'matched_dry_{donor}_central_{m}' for m in ('La','Dy')]
            missing=[i for i in ids if i not in data[method]]
            entry={'method':method,'donor':donor,'status':'unavailable','missing_terms':missing}
            if not missing:
                derivatives=[radial(data[method][i]['force_kcal_mol_A'],geometries[i]) for i in ids]
                wl,wd,dl,dd=derivatives
                entry.update(status='complete',wet_La_kcal_mol_A=wl,wet_Dy_kcal_mol_A=wd,
                             dry_La_kcal_mol_A=dl,dry_Dy_kcal_mol_A=dd,
                             water_addition_La_kcal_mol_A=wl-dl,water_addition_Dy_kcal_mol_A=wd-dd,
                             paired_water_addition_derivative_kcal_mol_A=(wd-dd)-(wl-dl))
            radials.append(entry)
    differences=[]
    for donor in DONORS:
        reference=next(r for r in hydration if r['method']=='ecp_pbe0' and r['donor']==donor)
        for method in ('MACE','AMOEBA_transfer'):
            model=next(r for r in hydration if r['method']==method and r['donor']==donor)
            available=reference['status']==model['status']=='complete'
            differences.append({'model':method,'reference':'ecp_pbe0','donor':donor,
                                'status':'complete' if available else 'unavailable',
                                'model_minus_reference_water_addition_kcal_mol':
                                model['value_kcal_mol']-reference['value_kcal_mol'] if available else None})
    return {'same_aquo8_exchange_contrasts':contrasts,'matched_water_addition_contrasts':hydration,
            'separate_original_dry_placement_contrasts':placement,'ligand_radial_response':radials,
            'water_addition_model_reference_differences':differences}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--out',required=True)
    args=ap.parse_args();start=time.monotonic();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
    result={'status':'started','new_molecular_calls':0,'local_threads':1,'config':pin(args.config),'implementation':pin(__file__)}
    try:
        cfg=load(result['config']);require(cfg['implementation']==pin(__file__),'Wrong collector');verify(cfg['specification'])
        dry=load(cfg['prepared_inventory']);dry_panel=load(cfg['native_panel']);wet_panel=load(cfg['wet_native_panel'])
        old=load(cfg['original_comparison']);wet=load(cfg['hydrated_comparison']);original_panel=load(old['panel'])
        mixed=load(wet_panel['source_inventory'])
        geometries={g['geometry_id']:g for g in dry['geometries']+mixed['geometries']+mixed['reused_controls']+original_panel['geometries']}
        needed={f'{prefix}{donor}_central_{metal}' for prefix in ('','hydrated_','matched_dry_')
                for donor in DONORS for metal in ('La','Dy')}|{'aquo8_central_La','aquo8_central_Dy'}
        geometries={k:v for k,v in geometries.items() if k in needed};require(set(geometries)==needed,'Geometry inventory differs')
        for g in geometries.values(): geometry(g['xyz'],g)
        for g in dry['geometries']:
            parent=geometries[g['source_geometry_id']];keep=g['retained_source_indices']
            require(g['symbols']==[parent['symbols'][i] for i in keep] and
                    g['positions_A']==[parent['positions_A'][i] for i in keep] and
                    g['charge']==parent['charge'] and g['physical_multiplicity']==parent['physical_multiplicity'],
                    'Water removal changes surviving geometry/state')
        panels={p['path']:(p,load(p)) for p in (cfg['native_panel'],cfg['wet_native_panel'],old['panel'])}
        data={m:{} for m in METHODS};source_status=[];model_metadata=[]
        # Cache denominators are admitted and explicitly pinned; only scope rows are read into data.
        cached_mace={old['mace']['path']:old['mace']};cached_amoeba={old['polarizable']['path']:old['polarizable']}
        cached_refs={p['path']:p for p in old['references']}
        for record in cfg.get('additional_original_references',[]):
            cached_refs[record['path']]=record
        for record in wet['inputs']:
            value=load(record)
            if 'force_units' in value and 'model' in value: cached_mace[record['path']]=record
            elif 'complete_rows' in value and 'rows' in value: cached_amoeba[record['path']]=record
            elif 'cell_id' in value: cached_refs[record['path']]=record
        for record in cached_mace.values():
            state=ingest_mace(record,geometries,data);model_metadata.append(state['model'])
        for record in cached_amoeba.values(): ingest_amoeba(record,geometries,data)
        for record in cached_refs.values():
            source_status.append({'result':record,**ingest_reference(record,panels,geometries,data)})
        dry_ids=[g['geometry_id'] for g in dry['geometries']];new_status=[]
        for method in ('MACE','AMOEBA_transfer'):
            binding=cfg['new_results'][method];path=Path(binding['path'])
            if not path.exists(): state={'status':'missing_result','row_statuses':{}}
            elif binding['config'] is None: state={'status':'awaiting_configuration_binding','row_statuses':{}}
            else:
                record=pin(path)
                state=(ingest_mace if method=='MACE' else ingest_amoeba)(record,geometries,data,binding['config'],dry_ids)
                state['result']=record
                if method=='MACE' and state.get('model') is not None: model_metadata.append(state['model'])
            new_status.append({'method':method,'path':str(path),**state})
        for binding in cfg['new_results']['ecp_pbe0']:
            path=Path(binding['path'])
            state=({'status':'missing_result'} if not path.exists() else
                   ingest_reference(pin(path),panels,geometries,data,binding['cell_id']))
            new_status.append({'method':'ecp_pbe0','path':str(path),'cell':binding['cell_id'],**state})
        require(all(m==model_metadata[0] for m in model_metadata),'Cached/new MACE Hamiltonian differs')
        slots=[{'method':method,'geometry':gid,'status':'admitted' if gid in data[method] else 'unavailable'}
               for method in METHODS for gid in dry_ids]
        result.update(evaluate(data,geometries))
        result.update(status='complete' if all(s['status']=='admitted' for s in slots) else 'partial',
                      new_endpoint_slots=slots,new_result_statuses=new_status,cached_reference_statuses=source_status,
                      endpoints=data,available_new={m:sum(g in data[m] for g in dry_ids) for m in METHODS},
                      available_all_scoped={m:len(v) for m,v in data.items()},
                      reviewed_pins=list(PINS.values()),physical_accuracy_qualified=False,
                      binding_free_energy_claimed=False,production_promoted=False,
                      units={'energy':'kcal/mol','force':'kcal/mol/Angstrom','MACE_eV_factor':EV,'Hartree_factor':HA,'bohr_A':BOHR})
        verify(result['config'])
    except Exception as exc:
        result.update(status='failed_validation',error_type=type(exc).__name__,error=str(exc),reviewed_pins=list(PINS.values()))
    result.update(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-start)
    write(out/'RESULT.json',result)
    lines=['# Matched water-addition comparison','',f'Status: **{result["status"]}**. Zero new molecular calls.','']
    for key,title in [('matched_water_addition_contrasts','Water addition at fixed metal–ligand geometry'),
                      ('separate_original_dry_placement_contrasts','Separate original-dry placement change')]:
        lines += ['## '+title,'','| Method | Donor | Status | Dy−La change, kcal/mol |','|---|---|---|---:|']
        for row in result.get(key,[]):
            value='—' if row['value_kcal_mol'] is None else f'{row["value_kcal_mol"]:.6f}'
            lines.append(f'| {row["method"]} | {row["donor"]} | {row["status"]} | {value} |')
        lines.append('')
    lines+=['Both sections use the same aquo8 reference. Missing results are unavailable, not zero.',
            'The first section is a total fixed-configuration many-body water-addition response.',
            'The second changes dry donor placement and is not part of the hydration-only contrast.',
            'Frozen-f/state/SCF/model/ensemble uncertainties remain; this is not binding free energy.']
    if 'error' in result: lines+=['','Validation error: '+result['error']]
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'status':result['status'],'available_new':result.get('available_new'),'result':str(out/'RESULT.json')}))
    return 1 if result['status']=='failed_validation' else 0


if __name__=='__main__': raise SystemExit(main())
