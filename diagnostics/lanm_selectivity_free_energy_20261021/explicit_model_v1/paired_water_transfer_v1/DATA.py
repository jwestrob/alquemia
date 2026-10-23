"""Exact consumed32 joins and water-only rendering. No native calls."""
import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
EP=HERE.parent/'paired_water_repair_v1/EVALUATE.py'
assert hashlib.sha256(EP.read_bytes()).hexdigest()=='72af9700048cc48e27ee1ad06bdd98cf9cefa974be083a286259e678d47d4003'
sp=importlib.util.spec_from_file_location('transfer_saved_util',EP);e=importlib.util.module_from_spec(sp);sp.loader.exec_module(e)
pin,verify,read,require,module=e.pin,e.verify,e.read,e.require,e.module
METALS=('La','Dy');FAMILIES=('acetate_mono','acetate_bi','acetamide');SIDES=('minus','central','plus')
CENTERS=tuple(q+'_central_'+m for q in ('water','aquo8') for m in METALS)
RADIAL=tuple('water_'+s+'_'+m for s in ('minus','plus') for m in METALS)
HYDRATED=tuple('hydrated_'+f+'_'+s+'_'+m for f in FAMILIES for s in SIDES for m in METALS)
DRY=tuple('matched_dry_'+f+'_central_'+m for f in FAMILIES for m in METALS)
CELLS=CENTERS+RADIAL+HYDRATED+DRY
FD=CENTERS+tuple(k for k in HYDRATED if '_central_' in k)
REF_HASH={'water_minus_La':'e889750771ca1db81a72bc02f492df9f82c2968e23bbe70745e7d63cfa2e1206',
 'water_plus_La':'404455607dcdfea0f1857f4da68cd6ba3ab4affcf4057275e6bde80c5e3c2ab8',
 'water_minus_Dy':'098d03cc2893463fd3c19a340f4bef6dadd0414cc9371d76c80a7e4912a659d9',
 'water_plus_Dy':'05dd6ce6305828bb17eae27e6c3712ce8427ae1b13f0c2db459766c6ae58b4f1'}

def default_config():
    old=e.default_config();b=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
    return {'data_source':pin(__file__),'specification':pin(HERE/'SPEC.md'),'evaluator_helper':pin(EP),
      'training':old['training_result'],'original_panel':pin(b/'explicit_model_v1/panel_v1/RESULT.json'),
      'radial_references':{k:pin(b/'reference_chemistry_v1/ecp_remaining_v1/tasks'/('ecp_pbe0__'+k)/'execution/RESULT.json') for k in RADIAL},
      'reference_helper':pin(HERE.parents[1]/'reference_chemistry_v1/REFERENCE.py'),
      'body_helper':pin(HERE.parent/'matched_dry_v1/OBJECTIVE.py'),'resolver':old['resolver'],
      'native_helper':old['native_helper'],'executables':old['executables'],
      'testgrad_source':pin(b/'explicit_model_v1/tinker_source/source/testgrad.f'),
      'initial_vector':list(e.INITIAL),'fixed_epsilons':e.EPSILON,
      'candidate_result':None,'root_release':None,'maximum_native_programs':48,'maximum_internal_FD_energies':1380}

def forces_ok(row,v):
    require(len(v)==row['geometry']['atoms'] and all(len(x)==3 and all(math.isfinite(t) for t in x) for x in v),'Invalid forces')

def prepare(cfg):
    require(cfg['data_source']==pin(__file__) and cfg['evaluator_helper']==pin(EP),'Data source differs')
    require(cfg['specification']['sha256']=='2dcecf34a7a5cf96beca2675effb956383a889fc502d03601ad4bfaf68bde78c','Frozen transfer SPEC changed')
    require(cfg['candidate_result'] is None and cfg['root_release'] is None,'Preparation is candidate-independent')
    require(cfg['maximum_native_programs']==48 and cfg['maximum_internal_FD_energies']==1380,'Finite ceiling differs')
    require(cfg['initial_vector']==list(e.INITIAL) and cfg['fixed_epsilons']==e.EPSILON,'Pair constants differ')
    require(cfg['training']['sha256']==e.TRAINING_SHA,'Wrong training authority')
    require(cfg['original_panel']['sha256']=='b256fdca923c0f476651e73250b93672eb49f7634fbceca66d2e6ea1d398a052','Wrong consumed panel')
    require(cfg['body_helper']['sha256']=='bf65a20a379811f8d35e2cb5ae78365a5d34022691b332bae5df223196078819','Body weights changed')
    require(cfg['testgrad_source']['sha256']=='5b3523b67e95441bbc0d45582b0055c263cbda13fd0a8b03dc81d8ce33924333','FD convention changed')
    for key in ('specification','resolver','native_helper','reference_helper','body_helper','testgrad_source'):verify(cfg[key])
    for p in cfg['executables'].values():verify(p)
    native=module(cfg['native_helper'],'transfer_native');body=module(cfg['body_helper'],'transfer_body');ref=module(cfg['reference_helper'],'transfer_reference')
    training=read(cfg['training']);require(training['status']=='training_complete' and training['complete']==28 and training['saved_controls_passed'],'Training incomplete')
    rows=copy.deepcopy(training['rows']);original=read(cfg['original_panel']);require(original['status']=='panel_complete','Consumed panel unavailable')
    for gid in RADIAL:
        require(cfg['radial_references'][gid]['sha256']==REF_HASH[gid],'Wrong radial QM result')
        old=next(x for x in original['rows'] if x['id']==gid);require(old['status']=='complete','Radial baseline missing')
        inp=old['input'];g=inp['source_geometry'];rr=read(cfg['radial_references'][gid])
        require(rr['status']=='complete' and rr['returncode']==0 and all(rr['checks'].values()),'Radial reference failed')
        require(rr['cell_id']=='ecp_pbe0__'+gid and rr['adapter']==cfg['reference_helper'],'Radial reference source differs')
        p=ref.load_panel(verify(rr['panel']));c=next(x for x in p['cells'] if x['cell_id']==rr['cell_id'])
        require(c['geometry']['symbols']==g['symbols'] and c['geometry']['positions_A']==g['positions_A'],'Radial geometry differs')
        for k in ('output','input','engrad','geometry'):verify(rr[k])
        replay=ref.collect_one(p,c,Path(rr['output']['path']).parent,0)
        require(all(rr[k]==v for k,v in replay.items()),'Radial raw reference differs')
        row={'status':'complete','input':inp,'geometry':g,'baseline_result':cfg['original_panel'],'reference_result':cfg['radial_references'][gid],
          'baseline':{'energy':old['native']['total_energy_kcal_mol'],'forces':old['native']['force_kcal_mol_A'],'native':old['native']},
          'reference':{'energy':rr['energy_hartree']*ref.HA,'energy_hartree':rr['energy_hartree'],'forces':rr['forces_kcal_mol_A'],
            **{k:rr[k] for k in ('output','engrad','input','panel')},
            'state':{k:c[k] for k in ('charge','multiplicity','physical_multiplicity','ecp_core','explicit_electrons')}}}
        row['bodies']=body.bodies(row);rows[gid]=row
    require(set(rows)==set(CELLS),'Exactly32 consumed cells required')
    for gid in CELLS:
        r=rows[gid];g=r['geometry'];i=r['input'];require(r['status']=='complete' and i['id']==g['geometry_id']==gid,'Cell identity')
        for k in ('xyz','key','parameters'):verify(i[k])
        verify(g['xyz']);verify(r['baseline_result']);rr=read(r['reference_result'])
        for k in ('input','output','engrad','panel'):verify(r['reference'][k])
        require(rr['status']=='complete' and rr['returncode']==0 and all(rr['checks'].values()),'Saved reference not admitted')
        require(rr['forces_kcal_mol_A']==r['reference']['forces'] and rr['energy_hartree']==r['reference']['energy_hartree'],'Reference arrays changed')
        require(rr['energy_hartree']*ref.HA==r['reference']['energy'],'Reference units changed')
        native_row=native.parse(verify(r['baseline']['native']['analyze_output']),verify(r['baseline']['native']['gradient_output']),i)
        require(native_row==r['baseline']['native'],'Raw baseline replay changed')
        require(native_row['total_energy_kcal_mol']==r['baseline']['energy'] and native_row['force_kcal_mol_A']==r['baseline']['forces'],'Baseline arrays changed')
        require(i['key']['sha256']=='6cdcc95d348953b8aafdb11b89cc03946eafee20eaa9a301622f8a8c3715cd89' and i['parameters']['sha256']=='e1cc2bcc6d7fb319aba7abaaf6b4483f2ef5ba730727212b5c157a18a290f38b','Original parameter model differs')
        lines=[x.split() for x in verify(i['xyz']).read_text().splitlines()[1:] if x.strip()]
        require(len(lines)==g['atoms']==i['atoms'] and [int(x[5]) for x in lines]==i['types'],'Native atom/type inventory')
        require([[int(v)-1 for v in x[6:]]for x in lines]==i['bonds_zero_based'],'Native bonds')
        require(max(abs(float(v[k+2])-g['positions_A'][j][k])for j,v in enumerate(lines)for k in range(3))<=5.1e-11,'Native geometry')
        q=2 if 'acetate_' in gid else 3;state=r['reference']['state'];metal=g['metal']
        require(i['charge']==g['charge']==state['charge']==q and state['multiplicity']==1 and state['ecp_core']==(46 if metal=='La' else 55),'Charge/ECP state')
        require(g['physical_multiplicity']==state['physical_multiplicity']==(1 if metal=='La' else 6),'Physical spin label')
        require(body.bodies(r)==r['bodies'],'Rigid body inventory')
        for src in ('baseline','reference'):forces_ok(r,r[src]['forces']);require(math.isfinite(r[src]['energy']),'Nonfinite energy')
    for q in sorted({k.rsplit('_',1)[0]for k in rows}):
        a,b=[rows[q+'_'+m] for m in METALS]
        require(a['geometry']['positions_A']==b['geometry']['positions_A'] and a['bodies']==b['bodies'] and a['input']['bonds_zero_based']==b['input']['bonds_zero_based'],'Paired coordinates/bodies differ')
    paths={}
    for family in ('water',)+tuple('hydrated_'+f for f in FAMILIES):
        for metal in METALS:
            center=rows[family+'_central_'+metal];g=center['geometry'];ids=[g['metal_index']]if family=='water'else g['ligand_indices'];u=g['radial_direction']
            require(abs(sum(v*v for v in u)-1)<1e-10,'Direction not unit')
            for side,h in (('minus',-.15),('central',0.),('plus',.15)):
                rg=rows[family+'_'+side+'_'+metal]['geometry'];require(rg['radial_direction']==u,'Direction changed')
                require(max(abs(x[k]-y[k]-(h*u[k]if j in ids else 0))for j,(x,y)in enumerate(zip(rg['positions_A'],g['positions_A']))for k in range(3))<2e-12,'Saved path does not match declared motion')
            paths[family+'_'+metal]={'moving_atoms':ids,'direction':u,'step_A':.15,'derivative':'minus sum(force of moving atoms) dot direction','small_step_qualification':False}
    require(sum(6*rows[k]['geometry']['atoms'] for k in FD)==1380,'FD energy count changed')
    return {'status':'consumed32_prepared','config':cfg,'rows':rows,'paths':paths,'inventory':{'fitted_centers':list(CENTERS),'water_radial':list(RADIAL),'hydrated':list(HYDRATED),'unchanged_dry':list(DRY),'full_cartesian':list(FD)},'maximum_native_programs':48,'maximum_internal_FD_energies':1380,'new_native_calls':0,'reserved_validation_opened':False}

def render(cfg,row,values):
    require(len(values)==4 and all(math.isfinite(v) for v in values),'Four finite values required')
    for j,v in enumerate(values):require((.05<=v<=.5)if j%2 else (.85*e.INITIAL[j]<=v<=1.15*e.INITIAL[j]),'Frozen parameter bounds')
    resolver=module(cfg['resolver'],'transfer_resolver');prm=verify(row['input']['parameters']).read_text();old=verify(row['input']['key']).read_text();before=resolver.parse_parameters(prm,old)
    mt=500 if row['geometry']['metal']=='La' else 501;offset=0 if mt==500 else 2
    require([t for t,a in before['atom'].items()if a['class']==34]==[36],'Unexpected water class alias')
    recipient=36 in row['input']['types'];require(recipient==(row['input']['id'] not in DRY),'Dry/wet recipient mismatch')
    original=not recipient or values[offset:offset+2]==list(e.INITIAL[offset:offset+2])
    metal=row['geometry']['metal']
    key=old if original else old.rstrip()+f'\n\n# Paired water development override: {metal} only\nvdwpr 34 {mt} {values[offset]:.17g} {e.EPSILON[metal]:.17g}\npolpair 36 {mt} {values[offset+1]:.17g}\n'
    after=resolver.parse_parameters(prm,key)
    for field in ('settings','atom','vdw','polarize','metal_multipoles'):require(before[field]==after[field],'Unrelated parameters changed')
    for field,allowed in (('vdwpr',tuple(sorted((34,mt)))),('polpair',tuple(sorted((36,mt))))):
        require({k:v for k,v in before[field].items()if k!=allowed}=={k:v for k,v in after[field].items()if k!=allowed},'Nonwater override changed')
    require(resolver.resolve_pair(before,mt,37)==resolver.resolve_pair(after,mt,37),'Water H changed')
    for t in set(row['input']['types'])-{mt,36}:require(resolver.resolve_pair(before,mt,t)==resolver.resolve_pair(after,mt,t),'Donor pair changed')
    pair=resolver.resolve_pair(after,mt,36)
    if recipient:require(pair['Rmin_pair_A']==values[offset] and pair['thole']==values[offset+1] and pair['epsilon_kcal_mol']==e.EPSILON[row['geometry']['metal']],'Rendered pair differs')
    return key,{'water_recipient':recipient,'original_key_reused':original,'resolved_water_pair':pair,'effective_values':values[offset:offset+2]if recipient else None}
