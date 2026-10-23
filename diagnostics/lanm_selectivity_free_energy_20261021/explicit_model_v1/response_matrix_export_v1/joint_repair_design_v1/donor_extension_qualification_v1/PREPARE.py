#!/usr/bin/env python3
"""Freeze the one declared donor numerical panel. File preparation only."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import shutil
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
REL = Path('lanm_selectivity_free_energy_20261021')
D = ROOT/'diagnostics'/REL
W = ROOT/'workspaces'/REL
JOINT = Path('explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1')
PROFILES = {'absent': None, 'zero': (0., 0.), 'drive': (.1, 0.),
            'stiffness': (0., .1), 'combined': (.1, .1)}
FIELD_T = [.7562812564509392, -.15594410634600103, -.6353897204369424]
FIELD_H = (.001, .0005)
TRANSLATION = [.37, -.21, .16]
STRAIN = [[1., .2, -.3], [.2, -.7, .4], [-.3, .4, .5]]

def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def verify(p):
    if pin(p['path']) != p:
        raise ValueError('Changed authority: '+p['path'])
    return Path(p['path'])

def read(p):
    return json.loads(verify(p).read_text())

def save(path, value):
    with Path(path).open('x') as f:
        f.write(json.dumps(value, indent=2, allow_nan=False)+'\n')
    return pin(path)

def require(value, message):
    if not value:
        raise ValueError(message)

def dot(a, b): return math.fsum(x*y for x,y in zip(a,b))
def norm(a): return math.sqrt(dot(a,a))
def mv(a, b): return [dot(row,b) for row in a]

def parse_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0]); atoms = []
    require(len(lines) == n+1, 'Literal XYZ line denominator')
    for i,line in enumerate(lines[1:],1):
        t = line.split()
        require(int(t[0]) == i, 'XYZ atom ordering')
        atoms.append({'index': i, 'name': t[1], 'position_A': list(map(float,t[2:5])),
                      'type': int(t[5]), 'bonds': list(map(int,t[6:])),
                      'coordinate_tokens': t[2:5]})
    require(all(math.isfinite(x) for a in atoms for x in a['position_A']), 'Finite XYZ')
    return atoms

def xyz_text(atoms, positions, title):
    lines = [f'{len(atoms)}  {title}']
    for i,(a,r) in enumerate(zip(atoms,positions),1):
        lines.append(f'{i:5d} {a["name"]:3s} '+ ' '.join(format(x,'.15f') for x in r)
                     +f' {a["type"]:5d}'+''.join(f' {b:5d}' for b in a['bonds']))
    return '\n'.join(lines)+'\n'

def rotation():
    n = [x/math.sqrt(6.) for x in (1.,2.,-1.)]
    c,s = math.cos(.23),math.sin(.23)
    cross = [[0.,-n[2],n[1]],[n[2],0.,-n[0]],[-n[1],n[0],0.]]
    return [[c*(i==j)+(1-c)*n[i]*n[j]+s*cross[i][j] for j in range(3)] for i in range(3)]

def prepare(out):
    start,cpu = time.monotonic(),time.process_time()
    out = Path(out).resolve(); out.mkdir(parents=True,exist_ok=False)
    sources = {
        'design': pin(D/JOINT/'donor_extension_design_v1/DECISION.md'),
        'design_sources': pin(D/JOINT/'donor_extension_design_v1/SOURCES.json'),
        'water_manifest': pin(W/JOINT/'qualification_v1/preparation_v3/MANIFEST.json'),
        'water_actual_review': pin(W/JOINT/'qualification_v1/independent_actual_review_v1/REVIEW.json'),
        'donor_manifest': pin(W/'explicit_model_v1/paired_donor_fit_v1/dry_acetate_component_derivatives_v1/preparation_v1/MANIFEST.json'),
        'donor_result': pin(W/'explicit_model_v1/paired_donor_fit_v1/dry_acetate_component_derivatives_v1/run_v1/execution/RESULT.json'),
        'donor_actual_review': pin(W/'explicit_model_v1/paired_donor_fit_v1/dry_acetate_component_derivatives_v1/independent_actual_review_v1/REVIEW.json'),
        'old_field_design': pin(D/'reference_chemistry_v1/dry_donor_response_design_v1/SPEC.md'),
        'old_field_inventory': pin(W/'reference_chemistry_v1/dry_donor_response_design_v1/SOURCE_INVENTORY.json'),
        'original_testgrad': pin(W/'explicit_model_v1/tinker_build/testgrad'),
        'native_parser': pin(D/'explicit_model_v1/RUN_PANEL.py'),
        'gradient_parser': pin(D/'explicit_model_v1/paired_water_transfer_v1/RUN.py'),
        'matrix_parser': pin(D/'explicit_model_v1/response_matrix_export_v1/scout/RUN.py'),
        'water_qualification_source': pin(D/JOINT/'qualification_v1/RUN_v2.py'),
    }
    require(sources['design']['sha256']=='fbd1ff8abba3eeb8b851ac27c988900a6d8afe48d74a3e03a13246146e52e7e9','Frozen donor design')
    require(sources['water_actual_review']['sha256']=='7be6786d9388753735885f32b4b7164de2397f1cf5c6bc9e5d128b12a0e4f16f','Qualified water authority')
    require(sources['donor_manifest']['sha256']=='45311d308493e78c3ecde8522e80e0835124235f4257ac6c96582517910d9715','Frozen donor inputs')
    wm,wr,dm,dr = [read(sources[k]) for k in ('water_manifest','water_actual_review','donor_manifest','donor_result')]
    sources['water_result'] = wr['result']; water_result = read(wr['result'])
    require(water_result['numerical_checks_passed'] and wr['zero_baseline_parity_passed']
            and wr['force_checks_passed'] and wr['export_checks_passed'],'Water numerical admission')
    require(dr['complete_cells']==16 and dr['error'] is None,'Donor sixteen-cell producer')
    require(read(sources['donor_actual_review'])['result']==sources['donor_result'],'Donor result/review join')
    rows = {}
    for metal in ('La','Dy'):
        c = dm['central'][metal]; native = read(c['source']['authority'])
        require(native['status']=='complete','Selected central native admission')
        rows[metal+'_central'] = {'kind':'donor_center','metal':metal,
            'input':c['source']['effective_input'],'baseline':{
                'authority':c['source']['authority'],'energy':native['energy'],
                'forces':native['forces'],'outputs':c['source']['native_outputs']},
            'rotation_frame':c['analytic']['frame'],'q_A':0.}
    for cid,r in dm['rows'].items():
        if r['direction']!='rotation': continue
        old = dr['rows'][cid]
        require(old['status']=='complete','Endpoint native admission')
        for k in ('xyz','key','parameters'):
            require(old['input'][k]['sha256']==r['input'][k]['sha256'],'Endpoint input join')
        rows[cid] = {'kind':'donor_rotation','metal':r['metal'],'input':r['input'],
            'baseline':{'authority':sources['donor_result'],'row_id':cid,
                        'energy':old['parsed']['native']['total_energy_kcal_mol'],
                        'forces':None,'outputs':{'analyze_output':old['receipts'][0]['output']}},
            'rotation_frame':dm['central'][r['metal']]['analytic']['frame'],'q_A':r['q_A']}
    require(len(rows)==10,'Exactly ten consumed donor geometries')
    for cid in ('ion_La','ion_Dy','aquo9_central_La__near__two_water',
                'aquo9_central_La__dimer_01','aquo9_central_Dy__dimer_01'):
        old = wm['rows'][cid]; active = '__dimer_' in cid
        if active:
            r = water_result['rows']['combined__'+cid]
            require(r['status']=='complete','Active water gradient admission')
            inp = r['runtime_input']; b = r['value']
            require('POLAR-LOCAL-RESPONSE 500 36 0.050000000000000003' in verify(inp['key']).read_text(), 'Actual nonzero qualified water profile')
            authority = sources['water_result']; rowid='combined__'+cid
        else:
            inp = old['baseline_gradient_input']; b = old['baseline_gradient']
            authority = sources['water_manifest']; rowid=cid
        rows[cid] = {'kind':'unaffected_active_water' if active else 'unaffected',
            'metal':old['metal'],'input':inp,'water_profile_active':active,
            'baseline':{'authority':authority,'row_id':rowid,'energy':b['energy'],
                'forces':b['forces'],'outputs':{'gradient_output':b['output']}}}
    # Literal deletion, unchanged ligand coordinates and topology. No relaxation.
    parent = rows['La_central']; atoms = parse_xyz(verify(parent['input']['xyz']))
    require(atoms[0]['type']==500 and not atoms[0]['bonds'],'Unbonded metal deletion')
    ligand = []
    for a in atoms[1:]:
        require(1 not in a['bonds'],'Ligand has no metal bonds')
        ligand.append(dict(a,index=a['index']-1,bonds=[b-1 for b in a['bonds']]))
    derived = out/'derived'; derived.mkdir()
    lines=['7  literal_acetate7_from_consumed_La_center']
    for a in ligand:
        lines.append(f'{a["index"]:5d} {a["name"]:3s} '+' '.join(a['coordinate_tokens'])
                     +f' {a["type"]:5d}'+''.join(f' {b:5d}' for b in a['bonds']))
    (derived/'acetate7.xyz').write_text('\n'.join(lines)+'\n')
    rows['acetate7']={'kind':'unaffected_literal_deletion','metal':None,
        'input':dict(parent['input'],xyz=pin(derived/'acetate7.xyz')),
        'baseline':None,'parent':parent['input']['xyz'],
        'parent_atom_indices_0based':list(range(1,8)),
        'baseline_to_be_measured_by':'old__acetate7__gradient'}
    require(len(rows)==16,'Exactly sixteen base fixtures')
    for cid,row in rows.items():
        for p in row['input'].values(): verify(p)
        atoms=parse_xyz(row['input']['xyz']['path']); row['atoms']=len(atoms)
        row['types']=[a['type'] for a in atoms]
        row['positions_A']=[a['position_A'] for a in atoms]
        row['bonds_1based']=[a['bonds'] for a in atoms]
        row['names']=[a['name'] for a in atoms]
        if row['baseline']:
            for p in row['baseline']['outputs'].values(): verify(p)
        if row['kind'].startswith('donor'):
            require(row['types']==[500 if row['metal']=='La' else 501,401,402,403,403,404,404,404], 'Exact donor type order')
    units=wm['units']; factor=units['MV_cm_per_field_au']
    require(abs(norm(FIELD_T)-1)<2e-15,'Frozen unit field direction')
    programs=[]
    def add(cid,profile,program,stage,suffix='',fd=None,field=None,export=False,positions=None,transform=None,old=False):
        row=rows[cid]; jid='__'.join(filter(None,[profile,cid,suffix,program]))
        inpdir=out/'inputs'/jid; inpdir.mkdir(parents=True)
        src=row['input']; original=verify(src['xyz'])
        if positions is None: shutil.copyfile(original,inpdir/'system.xyz')
        else: (inpdir/'system.xyz').write_text(xyz_text(parse_xyz(original),positions,jid))
        shutil.copyfile(verify(src['parameters']),inpdir/'amoeba09.prm')
        key=verify(src['key']).read_text()
        require(not any(line.split() and line.split()[0].upper() in ('POLAR-DONOR-RESPONSE','EXTERNAL-FIELD','RESPONSE-MATRIX-EXPORT') for line in key.splitlines()),'Source key has no donor/export/field overlay')
        overlay=[]
        if not old:
            present={line.split()[0].upper() for line in key.splitlines() if line.split()}
            if 'POLAR-PRINT' not in present: overlay.append('POLAR-PRINT')
            if program=='analyze' and (export or field is not None or stage=='E') and 'SAVE-UINDUCE' not in present:
                overlay.append('SAVE-UINDUCE')
            if PROFILES[profile] is not None:
                for mt in (500,501): overlay.append(f'POLAR-DONOR-RESPONSE {mt} 403 '+' '.join(format(v,'.17g') for v in PROFILES[profile]))
            if export: overlay.append('RESPONSE-MATRIX-EXPORT')
            if field is not None: overlay.append('EXTERNAL-FIELD '+' '.join(format(factor*field*x,'.17g') for x in FIELD_T))
        (inpdir/'system.key').write_text(key+('\n'+'\n'.join(overlay)+'\n' if overlay else ''))
        inputs={k:pin(inpdir/name) for k,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]}
        args=['system.xyz',('EM' if export or field is not None or stage=='E' else 'E')] if program=='analyze' else ['system.xyz','Y','Y',format(fd,'.8f')] if fd else ['system.xyz','Y','N']
        programs.append({'id':jid,'cell':cid,'profile':profile,'program':program,'stage':stage,
            'binary_role':'original_testgrad' if old else program,'argv_tail':args,'input':inputs,
            'source_input':src,'key_overlay':overlay,'field_au':field,'exports':2 if export else 0,
            'fd_total_separation_A':fd,'fd_each_side_A':None if fd is None else fd/2,
            'internal_FD_energies':6*row['atoms'] if fd else 0,'transform':transform})
    donor=[cid for cid,r in rows.items() if r['kind'].startswith('donor')]
    controls=[cid for cid in rows if cid not in donor]
    add('acetate7','absent','testgrad','B',suffix='old',old=True)
    programs[-1]['id']='old__acetate7__gradient'
    for cid in donor:
        for profile in ('absent','zero'):
            central=rows[cid]['kind']=='donor_center'
            add(cid,profile,'analyze','A',export=central)
            if central: add(cid,profile,'testgrad','A')
    for cid in controls:
        for profile in ('absent','zero','combined'): add(cid,profile,'testgrad','B')
    for metal in ('La','Dy'):
        cid=metal+'_central'
        for profile in ('drive','stiffness','combined'):
            for label,sep in [('fd100',.0001),('fd050',.00005)]: add(cid,profile,'testgrad','C',suffix=label,fd=sep)
            add(cid,profile,'analyze','C',suffix='zero',export=True)
            for label,h in [('h',.001),('half',.0005)]:
                for sign in (-1,1):
                    tag=label+('_minus' if sign<0 else '_plus')
                    add(cid,profile,'analyze','D',suffix=tag,field=sign*h)
                    if profile=='combined': add(cid,profile,'testgrad','D',suffix=tag,field=sign*h)
    for cid in donor:
        if rows[cid]['kind']=='donor_rotation':
            add(cid,'combined','analyze','E');add(cid,'combined','testgrad','E')
    rot=rotation()
    for metal in ('La','Dy'):
        cid=metal+'_central'; pos=rows[cid]['positions_A']
        add(cid,'combined','testgrad','F',suffix='translation',positions=[[x+t for x,t in zip(r,TRANSLATION)] for r in pos],transform={'kind':'translation','vector_A':TRANSLATION})
        add(cid,'combined','testgrad','F',suffix='rotation',positions=[mv(rot,r) for r in pos],transform={'kind':'rotation','matrix':rot,'origin_A':[0.,0.,0.]})
        z=[[x-m for x,m in zip(r,pos[0])] for r in pos]; raw=[mv(STRAIN,r) for r in z]
        length=math.sqrt(math.fsum(dot(t,t) for t in raw));tangent=[[x/length for x in t] for t in raw]
        require(abs(math.fsum(dot(t,t) for t in tangent)-1)<1e-14,'Unit strain tangent')
        for label,h in [('h',.0001),('half',.00005)]:
            for sign in (-1,1):
                q=sign*h
                add(cid,'combined','analyze','G',suffix='strain_'+label+('_minus' if sign<0 else '_plus'),
                    positions=[[x+q*t for x,t in zip(r,v)] for r,v in zip(pos,tangent)],
                    transform={'kind':'strain','S':STRAIN,'normalization_L_A':length,'tangent':tangent,'q_A':q})
    counts=dict(Counter(p['stage'] for p in programs)); kinds=dict(Counter(p['program'] for p in programs))
    require(counts=={'A':24,'B':19,'C':18,'D':32,'E':16,'F':4,'G':8},'Declared stage counts')
    require(kinds=={'analyze':66,'testgrad':55} and len(programs)==121,'Declared program count')
    require(len({p['id'] for p in programs})==121,'Unique program IDs')
    require(sum(p['internal_FD_energies'] for p in programs)==576,'Internal FD count')
    require(sum(p['exports'] for p in programs)==20,'Exactly twenty snapshots')
    require(sum(p['binary_role']=='original_testgrad' for p in programs)==1,'One old-binary control')
    m={'status':'donor_numerical_inputs_prepared_no_execution','implementation':pin(__file__),
        'plan':pin(HERE/'QUALIFICATION_PLAN.md'),'sources':sources,'rows':rows,'programs':programs,
        'profiles':PROFILES,'units':units,'field_direction':FIELD_T,'origin_A':[0.,0.,0.],
        'stage_counts':counts,'program_counts':kinds,'maximum_native_programs':121,
        'maximum_internal_FD_energies':576,'maximum_workers':12,'maximum_snapshots':20,
        'native_source_manifest':None,'build_receipt':None,'executables':None,'root_release':None,
        'runner_implementation':None,'new_native_calls':0,'new_QM_calls':0,
        'physical_accuracy_qualified':False,
        'cost':{'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    return save(out/'MANIFEST.json',m)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);args=ap.parse_args()
    print(json.dumps(prepare(args.out)))
