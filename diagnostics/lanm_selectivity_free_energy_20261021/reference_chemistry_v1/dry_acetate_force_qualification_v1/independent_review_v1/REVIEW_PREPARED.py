"""Independent saved input/raw-center checks; no molecular program imports."""
import hashlib
import json
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[5]
D = ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_force_qualification_v1'
W = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_force_qualification_v1'
PANEL = W/'preparation_v1/PANEL.json'
HA, BOHR = 627.509474, .529177210903
FLOAT = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?'
checks = []
seen = {}

def pin(p):
    p = Path(p).resolve()
    if str(p) not in seen:
        seen[str(p)] = {'path':str(p), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    return seen[str(p)]

def check(ok, name):
    if not ok:
        raise AssertionError(name)
    checks.append(name)

def verify(p):
    check(pin(p['path']) == p, 'pin '+p['path'])
    return Path(p['path'])

def sub(a,b): return [x-y for x,y in zip(a,b)]
def dot(a,b): return math.fsum(x*y for x,y in zip(a,b))
def norm(a): return math.sqrt(dot(a,a))
def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]

def main():
    p = json.loads(PANEL.read_text())
    verify(p['adapter']); verify(p['specification'])
    for helper in p['helpers'].values(): verify(helper)
    for metal in p['assets'].values():
        for asset in metal.values(): verify(asset)
    training = json.loads(verify(p['training_result']).read_text())
    check(p['maximum_native_calls']==16 and p['ranks']==4 and p['steps_A']==[.01,.005]
          and p['rotation_length_A']==1 and p['derivative_resolution_kcal_mol_A']==.1, 'scope and criteria')
    ids = {f'pbe0__{m}__{d}__{s}__complex' for m in ('La','Dy') for d in ('radial','rotation')
           for s in ('minus_010','minus_005','plus_005','plus_010')}
    check(len(p['cells'])==16 and {c['cell_id'] for c in p['cells']}==ids, 'exact sixteen')
    scout = {x for x in ids if '__La__radial__' in x}
    check(set(p['scout_cell_ids'])==scout and set(p['remaining_cell_ids'])==ids-scout, 'four then twelve')
    centers, directions, records = {}, {}, []
    for m in ('La','Dy'):
        old = p['central'][m]; g = old['original_cell']['geometry']; x = g['positions_A']
        for k in ('saved_result','saved_start','executed_input','output','engrad','original_panel','original_adapter'): verify(old[k])
        saved = json.loads(Path(old['saved_result']['path']).read_text())
        data = training['rows']['matched_dry_acetate_mono_central_'+m]
        check(saved['status']=='complete' and all(saved['checks'].values()), m+' saved admission')
        check(g==data['geometry'] and g['ligand_indices']==list(range(1,8)) and g['metal_index']==0,
              m+' exact center/ligand')
        check(g['symbols']==[m,'C','C','O','O','H','H','H'] and g['charge']==2, m+' identity/state')
        raw = Path(old['output']['path']).read_text()
        lines = [l.strip() for l in Path(old['engrad']['path']).read_text().splitlines() if l.strip() and not l.lstrip().startswith('#')]
        n = int(lines[0]); e = float(lines[1]); gr = [float(v) for v in lines[2:2+3*n]]
        check(n==8 and len(gr)==24 and all(math.isfinite(v) for v in [e]+gr), m+' raw gradient dimension/finite')
        forces = [[-v*HA/BOHR for v in gr[3*i:3*i+3]] for i in range(n)]
        check(forces==saved['forces_kcal_mol_A']==data['reference']['forces'], m+' exact force units/sign/replay')
        final = re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+0-9.Ee]+)',raw)
        check(len(final)==1 and float(final[0])==e==saved['energy_hartree'], m+' final/engrad energy')
        electronic = re.findall(r'^Total Energy\s*:\s*('+FLOAT+r')\s+Eh',raw,re.M)
        disp = re.findall(r'^Dispersion correction\s+('+FLOAT+r')',raw,re.M)
        check(len(electronic)==len(disp)==1 and abs(float(electronic[0])+float(disp[0])-e)<1e-8,
              m+' postintegration electronic plus D4 closure')
        check('DISPERSION GRADIENT' in raw and 'COSX-gradient' in raw and 'ORCA SCF GRADIENT CALCULATION' in raw,
              m+' native gradient components present')
        c = [math.fsum(x[i][k] for i in range(1,8))/7 for k in range(3)]
        radial = sub(x[3],x[0]); u = [t/norm(radial) for t in radial]
        axis = cross(u,sub(x[4],x[3])); check(norm(axis)>1e-12,m+' nondegenerate fixed axis')
        axis = [t/norm(axis) for t in axis]
        deriv = {'radial':-math.fsum(dot(forces[i],u) for i in range(1,8)),
                 'rotation':-math.fsum(dot(axis,cross(sub(x[i],c),forces[i])) for i in range(1,8))}
        for k,v in deriv.items(): check(abs(v-old['analytic'][k])<1e-10,m+' independent '+k+' derivative')
        directions[m] = (g,c,u,axis)
        stage = old['strict_native_replay']
        check(stage['dft_xc_grid_metadata']['PostSCFGGA']=='off' and bool(stage['dft_xc_grid_metadata']['settings']),m+' XC metadata')
        check(stage['native_energy_stages']['has_analytic_gradient_stage'] and stage['native_energy_stages']['recomputed_final_exchange'],m+' COSX stages')
        centers[m] = {'derivatives_kcal_mol_A':deriv,'rotation_axis':axis,'centroid_A':c,
            'net_force_kcal_mol_A':[math.fsum(v[k] for v in forces) for k in range(3)],
            'final_exchange_change_hartree':stage['native_energy_stages']['final_exchange_energy_change_hartree'],
            'PostSCFGGA':stage['dft_xc_grid_metadata']['PostSCFGGA'],
            'SCF_mode':stage['native_scf_mode'],'SCF_stop_signals':stage['native_scf_stop_signals'],
            'source_energy_gradient_surface_equivalence_proven':False}
    for cell in p['cells']:
        m=cell['metal']; g,c,u,axis=directions[m]; q=cell['q_A']; theta=q/1.
        raw=verify(cell['source_geometry']).read_text().splitlines(); x=[[float(z) for z in line.split()[1:4]] for line in raw[2:]]
        check([line.split()[0] for line in raw[2:]]==g['symbols'] and x==cell['positions_A'] and x[0]==g['positions_A'][0],cell['cell_id']+' raw labels/positions')
        expected=[g['positions_A'][0]]
        for origin in g['positions_A'][1:]:
            if cell['direction']=='radial': target=[a+q*b for a,b in zip(origin,u)]
            else:
                v=sub(origin,c); axial=[dot(axis,v)*t for t in axis]; perp=sub(v,axial); tangent=cross(axis,v)
                target=[c[k]+axial[k]+math.cos(theta)*perp[k]+math.sin(theta)*tangent[k] for k in range(3)]
            expected.append(target)
        err=max(abs(a-b) for aa,bb in zip(x,expected) for a,b in zip(aa,bb))
        internal=max(abs(norm(sub(x[i],x[j]))-norm(sub(g['positions_A'][i],g['positions_A'][j]))) for i in range(1,8) for j in range(1,8))
        check(err<=1e-12 and internal<=1e-9,cell['cell_id']+' independent rigid transform')
        deck=verify(cell['template_input']).read_text(); original=Path(p['central'][m]['executed_input']['path']).read_text()
        check(original.count(' EnGrad')==1 and deck==original.replace(' EnGrad',''),cell['cell_id']+' exact sole deck change')
        check(cell['expected_electrons']==40 and cell['expected_basis_functions']==189 and cell['expected_ecp_core']==(46 if m=='La' else 55)
              and cell['charge']==2 and cell['multiplicity']==1 and not any(cell['ghost_mask']),cell['cell_id']+' exact state/basis')
        records.append({'cell_id':cell['cell_id'],'transform_max_A':err,'rigid_distance_max_A':internal})
    check(all(next(c for c in p['cells'] if c['cell_id']==cid.replace('__La__','__Dy__'))['positions_A']==next(c for c in p['cells'] if c['cell_id']==cid)['positions_A'] for cid in ids if '__La__' in cid),'paired literal geometry')
    result={'status':'passed_independent_saved_preparation_checks','implementation':pin(D/'REFERENCE.py'),'reviewer':pin(__file__),'panel':pin(PANEL),
        'checks_passed':len(checks),'checks':checks,'centers':centers,'records':records,'pins_verified':list(seen.values()),
        'new_molecular_calls':0,'author_controls_rerun':False,'held_outcomes_read':False}
    out=W/'independent_review_v1';out.mkdir(parents=True,exist_ok=True)
    target=out/'PREPARATION_CHECKS.json';assert not target.exists()
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'result':pin(target),'checks_passed':len(checks),'centers':centers},indent=2))

if __name__=='__main__': main()
