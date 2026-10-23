"""Saved scalar/rigid-response objective only. No engine imports or calls."""
import math

FAMILIES = ('acetate_mono', 'acetate_bi', 'acetamide')
WATER = ('water_central', 'aquo8_central')
DONOR = tuple('hydrated_'+f+'_'+s for f in FAMILIES for s in ('minus','central','plus')) + tuple('matched_dry_'+f+'_central' for f in FAMILIES)

def dot(a,b): return sum(x*y for x,y in zip(a,b))
def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def mean(a): return sum(a)/len(a)

def components(bonds):
    unseen=set(range(len(bonds))); out=[]
    while unseen:
        todo=[min(unseen)]; group=[]
        while todo:
            i=todo.pop()
            if i not in unseen: continue
            unseen.remove(i); group.append(i); todo.extend(bonds[i])
        out.append(sorted(group))
    return out

def bodies(row):
    g=row['geometry']; out=[]
    for ids in components(row['input']['bonds_zero_based']):
        if ids==[g['metal_index']]: category='ion'
        elif sorted(g['symbols'][i] for i in ids)==['H','H','O']: category='water'
        else: category='ligand'
        out.append({'indices':ids,'category':category})
    assert sum(x['category']=='ion' for x in out)==1
    assert sum(x['category']=='ligand' for x in out)<=1
    assert sorted(i for b in out for i in b['indices'])==list(range(g['atoms']))
    return out

def response(row, forces):
    x=row['geometry']['positions_A']; result=[]
    for body in row['bodies']:
        ids=body['indices']; center=[mean([x[i][k] for i in ids]) for k in range(3)]
        net=[sum(forces[i][k] for i in ids) for k in range(3)]
        torque=[0.,0.,0.]
        for i in ids:
            a=cross([x[i][k]-center[k] for k in range(3)], forces[i])
            torque=[torque[k]+a[k] for k in range(3)]
        # Torque/1 Angstrom: numerical value is unchanged, units become force.
        result.append({'category':body['category'],'net_force':net,'torque_over_1A':torque,
                       'indices':ids,'center_A':center})
    return result

def rigid_residual(row, a, b):
    aa=response(row,a); bb=response(row,b); cats=sorted(set(v['category'] for v in aa)); out=[]; raw=[]
    for cat in cats:
        pairs=[(u,v) for u,v in zip(aa,bb) if u['category']==cat]
        for u,v in pairs:
            keys=['net_force'] if cat=='ion' else ['net_force','torque_over_1A']
            weight=1/(len(cats)*len(pairs)*len(keys)*3)
            for key in keys:
                for x,y in zip(u[key],v[key]):
                    raw.append(x-y); out.append((x-y)/2*math.sqrt(weight))
    return out,raw

def subtract(a,b): return [[x-y for x,y in zip(u,v)] for u,v in zip(a,b)]
def baseline(data): return {k:{'energy':r['baseline']['energy'],'forces':r['baseline']['forces']} for k,r in data['rows'].items()}
def references(data): return {k:{'energy':r['reference']['energy'],'forces':r['reference']['forces']} for k,r in data['rows'].items()}
def delta(values,q): return values[q+'_Dy']['energy']-values[q+'_La']['energy']
def contrast(values,q): return delta(values,q)-delta(values,'aquo8_central')

def evaluate(data,model,stage):
    assert stage in ('water','donor'); ref=references(data); bases=WATER if stage=='water' else DONOR
    energies={}; vectors={}; raw_force={'Dy':[],'paired':[],'La':[]}; control_vectors=[]
    energies['exchange']=[contrast(model,q)-contrast(ref,q) for q in bases if q!='aquo8_central']
    if stage=='donor':
        energies['addition']=[(delta(model,'hydrated_'+f+'_central')-delta(model,'matched_dry_'+f+'_central'))-
                              (delta(ref,'hydrated_'+f+'_central')-delta(ref,'matched_dry_'+f+'_central')) for f in FAMILIES]
        for metal in ('Dy','La'):
            energies[metal+'_radial']=[(model['hydrated_'+f+'_'+s+'_'+metal]['energy']-model['hydrated_'+f+'_central_'+metal]['energy'])-
                                      (ref['hydrated_'+f+'_'+s+'_'+metal]['energy']-ref['hydrated_'+f+'_central_'+metal]['energy'])
                                      for f in FAMILIES for s in ('minus','plus')]
    for kind in ('Dy','paired','La'):
        vec=[]
        for q in bases:
            row=data['rows'][q+'_Dy']
            if kind=='paired':
                a=subtract(model[q+'_Dy']['forces'],model[q+'_La']['forces'])
                b=subtract(ref[q+'_Dy']['forces'],ref[q+'_La']['forces'])
            else: a=model[q+'_'+kind]['forces']; b=ref[q+'_'+kind]['forces']
            r,raw=rigid_residual(row,a,b);vec.extend(v/math.sqrt(len(bases)) for v in r);raw_force[kind].extend(raw)
        if kind=='La': control_vectors=vec
        else: vectors[kind]=vec
    weights={'exchange':.50,'Dy':.25,'paired':.25} if stage=='water' else {'exchange':.25,'addition':.15,'Dy_radial':.10,'Dy':.25,'paired':.25}
    weighted=[]; blocks={}
    for name,w in weights.items():
        v=vectors[name] if name in vectors else [z/math.sqrt(len(energies[name])) for z in energies[name]]
        blocks[name]=sum(z*z for z in v);weighted.extend(math.sqrt(w)*z for z in v)
    controls={}
    for kind in ('La','Dy','paired'):
        v=control_vectors if kind=='La' else vectors[kind]
        rms=2*math.sqrt(sum(z*z for z in v)); mx=max(abs(z) for z in raw_force[kind])
        controls[kind]={'category_balanced_rigid_RMS_kcal_mol_A':rms,'max_body_component_error_kcal_mol_A':mx,
                        'force_passed':rms<=2 and mx<=5}
    controls['La']['radial_work_passed']=all(abs(v)<=1 for v in energies.get('La_radial',[]))
    controls['La']['passed']=controls['La']['force_passed'] and controls['La']['radial_work_passed']
    epass=all(abs(v)<=1 for k,a in energies.items() if k!='La_radial' for v in a)
    return {'objective':sum(v*v for v in weighted),'weighted_residuals':weighted,'blocks':blocks,
            'energy_errors_kcal_mol':energies,'controls':controls,'energy_passed':epass,
            'candidate_development_passed':epass and all(controls[k]['force_passed'] for k in ('Dy','paired')),
            'La_control_passed':controls['La']['passed'],'physical_accuracy_qualified':False}

def radial_screen(data):
    rows=data['rows']; out=[]
    for family in FAMILIES:
        metalrows=[]
        for metal in ('La','Dy'):
            rr=[rows['hydrated_'+family+'_'+s+'_'+metal] for s in ('minus','central','plus')]
            g=rr[1]['geometry']; ids=g['ligand_indices']; u=g['radial_direction']
            assert abs(dot(u,u)-1)<1e-10
            for r,h in zip(rr,(-.15,0,.15)):
                for i,(x,y) in enumerate(zip(r['geometry']['positions_A'],g['positions_A'])):
                    assert max(abs(x[k]-y[k]-(h*u[k] if i in ids else 0)) for k in range(3))<2e-12
            e=[r['reference']['energy'] for r in rr]
            slopes=[-sum(dot(r['reference']['forces'][i],u) for i in ids) for r in rr]
            fd=(e[2]-e[0])/.30; tol=max(1.,.1*max(abs(fd),abs(slopes[1])))
            v={'family':family,'metal':metal,'step_A':.15,'moving_atoms':ids,'direction':u,
               'endpoint_and_central_analytic_slopes_kcal_mol_A':slopes,'secant_kcal_mol_A':fd,
               'central_minus_secant_kcal_mol_A':slopes[1]-fd,'screen_tolerance_kcal_mol_A':tol,
               'screen_passed':abs(slopes[1]-fd)<=tol,'small_step_qualified':False}
            out.append(v);metalrows.append(v)
        a,b=metalrows; an=b['endpoint_and_central_analytic_slopes_kcal_mol_A'][1]-a['endpoint_and_central_analytic_slopes_kcal_mol_A'][1]
        fd=b['secant_kcal_mol_A']-a['secant_kcal_mol_A'];tol=max(1.,.1*max(abs(an),abs(fd)))
        out.append({'family':family,'metal':'Dy_minus_La','step_A':.15,'analytic_kcal_mol_A':an,
                    'secant_kcal_mol_A':fd,'central_minus_secant_kcal_mol_A':an-fd,
                    'screen_tolerance_kcal_mol_A':tol,'screen_passed':abs(an-fd)<=tol,'small_step_qualified':False})
    return out
