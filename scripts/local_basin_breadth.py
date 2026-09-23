"""Conditional one-dimensional rotor integrals on pinned real donor geometry.

No molecular engine: finite-candidate scoring is provided by the shared pool.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import shutil
import numpy as np
from scipy.constants import R
from scipy.special import logsumexp
from affordable_common import HA_TO_KCAL, InvalidArtifact, paired, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, write_xyz
from mace_site_kinematics import Kinematics

PROTOCOL = 'physical_single_rotor_conditional_breadth_v1'
SETTINGS = {'temperature_K': 300., 'fine_nodes': 65, 'coarse_stride': 2,
            'inner_domain_fraction': .75, 'maximum_heavy_displacement_A': .8,
            'angular_ceiling_radian': .8, 'measure': 'dq_radian/(2*pi)',
            'per_metal_quadrature_tolerance_kcal_mol': .05,
            'paired_quadrature_tolerance_kcal_mol': .10,
            'per_metal_extent_tolerance_kcal_mol': .05,
            'paired_extent_tolerance_kcal_mol': .10,
            'outer_quarter_probability_maximum': .01,
            'selector': 'first_sidechain_torsion_in_existing_common_four_mode_order'}
CASES = ('1H4I', '4MAE', 'q9z4j7-pqq-la_model',
         'a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4')
RT = R * SETTINGS['temperature_K'] / 4184.


def mode_domain(task, mapping):
    kin = Kinematics(mapping)
    by_id = {m['id']: (i, m) for i, m in enumerate(kin.modes)}
    selected = [(by_id[r['id']]) for r in task['selector']['selected']
                if by_id[r['id']][1]['kind'] == 'sidechain_torsion']
    if not selected: raise InvalidArtifact('no supported selected sidechain torsion')
    k, mode = selected[0]
    if mode['unit'] != 'radian' or kin.metal in mode['moving_indices']:
        raise InvalidArtifact('unsupported rotor or moving substituted metal')
    i, j = mode['axis_indices']; axis = kin.positions[j]-kin.positions[i]
    axis = axis/np.linalg.norm(axis)
    heavy = sorted(set(mode['moving_indices']) & set(kin.heavy))
    if not heavy: raise InvalidArtifact('empty physical heavy rotor')
    radius = np.linalg.norm(np.cross(kin.positions[heavy]-kin.positions[i], axis), axis=1)
    maximum = float(radius.max())
    if maximum <= 1e-10: raise InvalidArtifact('degenerate rotor')
    theta = min(SETTINGS['angular_ceiling_radian'],
                2*math.asin(min(1., SETTINGS['maximum_heavy_displacement_A']/(2*maximum))))
    return k, theta, {'mode': mode, 'moving_heavy_indices': heavy,
                      'maximum_axis_radius_A': maximum,
                      'domain_radian': [-theta, theta]}


def geometry(kin, mode_index, angle, symbols):
    """One-mode exact geometry/mapping checks; no scientific evaluations."""
    q = np.zeros(len(kin.modes)); q[mode_index] = angle
    p, core, j, cj = kin.evaluate(q)
    z = np.zeros(len(q)); z[mode_index] = 1e-6
    plus = kin.evaluate(q+z); minus = kin.evaluate(q-z)
    jac_error = max(float(np.max(np.abs((plus[0]-minus[0])/2e-6-j[mode_index]))),
                    float(np.max(np.abs((plus[1]-minus[1])/2e-6-cj[mode_index]))))
    initial_j = kin.evaluate(np.zeros(len(q)))[2][mode_index]
    # Each atom's squared speed is constant: any physical mass-weighted metric
    # is therefore constant too. Caps are explicitly not counted in this metric.
    metric_error = float(np.max(np.abs(np.sum(j[mode_index]**2, axis=1)-np.sum(initial_j**2, axis=1))))
    moving = set(kin.modes[mode_index]['moving_indices'])
    fixed = [i for i in range(len(p)) if i not in moving]
    fixed_ok = np.allclose(p[fixed], kin.positions[fixed], atol=1e-12, rtol=0)
    bond_error = max(abs(np.linalg.norm(p[a]-p[b])-np.linalg.norm(kin.positions[a]-kin.positions[b]))
                     for a,b in kin.data['bonds'])
    maximum = float(np.linalg.norm(p[kin.heavy]-kin.positions[kin.heavy], axis=1).max())
    limits = np.where(np.array(symbols)[:,None]=='H', .55, 1.)
    limits = np.minimum(limits, limits.T)
    d0 = np.linalg.norm(kin.core[:,None]-kin.core[None,:], axis=2)
    d1 = np.linalg.norm(core[:,None]-core[None,:], axis=2)
    clashes = np.argwhere(np.tril((d1<limits)&(d0>=limits), -1)).tolist()
    ok = (jac_error<=1e-7 and metric_error<=1e-8 and bond_error<=1e-9 and fixed_ok
          and maximum<=SETTINGS['maximum_heavy_displacement_A']+1e-10 and not clashes)
    return q, core, {'pass': bool(ok), 'jacobian_max_error_A_per_radian': jac_error,
                     'atomwise_squared_speed_max_error_A2': metric_error,
                     'physical_metric_constant': metric_error<=1e-8,
                     'bond_max_error_A': float(bond_error), 'fixed_atoms_unchanged': bool(fixed_ok),
                     'maximum_heavy_displacement_A': maximum, 'new_severe_overlaps': clashes}


def prepare(inputs, agreement, output):
    source = read_json(inputs)
    chosen = [c for c in source['cases'] if c.get('basin_four')]
    if tuple(c['case_id'] for c in chosen) != CASES:
        raise InvalidArtifact('frozen four-source membership/order differs')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    cache = {}; cases = []
    for c in chosen:
        for key in ('pool_manifest', 'pool_collection', 'proposal_collection'): verify(c[key])
        cp = c['proposal_collection']['path']
        if cp not in cache:
            collection = read_json(cp); pin = collection['manifest']
            cache[cp] = (pin, read_json(verify(pin)))
        source_pin, manifest = cache[cp]
        tasks = [next(t for t in manifest['tasks'] if t['case_id']==c['case_id'] and t['metal']==z)
                 for z in ('Ca','La')]
        ca, la = tasks
        paired(verify(la['xyz']), verify(ca['xyz']), la['charge'], ca['charge'])
        maps = [read_json(verify(t['mapping']))['context'] for t in tasks]
        if maps[0]!=maps[1] or ca['selector']!=la['selector']:
            raise InvalidArtifact('paired physical mapping/selector differs')
        k, theta, domain = mode_domain(ca, maps[0]); kin = Kinematics(maps[0])
        rows = xyz(verify(ca['xyz'])); symbols = [r[0] for r in rows]
        if np.max(np.abs(kin.core-np.array([r[1:] for r in rows])))>1e-12:
            raise InvalidArtifact('mapping/source origin differs beyond declared 1e-12 A precision')
        nodes = []
        for i, angle in enumerate(np.linspace(-theta, theta, SETTINGS['fine_nodes'])):
            q, coords, check = geometry(kin, k, float(angle), symbols)
            cid = 'origin' if i==32 else f'basin_{i:02d}'
            if i==32:
                pin = ca['xyz']  # exact archived q0, not a rounded rewrite
            else:
                path = out/'coordinates'/c['case_id']/(cid+'.xyz')
                path.parent.mkdir(parents=True, exist_ok=True)
                write_xyz(path, [(r[0], *p) for r,p in zip(rows,coords)]); pin=record(path)
            nodes.append({'id':cid,'index':i,'angle_radian':float(angle),'full_q':q.tolist(),
                          'xyz':pin,'geometry_checks':check,
                          'status':'prepared' if check['pass'] else 'unsupported_geometry'})
        cases.append({'case_id':c['case_id'], 'source':c, 'source_manifest':source_pin,
                      'source_task_ids':{t['metal']:t['task_id'] for t in tasks},
                      'source_tasks':tasks, 'mode_index':k, **domain,
                      'points':nodes, 'status':'prepared' if all(p['geometry_checks']['pass'] for p in nodes)
                      else 'incomplete_physical_domain', 'chemical_states_changed':False})
    impl = out/'implementation'; impl.mkdir()
    pins={}
    for name in ('local_basin_breadth.py','mace_site_kinematics.py'):
        target=impl/name; shutil.copyfile(Path(__file__).parent/name,target); pins[name]=record(target)
    design={'protocol_id':PROTOCOL,'settings':SETTINGS,'inputs':record(inputs),'agreement':record(agreement),
            'cases':cases,'implementation':pins,'new_DFT_calls':0,'new_optimizer_starts':0,
            'maximum_new_MACE_calls':512,'maximum_new_GFN2_calls':1024,'production_changed':False}
    write_new(out/'design.json',design)
    result=validate(out/'design.json');write_new(out/'PREFLIGHT.json',result)
    return result


def validate(design):
    d=read_json(design)
    if d['protocol_id']!=PROTOCOL or d['settings']!=SETTINGS:raise InvalidArtifact('protocol/settings changed')
    verify(d['inputs']);verify(d['agreement'])
    for pin in d['implementation'].values():verify(pin)
    if tuple(c['case_id'] for c in d['cases'])!=CASES:raise InvalidArtifact('case membership differs')
    unsupported=[]
    for c in d['cases']:
        verify(c['source_manifest'])
        ca,la=c['source_tasks'];paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
        maps=[read_json(verify(t['mapping']))['context'] for t in (ca,la)]
        if maps[0]!=maps[1]:raise InvalidArtifact('paired mapping differs')
        k,theta,domain=mode_domain(ca,maps[0]);kin=Kinematics(maps[0]);rows=xyz(verify(ca['xyz']))
        if c['mode_index']!=k or c['domain_radian']!=domain['domain_radian']:raise InvalidArtifact('domain changed')
        grid=np.linspace(-theta,theta,65)
        if len(c['points'])!=65:raise InvalidArtifact('incomplete declared grid')
        for i,p in enumerate(c['points']):
            if p['index']!=i or p['angle_radian']!=grid[i]:raise InvalidArtifact('grid order changed')
            q=np.zeros(len(kin.modes));q[k]=grid[i]
            if p['full_q']!=q.tolist():raise InvalidArtifact('nonselected motion')
            coords=xyz(verify(p['xyz']))
            if [r[0] for r in coords]!=[r[0] for r in rows]:raise InvalidArtifact('chemical composition changed')
            if np.max(np.abs(np.array([r[1:] for r in coords])-kin.evaluate(q)[1]))>1e-12:
                raise InvalidArtifact('coordinate mapping differs')
            if not p['geometry_checks']['pass']:unsupported.append([c['case_id'],p['id']])
    return {'status':'prepared_dry_run_pass','cases':4,'nodes_per_case':65,'declared_metal_cells':520,
            'unsupported_geometry_points':unsupported,'molecular_calls':0,'design':record(design)}


def specification(design, output):
    """Export the root-owned finite-candidate interface without a new runner."""
    validate(design); d=read_json(design)
    cases=[]
    for c in d['cases']:
        candidates=[]
        for p in c['points']:
            # New aliases must not shadow the archived origin identity.
            name='basin_32' if p['id']=='origin' else p['id']
            if p['geometry_checks']['pass']:
                candidates.append({'id':name,'full_q':p['full_q'],'coordinate':p['xyz'],
                                   'point_index':p['index'],'angle_radian':p['angle_radian']})
        cases.append({'case_id':c['case_id'],'status':'prepared' if candidates else 'unavailable',
                      'candidates':candidates,'quadrature_design':record(design),
                      'declared_points':65,'unsupported_points':[p['id'] for p in c['points']
                                                               if not p['geometry_checks']['pass']]})
    value={'branch':'local_basin_breadth','inputs':d['inputs'],'agreement':d['agreement'],
           'coordinate_limits':{'maximum_angle_radian':.8,'maximum_heavy_displacement_A':.8},
           'maximum_candidates_per_case':65,'cases':cases,'quadrature_design':record(design)}
    write_new(output,value);return {'specification':record(output),'new_molecular_calls':0}


def simpson_weights(q):
    q=np.asarray(q,dtype=float)
    if q.ndim!=1 or len(q)<3 or len(q)%2!=1 or not np.isfinite(q).all():
        raise InvalidArtifact('odd finite quadrature grid required')
    spacing=np.diff(q)
    if spacing[0]<=0 or not np.allclose(spacing,spacing[0],rtol=1e-12,atol=1e-14):
        raise InvalidArtifact('uniform increasing grid required')
    w=np.ones(len(q));w[1:-1:2]=4;w[2:-1:2]=2
    return w*spacing[0]/3.


def integral(q, relative_energy):
    """Pure algebra; input electronic energies must come from actual receipts."""
    q=np.asarray(q);e=np.asarray(relative_energy,dtype=float)
    if e.shape!=q.shape or not np.isfinite(e).all():raise InvalidArtifact('missing/nonfinite grid energy')
    w=simpson_weights(q);minimum=float(np.min(e))
    log_terms=np.log(w)-(e-minimum)/RT
    log_width=float(logsumexp(log_terms));free=minimum-RT*(log_width-math.log(2*math.pi))
    return {'F_relative_kcal_mol':free,'minimum_relative_kcal_mol':minimum,
            'minimum_angle_radian':float(q[np.argmin(e)]),
            'width_radian':float(math.exp(log_width)), 'log_width_radian':log_width,
            'simpson_node_probability':np.exp(log_terms-log_width).tolist()}


def paired_curves(q, ca, la):
    q=np.asarray(q,dtype=float)
    if len(q)!=65 or not np.allclose(q,-q[::-1],atol=1e-14,rtol=0):
        raise InvalidArtifact('frozen symmetric 65-node grid required')
    result={};failures=[]
    for metal,values in (('Ca',ca),('La',la)):
        e=np.asarray(values,dtype=float)
        if e.shape!=(65,) or not np.isfinite(e).all():raise InvalidArtifact('complete real curve required')
        e=e-e[32]
        fine=integral(q,e);coarse=integral(q[::2],e[::2]);inner=integral(q[8:57],e[8:57])
        inner_coarse=integral(q[8:57:2],e[8:57:2])
        df=fine['F_relative_kcal_mol']-coarse['F_relative_kcal_mol']
        de=fine['F_relative_kcal_mol']-inner['F_relative_kcal_mol']
        # This uses the identical dq measure and a complete quadrature on each
        # subdomain, not partial Simpson endpoint weights or width normalization.
        tail=float(-math.expm1(min(0.,de/RT)))
        checks={'positive_domain_measure':de<=1e-10,
                'quadrature':abs(df)<=SETTINGS['per_metal_quadrature_tolerance_kcal_mol'],
                'extent':abs(de)<=SETTINGS['per_metal_extent_tolerance_kcal_mol'],
                'outer_quarter_weight':tail<=SETTINGS['outer_quarter_probability_maximum'],
                'minimum_in_inner_domain':bool(abs(fine['minimum_angle_radian'])<abs(q[56]))}
        failures += [metal+':'+k for k,v in checks.items() if not v]
        result[metal]={'fine':fine,'coarse':coarse,'inner':inner,'inner_coarse':inner_coarse,'checks':checks,
                       'fine_minus_coarse_kcal_mol':df,'outer_minus_inner_kcal_mol':de,
                       'outer_quarter_weight':tail,'relative_energy_kcal_mol':e.tolist()}
    def paired(name,key):return result['Ca'][name][key]-result['La'][name][key]
    total=paired('fine','F_relative_kcal_mol');minima=paired('fine','minimum_relative_kcal_mol')
    width=-RT*(result['Ca']['fine']['log_width_radian']-result['La']['fine']['log_width_radian'])
    qerr=total-paired('coarse','F_relative_kcal_mol');derr=total-paired('inner','F_relative_kcal_mol')
    if abs(qerr)>SETTINGS['paired_quadrature_tolerance_kcal_mol']:failures.append('paired:quadrature')
    if abs(derr)>SETTINGS['paired_extent_tolerance_kcal_mol']:failures.append('paired:extent')
    if not math.isclose(total,minima+width,abs_tol=1e-10,rel_tol=0):raise InvalidArtifact('component algebra mismatch')
    result.update(status='conditional_integral_qualified' if not failures else 'conditional_correction_unavailable',
                  failed_checks=failures,raw_deltaR_kcal_mol=total,raw_minima_deltaR_kcal_mol=minima,
                  raw_breadth_deltaR_kcal_mol=width,paired_quadrature_error_kcal_mol=qerr,
                  paired_extent_effect_kcal_mol=derr,qualified_deltaR_kcal_mol=None if failures else total,
                  qualified_breadth_deltaR_kcal_mol=None if failures else width,
                  chemical_state_populations=None,whole_pocket_entropy=None)
    return result


def component_energy(components):
    return (components['MACE_eV']*EV_TO_KCAL+
            (components['GFN2_ALPB_hartree']-components['GFN2_vacuum_hartree'])*HA_TO_KCAL)


def roughness(q, values):
    """Descriptive actual-grid differences, not a fit or new acceptance gate."""
    q=np.asarray(q,dtype=float); e=np.asarray(values,dtype=float)
    first=np.diff(e);second=np.diff(e,n=2)
    i=int(np.argmax(np.abs(first)));j=int(np.argmax(np.abs(second)))+1
    return {'maximum_adjacent_change_kcal_mol':float(np.max(np.abs(first))),
            'maximum_adjacent_interval_radian':[float(q[i]),float(q[i+1])],
            'maximum_second_difference_kcal_mol':float(np.max(np.abs(second))),
            'maximum_second_difference_center_radian':float(q[j]),
            'interpretation':'descriptive; curvature or SCF branch changes are not distinguished',
            'smoothing_or_fitting':False}


def analyze(design, collection, output):
    d=read_json(design);validate(design);coll=read_json(collection)
    cm=read_json(verify(coll['manifest']));spec=read_json(verify(cm['specification']))
    if spec.get('quadrature_design')!=record(design):raise InvalidArtifact('collection belongs to another grid design')
    results=[]
    for c in d['cases']:
        row=next((r for r in coll['cases'] if r['case_id']==c['case_id']),None)
        missing=[];curves={};cells={}
        for z in ('Ca','La'):
            cells[z]=[];curves[z]=[]
            for p in c['points']:
                representative=(row or {}).get('aliases',{}).get(p['id'],{}).get('representative',p['id'])
                cell=(row or {}).get('matrix',{}).get(z,{}).get(representative,{})
                components=cell.get('components')
                if components is None:
                    native=read_json(verify(cell['MACE'])) if cell.get('MACE') else {}
                    components={'MACE_eV':native.get('energy_eV') if native.get('status')=='complete' else None,
                                'GFN2_ALPB_hartree':cell.get('low',{}).get('alpb',{}).get('energy_hartree'),
                                'GFN2_vacuum_hartree':cell.get('low',{}).get('vacuum',{}).get('energy_hartree')}
                cells[z].append({'point_id':p['id'],'representative':representative,'components':components,
                                 'status':cell.get('status','missing'),
                                 'unavailable_reasons':{m:v.get('reason') for m,v in cell.get('low',{}).items()
                                                        if 'energy_hartree' not in v}})
                if cell.get('status')!='complete' or not p['geometry_checks']['pass']:
                    missing.append({'metal':z,'point':p['id'],'status':cell.get('status','missing'),
                                    'reasons':cells[z][-1]['unavailable_reasons']});continue
                curves[z].append(component_energy(components))
        if missing:
            result={'status':'conditional_correction_unavailable','missing':missing,
                    'qualified_deltaR_kcal_mol':None,'qualified_breadth_deltaR_kcal_mol':None}
        else:
            result=paired_curves([p['angle_radian'] for p in c['points']],curves['Ca'],curves['La'])
            result['origin_R_kcal_mol']=curves['Ca'][32]-curves['La'][32]
            result['raw_conditional_R_kcal_mol']=result['origin_R_kcal_mol']+result['raw_deltaR_kcal_mol']
            result['raw_sampled_minima_R_kcal_mol']=result['origin_R_kcal_mol']+result['raw_minima_deltaR_kcal_mol']
            angles=[p['angle_radian'] for p in c['points']]
            for z in ('Ca','La'):
                comps=[v['components'] for v in cells[z]]
                native=[v['MACE_eV']*EV_TO_KCAL for v in comps]
                solvent=[(v['GFN2_ALPB_hartree']-v['GFN2_vacuum_hartree'])*HA_TO_KCAL for v in comps]
                result[z]['grid_roughness']={name:roughness(angles,values) for name,values in
                    (('native_MACE',native),('ALPB_minus_vacuum',solvent),('composite',curves[z]))}
        results.append({'case_id':c['case_id'],'mode':c['mode'],'domain_radian':c['domain_radian'],
                        'cells':cells,'result':result})
    value={'protocol_id':PROTOCOL,'design':record(design),'collection':record(collection),'settings':SETTINGS,
           'cases':results,'analysis_implementation':record(__file__),
           'production_changed':False,'new_molecular_calls':0}
    write_new(output,value);return {'output':record(output),'cases':[{ 'case_id':r['case_id'],
        'status':r['result']['status'],'qualified_deltaR_kcal_mol':r['result']['qualified_deltaR_kcal_mol']} for r in results]}


def plot(result, output):
    """Editable reporting artifacts from actual collected grid components only."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.size':9})
    data=read_json(result);design=read_json(verify(data['design']))
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    titles={'1H4I':'1H4I','4MAE':'4MAE','q9z4j7-pqq-la_model':'Q9Z4J7',
            CASES[3]:'A0ACD6B9F2 · Ca-conditioned sample 4'}
    handles=[]
    for ax,case in zip(axes.flat,data['cases']):
        source=next(c for c in design['cases'] if c['case_id']==case['case_id'])
        q=np.array([p['angle_radian'] for p in source['points']]);theta=q[-1]
        ax.axvspan(-.75*theta,.75*theta,color='#edf1f4',zorder=0)
        ax.axhline(0,color='#cccccc',lw=.7);ax.axvline(0,color='#cccccc',lw=.7)
        for metal,style in (('Ca','-'),('La','--')):
            values={c['point_id']:c['components'] for c in case['cells'][metal]}
            origin=values.get('origin')
            if origin is None:continue
            native=[];solvent=[]
            for point in source['points']:
                comp=values.get(point['id'])
                native.append(np.nan if comp is None or comp.get('MACE_eV') is None else (comp['MACE_eV']-origin['MACE_eV'])*EV_TO_KCAL)
                solvent.append(np.nan if comp is None or any(comp.get(k) is None for k in ('GFN2_ALPB_hartree','GFN2_vacuum_hartree')) else ((comp['GFN2_ALPB_hartree']-comp['GFN2_vacuum_hartree'])-
                                  (origin['GFN2_ALPB_hartree']-origin['GFN2_vacuum_hartree']))*HA_TO_KCAL)
            native=np.array(native);solvent=np.array(solvent)
            for name,e,color,width in (('Native MACE',native,'#888888',1.),('ALPB − vacuum',solvent,'#337ab7',1.3),
                                        ('Composite',native+solvent,'#b44031',1.8)):
                line,=ax.plot(q,e,color=color,ls=style,lw=width,label=metal+' · '+name)
                if ax is axes.flat[0]:handles.append(line)
            composite=native+solvent
            if np.isfinite(composite).any():
                i=int(np.nanargmin(composite));ax.scatter(q[i],composite[i],s=28,marker='o' if metal=='Ca' else 's',
                                                         color='#b44031',zorder=5)
        r=case['result'];flags=r.get('failed_checks',[])
        subtitle=('quadrature/domain checks pass' if r['status']=='conditional_integral_qualified'
                  else 'correction unavailable')
        ax.set_title(titles[case['case_id']]+'\n'+source['mode']['id']+' · '+subtitle)
        ax.set_xlabel('Donor torsion displacement (radian)');ax.set_ylabel('Energy relative to own q0 (kcal/mol)')
        ax.set_xlim(-theta,theta)
        if flags:
            annotation='Failed quadrature and domain checks\nLa minimum outside inner domain'
        elif r.get('missing'):
            annotation=str(len(r['missing']))+' missing SCF points\nComplete integral unavailable'
        else:
            annotation='Conditional numerical result\nElectronic surface not independently qualified'
        ax.text(.02,.98,annotation,transform=ax.transAxes,ha='left',va='top',fontsize=7,
                bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
    ordered=[handles[i] for i in (0,3,1,4,2,5)]
    fig.legend(handles=ordered,loc='outside lower center',ncol=3,frameon=False)
    fig.suptitle('Direct donor-torsion basin profiles · fixed spectators and chemical states · 300 K',fontsize=12)
    for ext in ('svg','pdf'):fig.savefig(out/('basin_profiles.'+ext))
    plt.close(fig)
    caption=("Native MACE, native GFN2 ALPB-minus-vacuum, and their composite are each "
             "shifted by that metal's value at the identical source q0. Solid: Ca; dashed: La. "
             "Circles/squares mark minima among computed composite points (incomplete grids remain unqualified). Grey shading is the fixed inner "
             "75% domain; full axes show the complete 65-node integration domain. Curves join "
             "computed grid points without fitting. Missing nodes break curves and prevent "
             "qualification. These are conditional single-coordinate electronic/configurational "
             "diagnostics on consumed sources, not chemical-state populations or DFT-validated "
             "adaptive candidates. Exact numeric gates and limitations are in the frozen plan.\n")
    (out/'CAPTION.md').write_text(caption)
    write_new(out/'PROVENANCE.json',{'result':record(result),'plot_implementation':record(__file__),
              'artifacts':[record(out/n) for n in ('basin_profiles.svg','basin_profiles.pdf','CAPTION.md')],
              'new_molecular_calls':0,'curve_fitting':False})
    return {'directory':str(out),'provenance':record(out/'PROVENANCE.json')}


def main():
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='op',required=True)
    p=sub.add_parser('prepare');p.add_argument('--inputs',required=True);p.add_argument('--agreement',required=True);p.add_argument('--output',required=True)
    p=sub.add_parser('dry-run');p.add_argument('--design',required=True)
    p=sub.add_parser('specification');p.add_argument('--design',required=True);p.add_argument('--output',required=True)
    p=sub.add_parser('plot');p.add_argument('--result',required=True);p.add_argument('--output',required=True)
    p=sub.add_parser('analyze');p.add_argument('--design',required=True);p.add_argument('--collection',required=True);p.add_argument('--output',required=True)
    a=ap.parse_args()
    if a.op=='prepare':r=prepare(a.inputs,a.agreement,a.output)
    elif a.op=='dry-run':r=validate(a.design)
    elif a.op=='specification':r=specification(a.design,a.output)
    elif a.op=='plot':r=plot(a.result,a.output)
    else:r=analyze(a.design,a.collection,a.output)
    print(json.dumps(r,indent=2,allow_nan=False))

if __name__=='__main__':main()
