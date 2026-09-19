"""Source-defined water networks and constrained native-DFT orientation trials.

Research paths only. Occupancy probabilities require free-energy terms which
are not supplied by these electronic calculations.
"""
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path
import re
import shutil

import numpy as np

import affordable_peptide as peptide
import hydration_square as square
from affordable_common import (HA_TO_KCAL, InvalidArtifact, cache_key, paired,
                               read_json, record, verify, write_new, xyz)

PROTOCOL = 'native_r2scan3c_cpcm_water_network_orientations_v1'
METHOD = '! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 Opt'


def residue_id(meta):
    return tuple(meta[k] for k in ('chain', 'resnum', 'insertion_code', 'resname'))


def source_id(meta):
    return residue_id(meta) + (meta['atom'],)


def discover(config):
    """One physical neighbor shell; chemistry is completed through SourceGraph."""
    cfg = read_json(config)
    source = read_json(verify(cfg['source_config']))
    parents = []
    union = set()
    sequences = []
    for item in source['parents']:
        parent, atoms, waters = square.checked_parent(item)
        graph = peptide.SourceGraph(verify(parent['source_structure']),
                                    verify(parent['topology_definition']))
        site_chains = {w['source']['chain'] for w in waters}
        sequences.append({(r.chain, r.resnum, r.insertion_code): r.canonical_resname
                          for r in graph.residues.values()
                          if r.key in graph.supported_topology_residues and r.chain in site_chains})
        selected = {source_id(a['source']) for a in parent['atom_graph']['source_to_qm']
                    if a['kind'] == 'source'}
        contacts, fragments, outer = [], set(), set()
        variable_ids = {residue_id(w['source']) for w in waters}
        for water in waters:
            origin = np.array(atoms[water['oxygen_index']][1:])
            for key, atom in graph.atoms.items():
                if atom.element.name not in ('N', 'O', 'S'):
                    continue
                meta = graph.meta[key]
                if residue_id(meta) == residue_id(water['source']):
                    continue
                distance = float(np.linalg.norm(np.array(tuple(atom.pos)) - origin))
                if distance > cfg['contact_cutoff_A']:
                    continue
                contacts.append({'water': water['source'], 'neighbor': meta,
                                 'distance_A': distance,
                                 'in_original_core': source_id(meta) in selected})
                if meta['canonical_resname'] == 'HOH':
                    if residue_id(meta) not in variable_ids:
                        outer.add(residue_id(meta))
                elif key[:2] not in graph.supported_topology_residues:
                    raise InvalidArtifact(f'unsupported nearby polar atom: {meta}')
                elif meta['atom'] in ('C', 'O', 'OXT'):
                    fragments.add(('peptide', residue_id(meta)))
                elif meta['atom'] == 'N':
                    previous = graph.amide_prev.get(key[:2])
                    if previous is None:
                        raise InvalidArtifact(f'unsupported terminal amine contact: {meta}')
                    fragments.add(('peptide', residue_id(graph.meta[previous])))
                else:
                    fragments.add(('sidechain', residue_id(meta)))
        union.update(fragments)
        parents.append({'item': item, 'parent': parent, 'graph': graph,
                        'contacts': contacts, 'outer': outer,
                        'variable': variable_ids, 'waters': waters})
    for sequence in sequences[1:]:
        common = set(sequence) & set(sequences[0])
        if not common or any(sequence[k] != sequences[0][k] for k in common):
            raise InvalidArtifact('replicate union requires matching indexed protein sequences')
    for data, sequence in zip(parents, sequences):
        data['sequence_audit'] = {'site_chain_modeled_residues': len(sequence),
            'common_residue_count': len(set.intersection(*(set(s) for s in sequences))),
            'common_sequence_identity': True,
            'scope': 'site chains; remote unresolved termini and other deposited chains may differ'}
    return cfg, source, parents, sorted(union)


def selector(rid):
    return dict(zip(('chain', 'resnum', 'insertion_code', 'resname'), rid))


def materialize(data, union, reference):
    parent, graph = data['parent'], data['graph']
    selected = set()
    for a in parent['atom_graph']['source_to_qm']:
        if a['kind'] == 'source' and a['source']['element'] not in ('H', 'D') and a['source']['canonical_resname'] != 'HOH':
            s = a['source']
            selected.add((s['chain_index'], s['residue_index'], s['atom']))
    ledger = [dict(a) for a in parent['charge_ledger'] if a['kind'] != 'water']
    charged = {residue_id(a['source']) for a in ledger if a['kind'] == 'sidechain'}
    for kind, rid in union:
        r = graph.locate(selector(rid))
        if kind == 'peptide':
            selected.update(graph.peptide(r.key))
        else:
            _, checked = peptide.old.sidechain_fragment(r, [], 0)
            selected.update(graph.key(r.key, n) for n in peptide.old.SIDECHAIN_QM_HEAVY_ATOMS[r.canonical_resname])
            if rid not in charged:
                ledger.append({'kind': 'sidechain', 'source': r.source_dict(),
                               'formal_charge': checked['formal_charge']})
                charged.add(rid)
    atoms, mapping = graph.materialize(selected)
    atoms = [('La', *parent['selected_site']['xyz_A'])] + atoms
    water_inventory = []
    for rid in sorted(data['variable'] | data['outer']):
        r = graph.locate(selector(rid))
        rows, checked = peptide.old.water_fragment(r, [], 0)
        if checked['formal_charge'] != 0:
            raise InvalidArtifact('non-neutral source water')
        water_inventory.append(r.source_dict())
        ledger.append({'kind': 'water', 'source': r.source_dict(), 'formal_charge': 0})
        for atom in rows:
            matches = [k for k in graph.by_residue[r.key].values()
                       if peptide.old.atom_to_qm(graph.atoms[k]) == atom]
            if len(matches) != 1:
                raise InvalidArtifact('ambiguous source water atom')
            mapping['source_to_qm'].append({'qm_index': len(atoms), 'kind': 'source',
                                            'source': graph.meta[matches[0]], 'xyz_A': list(atom[1:])})
            atoms.append(atom)
    pseudo = {'atom_graph': mapping, 'charge_ledger': ledger,
              'explicit_water_inventory': water_inventory}
    groups = square.water_groups(pseudo, atoms)
    normalized = square.reference_geometry(atoms, groups, reference)
    # These are prepared coordinates; retain the unmodified source coordinates too.
    for a in mapping['source_to_qm']:
        a['source_xyz_A'] = a['xyz_A']
        a['xyz_A'] = list(normalized[a['qm_index']][1:])
    for w in groups:
        w['role'] = 'variable' if residue_id(w['source']) in data['variable'] else 'frozen_outer'
    return normalized, mapping, ledger, groups


def radial_seed(atoms, groups, mapping):
    """Label-free second orientation: bisector away from metal, source-defined plane."""
    result = list(atoms)
    center = np.array(atoms[0][1:])
    oxygens = [a['qm_index'] for a in mapping['source_to_qm']
               if a['kind'] == 'source' and a['source']['element'] == 'O'
               and a['source']['canonical_resname'] != 'HOH']
    for water in groups:
        if water['role'] != 'variable':
            continue
        o = np.array(atoms[water['oxygen_index']][1:])
        bisector = o-center; bisector /= np.linalg.norm(bisector)
        axis = None
        for i in sorted(oxygens, key=lambda i: (np.linalg.norm(np.array(atoms[i][1:])-o), i)):
            v = np.array(atoms[i][1:])-o
            v -= (v @ bisector)*bisector
            if np.linalg.norm(v) > 1e-6:
                axis = v / np.linalg.norm(v)
                break
        if axis is None:
            raise InvalidArtifact('no source-defined water orientation plane')
        old = np.array([atoms[i][1:] for i in water['hydrogen_indices']])-o
        axial = np.linalg.norm(old.sum(axis=0))/2
        radial = np.linalg.norm(old[0]-old[1])/2
        for i, sign in zip(water['hydrogen_indices'], (1, -1)):
            result[i] = ('H', *(o + axial*bisector + sign*radial*axis))
    return result


def write_xyz(path, atoms, comment):
    Path(path).write_text(str(len(atoms))+'\n'+comment+'\n'+
                         '\n'.join(f'{a[0]} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}' for a in atoms)+'\n')


def constrained_input(atoms, groups, charge, maxiter):
    mobile = sorted(i for w in groups if w['role'] == 'variable' for i in w['hydrogen_indices'])
    if not mobile:
        return '! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF\n* xyzfile '+str(charge)+' 1 core.xyz\n'
    lines = [METHOD, '%geom', f' MaxIter {maxiter}', ' Constraints']
    lines.extend(f'  {{ C {i} C }}' for i in range(len(atoms)) if i not in mobile)
    for w in groups:
        if w['role'] == 'variable':
            o = w['oxygen_index']; a, b = w['hydrogen_indices']
            lines.extend((f'  {{ B {o} {a} C }}', f'  {{ B {o} {b} C }}', f'  {{ A {a} {o} {b} C }}'))
    lines.extend((' end', 'end', f'* xyzfile {charge} 1 core.xyz'))
    return '\n'.join(lines)+'\n'


def prepare(config, output, occupancy=False, case=None, exclude_full=False):
    cfg, source, parents, union = discover(config)
    if exclude_full and not occupancy:
        raise InvalidArtifact('exclude-full applies only to an occupancy continuation')
    verify(cfg['agreement'])
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    reference = xyz(verify(source['water_reference_geometry']))
    tasks, prepared = [], []
    case_filter = case
    if case is not None and case not in [p['item']['case'] for p in parents]:
        raise InvalidArtifact('unknown case filter')
    for data in parents:
        if case_filter is not None and data['item']['case'] != case_filter:
            continue
        atoms, mapping, ledger, groups = materialize(data, union, reference)
        case = data['item']['case']; case_dir = out/case; case_dir.mkdir()
        protein_charge = sum(r['formal_charge'] for r in ledger)
        for metal in ('La', 'Ca'):
            q = protein_charge + (3 if metal == 'La' else 2)
            peptide.old.require_closed_shell(metal, [(metal, *atoms[0][1:])]+atoms[1:], q)
        variable = [w for w in groups if w['role'] == 'variable']
        subsets = [tuple(range(len(variable)))]
        if occupancy:
            subsets = [s for n in range(len(variable)+1) for s in combinations(range(len(variable)), n)]
        if exclude_full:
            subsets = [s for s in subsets if len(s) != len(variable)]
        states = []
        for subset in subsets:
            omitted = {i for j,w in enumerate(variable) if j not in subset for i in w['indices']}
            keep = [i for i in range(len(atoms)) if i not in omitted]
            new_index = {old: new for new,old in enumerate(keep)}
            state_atoms = [atoms[i] for i in keep]
            state_groups = []
            for w in groups:
                if w['oxygen_index'] not in omitted:
                    state_groups.append(dict(w, indices=[new_index[i] for i in w['indices']],
                        oxygen_index=new_index[w['oxygen_index']], hydrogen_indices=[new_index[i] for i in w['hydrogen_indices']]))
            state_mapping = dict(mapping, source_to_qm=[dict(a, qm_index=new_index[a['qm_index']])
                                for a in mapping['source_to_qm'] if a['qm_index'] not in omitted])
            pattern = ''.join('1' if j in subset else '0' for j in range(len(variable)))
            seeds = cfg['orientation_seeds'] if subset else ['source']
            ids = []
            for seed in seeds:
                rows = state_atoms if seed == 'source' else radial_seed(state_atoms, state_groups, state_mapping)
                pair_paths = {}
                for metal in ('La', 'Ca'):
                    tid = f'{case}__{pattern}__{seed}__{metal}'; ids.append(tid)
                    md = out/'tasks'/tid; md.mkdir(parents=True)
                    rows_m = [(metal,*rows[0][1:])]+rows[1:]
                    q = protein_charge+(3 if metal == 'La' else 2)
                    xp = md/'core.xyz'; write_xyz(xp, rows_m, tid)
                    ip = md/'endpoint.inp'; ip.write_text(constrained_input(rows_m,state_groups,q,cfg['optimization_maxiter']))
                    task = {'task_id':tid,'case':case,'pattern':pattern,'seed':seed,'metal':metal,
                            'charge':q,'multiplicity':1,'input':record(ip),'xyz':record(xp),
                            'output_path':str(md/'endpoint.out'),'protocol_id':PROTOCOL,
                            'task_type':'constrained_optimization' if subset else 'single_point',
                            'variable_water_count':len(subset),'groups':state_groups,
                            'mobile_indices':sorted(i for w in state_groups if w['role']=='variable' for i in w['hydrogen_indices']),
                            'all_electron_count_for_parity':peptide.old.require_closed_shell(metal,rows_m,q)}
                    # Opt prints analytic gradients at its optimization iterates;
                    # its final energy-only evaluation need not write .engrad.
                    task['cache_key'] = cache_key(task | {'configuration': record(config)})
                    tasks.append(task); pair_paths[metal] = xp
                paired(pair_paths['La'],pair_paths['Ca'],protein_charge+3,protein_charge+2)
            states.append({'pattern':pattern,'kept_variable_ids':[variable[j]['source'] for j in subset],
                           'tasks':ids,'source_mapping':state_mapping,'retained_full_indices':keep})
        full = {'case':case,'parent':data['item']['preparation'],'source_structure':data['parent']['source_structure'],
                'protein_fragments_union':union,'sequence_audit':data['sequence_audit'],
                'contacts':data['contacts'],'atom_graph':mapping,
                'charge_ledger':ledger,'water_groups':groups,'atoms':atoms,'states':states,
                'atom_count':len(atoms),'charges':{'La':protein_charge+3,'Ca':protein_charge+2}}
        write_new(case_dir/'preparation.json',full); prepared.append(record(case_dir/'preparation.json'))
    implementation = out/'implementation'; implementation.mkdir()
    pins = {}
    for name in ('hydration_network.py','hydration_square.py','affordable_peptide.py'):
        dst = implementation/name; shutil.copyfile(Path(__file__).with_name(name),dst); pins[name]=record(dst)
    result = {'schema_version':'alquemia.hydration_network.v1','protocol_id':PROTOCOL,
              'configuration':record(config),'agreement':cfg['agreement'],'source_config':cfg['source_config'],
              'preparations':prepared,'tasks':tasks,'implementation':pins,'orca':source['orca'],
              'execution_policy':source['execution_policy'],'occupancy_enumerated':occupancy,'case_filter':case_filter,
              'full_patterns_excluded_for_separate_collection':exclude_full,
              'physical_scope':'variable-water orientations; oxygen/protein/outer-water coordinates fixed',
              'free_energy_status':'missing_bound_basin_and_nonelectrostatic_terms',
              'absolute_occupancy':None,'classifier':None,'baseline_changed':False,
              'cost_tracking':{'compute_budget':None,'wall_time_limit':None,'task_count':len(tasks)}}
    write_new(out/'manifest.json',result)
    return {'manifest':record(out/'manifest.json'),'tasks':len(tasks),
            'preparations':[{'case':read_json(p['path'])['case'],'atoms':read_json(p['path'])['atom_count']} for p in prepared]}


def optimization_trace(text, initial):
    """Read actual printed total gradients, each with its evaluated coordinates.

    ORCA's final energy-only evaluation is not the last gradient evaluation.
    Printed coordinates have six decimals; retain that precision limitation.
    """
    n=len(initial); elements=[a[0] for a in initial]
    coordinate_tables=[]
    for match in re.finditer(r'^CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n-+\s*\n',text,re.M):
        rows=[]
        for line in text[match.end():].splitlines():
            fields=line.split()
            if len(fields)!=4:break
            try:rows.append((fields[0],*[float(v) for v in fields[1:]]))
            except ValueError:break
        if [a[0] for a in rows]!=elements:
            raise InvalidArtifact('printed optimization coordinate ordering/count differs')
        coordinate_tables.append((match.start(),rows))
    energy_matches=list(re.finditer(r'FINAL SINGLE POINT ENERGY\s+([-+\d.Ee]+)',text))
    gradients=[]
    for match in re.finditer(r'^CARTESIAN GRADIENT\s*\n-+\s*\n\s*\n',text,re.M):
        rows=[]
        for line in text[match.end():].splitlines():
            row=re.fullmatch(r'\s*(\d+)\s+(\w+)\s*:\s*([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s*',line)
            if not row:break
            rows.append((int(row[1]),row[2],[float(row[i]) for i in (3,4,5)]))
        if len(rows)!=n or [r[0] for r in rows]!=list(range(1,n+1)) or [r[1] for r in rows]!=elements:
            raise InvalidArtifact('printed total gradient ordering/count differs')
        energies=[e for e in energy_matches if e.start()<match.start()]
        coordinates=[c for c in coordinate_tables if c[0]<match.start()]
        if not energies or not coordinates or coordinates[-1][0]>energies[-1].start():
            raise InvalidArtifact('gradient lacks a matched preceding energy/geometry')
        values=np.array([r[2] for r in rows]);coords=np.array([r[1:] for r in coordinates[-1][1]])
        if not np.isfinite(values).all() or not np.isfinite(coords).all():
            raise InvalidArtifact('nonfinite optimization trace')
        gradients.append({'energy_hartree':float(energies[-1][1]),'coordinates_A':coords.tolist(),
            'gradient_Ha_per_bohr':values.tolist(),'output_line':text.count('\n',0,match.start())+1})
    if not gradients or not coordinate_tables:
        raise InvalidArtifact('no actual analytic optimization gradient available')
    return {'analytic_gradient_evaluations':gradients,
            'coordinate_print_resolution_A':1e-6,'gradient_print_resolution_Ha_per_bohr':1e-9,
            'last_printed_coordinates_A':[r[1:] for r in coordinate_tables[-1][1]],
            'final_endpoint_gradient':None,
            'gradient_scope':'evaluated optimization iterates; not assigned to the later final energy-only geometry'}


def geometry_checks(initial, final, groups, mobile, *, enforce=True):
    if [a[0] for a in initial]!=[a[0] for a in final]:
        raise InvalidArtifact('final XYZ atom ordering differs')
    coords=np.array([a[1:] for a in final]);fixed=[i for i in range(len(initial)) if i not in mobile]
    if not np.isfinite(coords).all():raise InvalidArtifact('nonfinite final coordinates')
    frozen_error=float(np.max(np.abs(coords[fixed]-np.array([initial[i][1:] for i in fixed])))) if fixed else 0.
    if enforce and frozen_error>2e-5:raise InvalidArtifact('frozen coordinates moved')
    geometry_error=0.
    for w in groups:
        ids=w['indices'];before=np.array([initial[i][1:] for i in ids]);after=coords[ids]
        geometry_error=max(geometry_error,float(np.max(np.abs(
            np.linalg.norm(before[:,None]-before[None,:],axis=2)-np.linalg.norm(after[:,None]-after[None,:],axis=2)))))
    if enforce and geometry_error>2e-4:raise InvalidArtifact('rigid water geometry changed')
    return frozen_error,geometry_error


def collect(manifest, output):
    m=read_json(manifest); results=[]
    for t in m['tasks']:
        e=None
        try:
            op=Path(t['output_path']); rp=Path(str(op)+'.execution.json')
            e=square.endpoint(record(op),record(rp),t['xyz'],t['input'])
            text=op.read_text()
            if re.search(r'numerical gradient|numerical differentiation',text,re.I):
                raise InvalidArtifact('numerical derivative was used')
            opt=t['task_type']=='constrained_optimization'
            if opt and 'THE OPTIMIZATION HAS CONVERGED' not in text:
                raise InvalidArtifact('geometry optimization did not converge')
            e['optimization_converged']=opt
            initial=xyz(verify(t['xyz']))
            if opt:
                final=op.parent/'endpoint.runtime.xyz'
                final_atoms=xyz(final)
                trace=optimization_trace(text,initial)
                e.update(optimization_trace=trace,final_xyz=record(final))
                frozen_error,geometry_error=geometry_checks(initial,final_atoms,t['groups'],t['mobile_indices'],enforce=False)
                e.update(frozen_coordinate_error_A=frozen_error,water_geometry_error_A=geometry_error)
                geometry_checks(initial,final_atoms,t['groups'],t['mobile_indices'])
                if np.max(np.abs(np.array(trace['last_printed_coordinates_A'])-np.array([a[1:] for a in final_atoms])))>2e-6:
                    raise InvalidArtifact('final XYZ differs from last printed geometry')
                comment_energy=re.search(r'\bE\s+([-+\d.Ee]+)',final.read_text().splitlines()[1])
                if not comment_energy or abs(float(comment_energy[1])-e['energy_hartree'])>1e-9:
                    raise InvalidArtifact('final XYZ energy label differs')
                receipt=read_json(rp)
                e.update(optimization_trace=trace,final_xyz=record(final),frozen_coordinate_error_A=frozen_error,
                         water_geometry_error_A=geometry_error,
                         runner_expected_engrad_complete=receipt.get('expected_engrad_complete'),
                         final_endpoint_gradient_status='not_requested_by_native_Opt')
            energies=[float(x) for x in re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+\d.Ee]+)',text)]
            e.update(initial_energy_hartree=energies[0],scf_energy_evaluation_count=len(energies),
                     relaxation_kcal_mol=(e['energy_hartree']-energies[0])*HA_TO_KCAL,
                     optimization_converged=opt)
            row={'task_id':t['task_id'],'case':t['case'],'pattern':t['pattern'],'seed':t['seed'],
                 'metal':t['metal'],'variable_water_count':t['variable_water_count'],'status':'complete','result':e}
        except (InvalidArtifact,FileNotFoundError,KeyError,ValueError) as err:
            row={'task_id':t['task_id'],'case':t['case'],'pattern':t['pattern'],'seed':t['seed'],
                 'metal':t['metal'],'status':'unavailable','failure':str(err),
                 'observed_endpoint_not_geometry_qualified':e}
        results.append(row)
    result={'manifest':record(manifest),'rows':results,'status':'complete' if all(r['status']=='complete' for r in results) else 'incomplete',
            'occupancy_probabilities':None,'calibrated_score':None,'baseline_changed':False,'collector':record(__file__)}
    write_new(output,result)
    return {'status':result['status'],'completed':sum(r['status']=='complete' for r in results),'tasks':len(results)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='operation',required=True)
    q=sub.add_parser('prepare'); q.add_argument('--config',required=True);q.add_argument('--output',required=True);q.add_argument('--occupancy',action='store_true');q.add_argument('--case');q.add_argument('--exclude-full',action='store_true')
    q=sub.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    args=vars(p.parse_args());op=args.pop('operation')
    import json
    print(json.dumps(globals()[op](**args),indent=2))
