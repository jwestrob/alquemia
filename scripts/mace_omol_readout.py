"""Native atomic-energy bookkeeping on four unchanged, previously scored GGR states."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz, paired
from affordable_response import source_key
from mace_hybrid import EV_TO_KCAL

CASES = ('GGR_extended', 'GGR_connected')
DEFINITION = 'native_node_energy_plus_optional_embedding_readout_v1'


def source(collection):
    from mace_omol import validate, collect
    saved = read_json(collection); mp = verify(saved['manifest']); m = read_json(mp)
    if m['stage'] != 'benchmark':
        raise InvalidArtifact('readout requires the completed original OMOL benchmark')
    validate(mp)
    if saved != collect(mp) or saved['status'] != 'complete':
        raise InvalidArtifact('source benchmark differs from executed receipts or is incomplete')
    return saved, m


def tasks(parent):
    result = []
    for name in CASES:
        for metal in ('La', 'Ca'):
            candidates = [t for t in parent['tasks'] if t['task_id'] == name+'_'+metal]
            if len(candidates) != 1:
                raise InvalidArtifact('missing/duplicate frozen GGR endpoint')
            t = copy.deepcopy(candidates[0]); t.pop('cache_key')
            t['capture_native_readout'] = True
            result.append(t)
    return result


def prepare(collection, agreement, output):
    from mace_omol import common, seal
    _, parent = source(collection)
    _, out, m = common(verify(parent['inventory']), verify(parent['software']), agreement, output, 'readout')
    m.update(source_collection=record(collection), tasks=tasks(parent))
    # Verify the already prepared coordinate/source identities before scheduling.
    physical_mapping(m['tasks'])
    return seal(out, m)


def physical_mapping(task_list):
    maps = {}; sources = []
    for name in CASES:
        t = {r['metal']: r for r in task_list if r['case_id'] == name}
        if set(t) != {'La', 'Ca'} or t['La']['preparation'] != t['Ca']['preparation']:
            raise InvalidArtifact('paired GGR preparation differs')
        mapping = read_json(verify(t['La']['preparation']))
        prep = read_json(verify(mapping['source_preparation']))
        sources.append(prep['source_structure']); verify(prep['source_structure'])
        paired(verify(t['La']['xyz']), verify(t['Ca']['xyz']), t['La']['charge'], t['Ca']['charge'])
        atoms = xyz(verify(t['La']['xyz'])); rows = prep['atom_graph']['source_to_qm']
        if (len(rows)+1 != len(atoms) or sorted(r['qm_index'] for r in rows) != list(range(1,len(atoms)))):
            raise InvalidArtifact('source map is not bijective')
        all_rows = [{'qm_index':0, 'kind':'metal', 'id':'metal', 'element':'M',
                     'xyz_A':list(atoms[0][1:]), 'source':prep['selected_site']}]
        for row in sorted(rows, key=lambda r:r['qm_index']):
            i = row['qm_index']; r = copy.deepcopy(row)
            if row['kind'] == 'source':
                r.update(id=source_key(row['source']), element=row['source']['element'])
            elif row['kind'] == 'sigma_link_H':
                r.update(id='cap:'+source_key(row['retained'])+'->'+source_key(row['omitted']), element='H')
            else:
                raise InvalidArtifact('unsupported source/cap mapping')
            if atoms[i][0] != r['element'] or np.max(np.abs(np.array(atoms[i][1:])-row['xyz_A'])) > 1e-6:
                raise InvalidArtifact('native coordinate no longer matches its source/cap map')
            all_rows.append(r)
        if len({r['id'] for r in all_rows}) != len(all_rows):
            raise InvalidArtifact('duplicate source identity')
        maps[name] = all_rows
    if sources[0] != sources[1]:
        raise InvalidArtifact('representations have different physical sources')
    lookup = [{r['id']:r for r in maps[name]} for name in CASES]
    shared = set(lookup[0]) & set(lookup[1])
    for key in shared:
        a,b = (d[key] for d in lookup)
        if a['kind'] != b['kind'] or a['element'] != b['element'] or np.max(np.abs(np.array(a['xyz_A'])-b['xyz_A'])) > 1e-6:
            raise InvalidArtifact('shared atom identity has incompatible coordinates or chemistry')
    return maps, shared


class NativeCapture:
    """Observe one unmodified model forward; no additional inference or forces."""
    def __init__(self, model):
        self.outputs = []; self.embedding = []
        self.has_embedding = hasattr(model, 'embedding_readout')
        self.handles = [model.register_forward_hook(self.capture)]
        if self.has_embedding:
            self.handles.append(model.embedding_readout.register_forward_hook(self.capture_embedding))

    @staticmethod
    def array(tensor):
        return tensor.detach().cpu().numpy().copy()

    def capture(self, module, inputs, output):
        self.outputs.append({k:self.array(output[k]) for k in ('energy','node_energy','interaction_energy')})

    def capture_embedding(self, module, inputs, output):
        self.embedding.append(self.array(output))

    def save(self, calc, total, output):
        from mace_omol import TOL
        for h in self.handles: h.remove()
        if len(self.outputs) != 1 or len(self.embedding) != int(self.has_embedding):
            raise InvalidArtifact('unexpected number of native or embedding forwards')
        native = self.outputs[0]; node = native['node_energy'].reshape(-1)
        embedding = self.embedding[0].reshape(-1) if self.has_embedding else np.zeros(len(node))
        ase_node = np.asarray(calc.results['energies'])
        e0 = ase_node - np.asarray(calc.results['node_energy'])
        if (embedding.shape != node.shape or not np.isfinite(node).all() or not np.isfinite(embedding).all()
                or not np.isfinite(e0).all() or not np.allclose(ase_node,node,rtol=0,atol=1e-10)):
            raise InvalidArtifact('unsupported native node/embedding shape, units or values')
        error = (float(node.sum()+embedding.sum())-total)*EV_TO_KCAL
        native_error = (float(native['energy'].item())-total)*EV_TO_KCAL
        if max(abs(error),abs(native_error)) > TOL['energy_kcal_mol']:
            raise InvalidArtifact('native atomic/graph readouts do not reproduce total energy')
        path = Path(output)/'native_readout_eV.npz'
        np.savez(path, node_energy_eV=node, embedding_energy_eV=embedding, atomic_reference_eV=e0)
        return {'definition':DEFINITION, 'arrays':record(path), 'forward_count':len(self.outputs),
                'embedding_readout_present':self.has_embedding,
                'embedding_absence_means':'architecturally_absent_term_not_missing_correction',
                'component_sum_error_kcal_mol':error, 'native_total_error_kcal_mol':native_error,
                'interaction_energy_eV':float(native['interaction_energy'].item()),
                'native_total_energy_eV':float(native['energy'].item())}


def checked_arrays(result, task):
    from mace_omol import TOL
    ref = result['native_readout']
    if ref['definition'] != DEFINITION or ref['forward_count'] != 1:
        raise InvalidArtifact('unsupported readout definition or forward count')
    with np.load(verify(ref['arrays']), allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    if set(arrays) != {'node_energy_eV','embedding_energy_eV','atomic_reference_eV'}:
        raise InvalidArtifact('missing native readout arrays')
    n = len(xyz(verify(task['xyz'])))
    if any(a.shape != (n,) or not np.isfinite(a).all() for a in arrays.values()):
        raise InvalidArtifact('readout atom count/nonfinite failure')
    if ref['embedding_readout_present'] is False and np.any(arrays['embedding_energy_eV'] != 0):
        raise InvalidArtifact('absent embedding term has nonzero contributions')
    total = float(arrays['node_energy_eV'].sum()+arrays['embedding_energy_eV'].sum())
    if abs((total-result['energy_eV'])*EV_TO_KCAL) > TOL['energy_kcal_mol']:
        raise InvalidArtifact('saved readouts no longer sum to total')
    return arrays


def report(collection, output):
    from mace_omol import validate, collect, TOL
    saved = read_json(collection); mp = verify(saved['manifest']); m = read_json(mp)
    validate(mp)
    if m['stage'] != 'readout' or saved != collect(mp) or saved['status'] != 'complete':
        raise InvalidArtifact('complete actual readout collection required')
    old, _ = source(verify(m['source_collection'])); maps, shared = physical_mapping(m['tasks'])
    checks = []; data = {}; endpoint_terms = {}
    for t in m['tasks']:
        key=t['task_id']; new=saved['rows'][key]; baseline=old['rows'][key]
        arrays=checked_arrays(new,t); data[key]=arrays
        de=(new['energy_eV']-baseline['energy_eV'])*EV_TO_KCAL
        df=float(np.max(np.abs(np.load(verify(new['forces']))-np.load(verify(baseline['forces'])))))
        checks.append({'name':key+'_unchanged_endpoint','energy_error_kcal_mol':de,
                       'force_error_eV_A':df,'pass':abs(de)<=TOL['energy_kcal_mol'] and df<=TOL['force_max_eV_A']})
        endpoint_terms[key] = {'energy_eV':new['energy_eV'],'node_sum_eV':float(arrays['node_energy_eV'].sum()),
                              'embedding_sum_eV':float(arrays['embedding_energy_eV'].sum()),
                              'atomic_reference_sum_eV':float(arrays['atomic_reference_eV'].sum()),
                              'capture':new['native_readout']}
    per_atom=[]; groups={}; contrasts={}; embedding_deltas={}
    categories=('metal','shared_source_atoms','unique_source_atoms','synthetic_caps','embedding_graph_term')
    for name in CASES:
        ca,la=(data[name+'_'+metal] for metal in ('Ca','La'))
        embedding_deltas[name]=(ca['embedding_energy_eV']-la['embedding_energy_eV'])*EV_TO_KCAL
        delta=(ca['node_energy_eV']-la['node_energy_eV'])*EV_TO_KCAL
        group={key:0. for key in categories}
        for row, value in zip(maps[name],delta):
            category=('metal' if row['kind']=='metal' else 'synthetic_caps' if row['kind']=='sigma_link_H'
                      else 'shared_source_atoms' if row['id'] in shared else 'unique_source_atoms')
            group[category]+=float(value)
            per_atom.append({'case_id':name,**row,'category':category,'R_node_kcal_mol':float(value)})
        group['embedding_graph_term']=float((ca['embedding_energy_eV']-la['embedding_energy_eV']).sum()*EV_TO_KCAL)
        contrasts[name]=(saved['rows'][name+'_Ca']['energy_eV']-saved['rows'][name+'_La']['energy_eV'])*EV_TO_KCAL
        error=sum(group.values())-contrasts[name]
        checks.append({'name':name+'_contrast_sum','error_kcal_mol':error,'pass':abs(error)<=TOL['energy_kcal_mol']})
        groups[name]=group
    delta_groups={key:groups[CASES[1]][key]-groups[CASES[0]][key] for key in categories}
    shift=contrasts[CASES[1]]-contrasts[CASES[0]]
    error=sum(delta_groups.values())-shift
    checks.append({'name':'representation_shift_sum','error_kcal_mol':error,'pass':abs(error)<=TOL['energy_kcal_mol']})
    # Linear embedding readout of h_Z + f(Q,spin): sum_i a_Zi + N*b(Q,spin).
    # This is an architectural identity, not a fitted correction or physical energy claim.
    nonmetal=np.concatenate([embedding_deltas[name][1:] for name in CASES])
    delta_n=len(maps[CASES[1]])-len(maps[CASES[0]])
    b=float(nonmetal[0]); predicted=delta_n*b
    algebra={'definition':'E_embedding=sum_i(a_element_i)+N*b(total_charge,multiplicity)',
             'atom_count_difference':delta_n,'per_atom_charge_state_contrast_kcal_mol':b,
             'nonmetal_contrast_range_kcal_mol':[float(nonmetal.min()),float(nonmetal.max())],
             'metal_embedding_contrast_difference_kcal_mol':float(embedding_deltas[CASES[1]][0]-embedding_deltas[CASES[0]][0]),
             'predicted_embedding_shift_kcal_mol':predicted,
             'observed_minus_predicted_kcal_mol':delta_groups['embedding_graph_term']-predicted,
             'fitted_parameters':False,'corrected_score_produced':False}
    result={'status':'complete','physical_cause_established':False,'score_correction_kcal_mol':None,
            'protocol_id':m['protocol_id'],'definition':DEFINITION,'sources':{'collection':record(collection),'original':m['source_collection']},
            'checks':checks,'all_checks_pass':all(c['pass'] for c in checks),'endpoint_terms':endpoint_terms,
            'R_kcal_mol':contrasts,'representation_shift_kcal_mol':shift,'grouped_R_kcal_mol':groups,
            'shift_by_group_kcal_mol':delta_groups,'per_atom':per_atom,'embedding_algebra':algebra,'baseline_changed':False,
            'interpretation':'model energy bookkeeping; atomic terms are not unique physical energies or a new score'}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/Path(__file__).name);result['implementation']=record(out/Path(__file__).name)
    write_new(out/'result.json',result)
    lines=['# GGR native OMOL readout diagnostic','',f'All replay/accounting checks pass: **{result["all_checks_pass"]}**.',
           f'Connected minus extended R: **{shift:.9f} kcal/mol**. No score has changed.','',
           '| Native readout group | Extended R | Connected R | Contribution to shift |','|---|---:|---:|---:|']
    for key in categories:
        lines.append(f'| {key} | {groups[CASES[0]][key]:.9f} | {groups[CASES[1]][key]:.9f} | {delta_groups[key]:.9f} |')
    lines+=['','Values are kcal/mol. Each graph term is retained separately; source hydrogens remain source atoms, not caps.',
            '',f'The explicit linear embedding readout has the form sum(a_element) + N*b(charge, spin). Its nonmetal Ca-minus-La term is {b:.12f} kcal/mol per atom. The extra {delta_n} atoms predict {predicted:.9f} kcal/mol of embedding shift; observed minus predicted is {algebra["observed_minus_predicted_kcal_mol"]:.3g}. This term does not depend on geometry. The remaining readouts still depend on charge-conditioned features and local context.',
            '', 'Atomic readouts are learned bookkeeping, not a unique physical energy partition. Identifying this extensive term does not establish that simply deleting it gives correct affinity. Both original scores remain in use; no corrected score is produced.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('collection','output'):p.add_argument('--'+key,required=True)
    result=report(**vars(p.parse_args()))
    print(json.dumps({k:result[k] for k in ('all_checks_pass','representation_shift_kcal_mol','shift_by_group_kcal_mol')},indent=2))
