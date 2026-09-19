"""Fixed-coordinate donor-neighbor context pilot; production cores stay immutable."""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import shutil
import numpy as np
import affordable_peptide as peptide
from affordable_common import (HA_TO_KCAL, InvalidArtifact, cache_key, read_json,
                               record, verify, write_new, xyz)
from hydration_network import source_id
from hydration_square import endpoint
from mace_hybrid import check_atoms, write_xyz

PROTOCOL = 'native_r2scan3c_cpcm_fixed_donor_second_shell_v1'
POLICY = {'direct_donor_cutoff_A': 3.3, 'polar_neighbor_cutoff_A': 3.5,
          'polar_elements': ['N', 'O', 'S'], 'new_water_policy': 'exclude_record_keep_original_inventory',
          'geometry': 'fixed_all_shared_source_and_surviving_cap_coordinates',
          'selection': 'direct_nonwater_ON_and_same_carbonyl_carboxylate_partners',
          'overlaps': 'source_graph_union_complete_amide_sigma_caps'}
CASES = ('1H4I', '4MAE', '1GLG', '1F6S', '6IP9')
EXTRA_SIDECHAINS = {
    'ARG': (('CB','CG','CD','NE','CZ','NH1','NH2'), {'CB':2,'CG':2,'CD':2,'NE':1,'NH1':2,'NH2':2}, 1),
    'LYS': (('CB','CG','CD','CE','NZ'), {'CB':2,'CG':2,'CD':2,'CE':2,'NZ':3}, 1),
    'TRP': (('CB','CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'),
            {'CB':2,'CD1':1,'NE1':1,'CE3':1,'CZ2':1,'CZ3':1,'CH2':1}, 0)}


def atom_key(meta):
    return meta['chain_index'], meta['residue_index'], meta['atom']


def cap_key(item):
    return atom_key(item['retained']), atom_key(item['omitted'])


def sidechain(graph, rkey):
    r = graph.residues[rkey]; rn = r.canonical_resname
    if rn in peptide.old.SIDECHAIN_QM_HEAVY_ATOMS:
        _, audit = peptide.old.sidechain_fragment(r, [], 0)
        names = peptide.old.SIDECHAIN_QM_HEAVY_ATOMS[rn]; charge = audit['formal_charge']
    elif rn in EXTRA_SIDECHAINS:
        names, expected, charge = EXTRA_SIDECHAINS[rn]
        observed = Counter(p[2] for h,p in graph.hparents.items() if p[:2] == rkey and p[2] in names)
        if dict(observed) != expected:
            raise InvalidArtifact(f'unsupported sidechain H/protonation: {r.source_dict()} {dict(observed)}')
    else:
        raise InvalidArtifact(f'unsupported second-shell sidechain {r.source_dict()}')
    return {graph.key(rkey, n) for n in names}, charge


def default_config(root, agreement):
    root = Path(root).resolve()
    hold = root/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/prepared'
    paths = {c:next(hold.glob('*'+c+'/*carve_manifest.json')) for c in CASES[:2]}
    paths['1GLG'] = root/'workspaces/affordable_challenger_20260915/verified_repairs/ggr_1glg_GGR/repair_manifest.json'
    for c in CASES[3:]:
        paths[c] = root/f'workspaces/benchmark_set_20260915/prepared/alacta_{c.lower()}_v2/strong_site/amide_v3/repair_manifest.json'
    transfer = root/'workspaces/hydration_network_20260918/core_transfer_v1/manifest.json'
    tm = read_json(transfer)
    ggr_release = read_json(verify(tm['ggr_release']))
    ggr = next(r for r in ggr_release['scores'] if r['case']=='ggr_1glg_GGR' and r['lane']=='repaired')
    entries=[]
    for c,p in paths.items():
        d = read_json(p); eps={}
        for metal in ('Ca','La'):
            if c in CASES[:2]:
                xp=d['outputs'][metal+'_xyz']; ip=d['outputs'][metal+'_input']; op=Path(ip['path']).with_suffix('.out')
                rr=record(str(op)+'.execution.json'); q=d['charge_ledger'][metal+'_total']
            elif c=='1GLG':
                xp=d['outputs'][metal]['xyz']; ip=d['outputs'][metal]['input']; op=verify(ggr['endpoints'][metal]['output'])
                rr=ggr['endpoints'][metal]['receipt'];q=d['outputs'][metal]['charge']
            else:
                t=next(t for t in tm['tasks'] if t['case']==c and t['metal']==metal)
                xp=t['xyz'];ip=t['input'];op=Path(t['output_path']);rr=record(str(op)+'.execution.json');q=t['charge']
            parsed=endpoint(record(op),rr,xp,ip)
            eps[metal]={'xyz':xp,'input':ip,'output':record(op),'receipt':rr,'charge':q,'energy_hartree':parsed['energy_hartree']}
        entries.append({'case':c,'parent':record(p),'endpoints':eps,
                        'group':'alpha_lactalbumin' if c in CASES[3:] else {'1GLG':'GGR','1H4I':'MxaF','4MAE':'XoxF'}[c],
                        'evidence_stratum':'functional_PQQ_class' if c in CASES[:2] else 'condition_qualified_affinity_direction'})
    return {'protocol_id':PROTOCOL,'policy':POLICY,'cases':entries,'agreement':record(agreement),
            'topology':read_json(paths['1GLG'])['topology_definition'],
            'source_water_preparation':record(transfer),
            'orca':tm['orca'],'execution_policy':tm['execution_policy'],
            'software':record(root/'workspaces/mace_omol_20260917/software_v1/software_manifest.json')}


def parent_state(item, topology):
    p=read_json(verify(item['parent'])); isp='fixed_core' in p
    source=read_json(verify(p['protonation_manifest']))['output'] if isp else p['source_structure']
    graph=peptide.SourceGraph(verify(source), verify(topology))
    original={m:xyz(verify(e['xyz'])) for m,e in item['endpoints'].items()}
    # Exact shared core identity; alpha endpoints differ in declared water H only.
    if [a[0] for a in original['Ca'][1:]] != [a[0] for a in original['La'][1:]]:
        raise InvalidArtifact('paired core composition differs')
    maps=[];opaque=[];charges={};selected=set()
    if isp:
        index=1
        for f in p['qm_fragments']:
            match=re.fullmatch(r'([^:]+):([A-Z]+)(-?\d+)([A-Za-z]?)',f['id'])
            if not match:raise InvalidArtifact('unrecognized pinned PQQ fragment selector')
            ch,rn,num,ic=match.groups();r=graph.locate(dict(chain=ch,resname=rn,resnum=int(num),insertion_code=ic))
            if rn!='PQQ':
                heavy,charge=sidechain(graph,r.key);selected |= heavy
                if charge!=f['formal_charge']:raise InvalidArtifact('parent sidechain charge mismatch')
                charges[r.key]=charge
            for a in f['atom_records']:
                row=original['La'][index]
                if row[0]!=a['element'] or np.max(np.abs(np.array(row[1:])-a['xyz_A']))>1e-8:
                    raise InvalidArtifact('PQQ fragment coordinate/order mismatch')
                if rn=='PQQ':
                    entry={'qm_index':index,'kind':'opaque_cofactor','fragment':f['id'],'name':a['name'],'element':a['element']}
                    if a['origin'].startswith('source_'):entry.update(source=graph.meta[graph.key(r.key,a['name'])])
                    opaque.append(entry)
                elif a['origin'].startswith('source_'):
                    maps.append({'qm_index':index,'kind':'source','source':graph.meta[graph.key(r.key,a['name'])]})
                elif a['origin']=='generated_Cbeta_link_cap':
                    maps.append({'qm_index':index,'kind':'sigma_link_H','retained':graph.meta[graph.key(r.key,'CB')],
                                 'omitted':graph.meta[graph.key(r.key,'CA')]})
                else:raise InvalidArtifact('unrecognized PQQ-core hydrogen representation')
                index+=1
        if index != len(original['La']):raise InvalidArtifact('incomplete parent PQQ map')
    else:
        maps=p['atom_graph']['source_to_qm']
        for a in maps:
            if a['kind']=='source' and a['source']['element'] not in ('H','D') and a['source']['canonical_resname']!='HOH':
                selected.add(atom_key(a['source']))
        for row in p['charge_ledger']:
            if row['kind']=='sidechain':
                r=graph.locate(row['source']);_,charge=sidechain(graph,r.key)
                if charge!=row['formal_charge']:raise InvalidArtifact('parent sidechain charge mismatch')
                charges[r.key]=charge
    old_sources={atom_key(a['source']):a['qm_index'] for a in maps+opaque if 'source' in a}
    if len(old_sources)!=sum('source' in a for a in maps+opaque):raise InvalidArtifact('duplicate parent source mapping')
    water_keys={k for k in old_sources if graph.meta[k]['canonical_resname']=='HOH'}
    for metal in ('Ca','La'):
        e=item['endpoints'][metal]
        actual=endpoint(e['output'],e['receipt'],e['xyz'],e['input'])
        if actual['energy_hartree']!=e['energy_hartree']:raise InvalidArtifact('archived energy differs')
        header=[s.strip() for s in verify(e['input']).read_text().splitlines() if s.strip().startswith('!')]
        if len(header)!=1 or set(header[0].lower().split()[1:])!={'r2scan-3c','noautostart','cpcm(water)','defgrid3'}:
            raise InvalidArtifact('parent native SP recipe unsupported')
    return dict(parent=p,source=source,graph=graph,original=original,maps=maps,opaque=opaque,
                selected=selected,charges=charges,old_sources=old_sources,water_keys=water_keys)


def discover(state):
    g=state['graph'];old=state['old_sources'];rows=state['original']['La'];center=np.array(rows[0][1:])
    anchors={k for k,i in old.items() if rows[i][0] in ('N','O') and k not in state['water_keys']
             and np.linalg.norm(np.array(rows[i][1:])-center)<=POLICY['direct_donor_cutoff_A']}
    for k in list(anchors):
        for mid in g.required[k]:
            if mid in g.atoms and g.atoms[mid].element.name=='C':
                anchors |= {z for z in g.required[mid] if z in old and g.atoms[z].element.name in ('O','N')}
        # PQQ is opaque to the standard-protein topology; its named carboxyl pair
        # is nevertheless explicitly present in the pinned complete cofactor.
        if g.meta[k]['canonical_resname']=='PQQ' and re.fullmatch(r'O[279][AB]',k[2]):
            sibling=(*k[:2], k[2][:-1]+('B' if k[2][-1]=='A' else 'A'))
            if sibling not in old:raise InvalidArtifact('incomplete frozen PQQ carboxyl group')
            anchors.add(sibling)
    if not anchors:raise InvalidArtifact('no source-mapped donor functional group')
    fragments=set();contacts=[];excluded=[]
    for k,atom in g.atoms.items():
        if atom.element.name not in POLICY['polar_elements'] or k in old:continue
        distances={a:float(atom.pos.dist(g.atoms[a].pos)) for a in anchors}
        nearest=min(distances,key=distances.get);distance=distances[nearest]
        if distance>POLICY['polar_neighbor_cutoff_A']:continue
        meta=g.meta[k];contact={'atom':meta,'anchor':g.meta[nearest],'distance_A':distance}
        if meta['canonical_resname']=='HOH':
            excluded.append(dict(contact,reason='fixed_explicit_water_inventory'));continue
        if k[:2] not in g.supported_topology_residues:
            raise InvalidArtifact(f'unsupported nearby polar species: {meta}')
        if k[2] in ('O','OXT','C'):frag=('peptide',k[:2])
        elif k[2]=='N':
            prev=g.amide_prev.get(k[:2])
            if prev is None:raise InvalidArtifact(f'absent bonded predecessor of {meta}')
            frag=('peptide',prev[:2])
        else:frag=('sidechain',k[:2])
        fragments.add(frag);contacts.append(contact)
    return anchors,fragments,contacts,excluded


def expansion(state):
    g=state['graph'];anchors,fragments,contacts,excluded=discover(state)
    selected=set(state['selected']);charges=dict(state['charges']);added=[]
    for kind,rkey in sorted(fragments):
        if kind=='peptide':selected |= g.peptide(rkey);charge=0
        else:
            nodes,charge=sidechain(g,rkey);selected |= nodes;charges[rkey]=charge
        added.append({'kind':kind,'source':g.residues[rkey].source_dict(),'formal_charge':charge})
    atoms,mapping=g.materialize(selected)
    # materialize numbers atoms starting at1 for metal0. Copy original cofactor
    # and complete existing waters after graph materialization, keeping exact H.
    for a in state['opaque']:
        entry=dict(a,qm_index=len(atoms)+1);mapping['source_to_qm'].append(entry)
        atoms.append(state['original']['La'][a['qm_index']])
    for a in state['maps']:
        if a['kind']=='source' and atom_key(a['source']) in state['water_keys']:
            mapping['source_to_qm'].append(dict(a,qm_index=len(atoms)+1))
            atoms.append(state['original']['La'][a['qm_index']])
    src_lookup=state['old_sources'];cap_lookup={cap_key(a):a['qm_index'] for a in state['maps'] if a['kind']=='sigma_link_H'}
    opaque_lookup={(a['fragment'],a['name']):a['qm_index'] for a in state['opaque']}
    core_to_context={0:0};added_source=[]
    for a in mapping['source_to_qm']:
        if a['kind']=='opaque_cofactor':oi=opaque_lookup[(a['fragment'],a['name'])]
        elif a['kind']=='source':
            key=atom_key(a['source']);oi=src_lookup.get(key)
            if oi is None:added_source.append(a['source'])
        else:oi=cap_lookup.get(cap_key(a))
        if oi is not None:core_to_context[oi]=a['qm_index']
    removed_caps=[a for a in state['maps'] if a['kind']=='sigma_link_H' and a['qm_index'] not in core_to_context]
    missing_source=set(src_lookup.values())-set(core_to_context)
    if missing_source:raise InvalidArtifact('expanded graph lost parent source atoms')
    paired_rows={}
    for metal,original in state['original'].items():
        rows=[original[0]]+list(atoms)
        for old,new in core_to_context.items():rows[new]=original[old]
        for a in mapping['source_to_qm']:a['xyz_A']=list(rows[a['qm_index']][1:])
        if any(rows[new]!=original[old] for old,new in core_to_context.items()):
            raise InvalidArtifact('shared core coordinate changed')
        coords=np.array([a[1:] for a in rows]);dist=np.linalg.norm(coords[:,None]-coords[None,:],axis=2)
        np.fill_diagonal(dist,np.inf)
        if float(dist.min())<.45:raise InvalidArtifact('overlapping expanded atoms/caps')
        paired_rows[metal]=rows
    dq=sum(charges.values())-sum(state['charges'].values())
    return paired_rows, {'anchors':[g.meta[k] for k in sorted(anchors)],'added_fragments':added,'contacts':contacts,
        'excluded_contacts':excluded,'mapping':mapping,'core_to_context':core_to_context,'removed_caps':removed_caps,
        'added_source_atoms':added_source,'added_formal_charge':dq,'source':state['source'],
        'new_atom_count':len(paired_rows['La']),'original_atom_count':len(state['original']['La']),
        'core_source_coordinates_unchanged':True,'water_inventory_unchanged':True,'mapping_coordinate_endpoint':'La',
        'paired_differing_nonmetal_indices':[i for i in range(1,len(paired_rows['La'])) if paired_rows['Ca'][i]!=paired_rows['La'][i]],
        'coordinate_note':'core source atoms and surviving caps exact; merging fragments replaces cut-bond caps with source atoms'}


def inventory(config, output):
    cfg=read_json(config);rows=[]
    if cfg['policy']!=POLICY or tuple(r['case'] for r in cfg['cases'])!=CASES:raise InvalidArtifact('scope changed')
    for item in cfg['cases']:
        try:
            s=parent_state(item,cfg['topology']);_,audit=expansion(s)
            rows.append({'case':item['case'],'status':'supported','audit':audit})
        except (ValueError,KeyError,OSError) as exc:
            rows.append({'case':item['case'],'status':'unsupported','reason':str(exc)})
    result={'configuration':record(config),'policy':POLICY,'rows':rows,'new_energy_calculations':0}
    write_new(output,result);return result


def prepare(config, output):
    cfg=read_json(config);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    if cfg['policy']!=POLICY or tuple(r['case'] for r in cfg['cases'])!=CASES:raise InvalidArtifact('scope changed')
    impl=out/'implementation';impl.mkdir();pins={}
    # Keep execution imports self-contained, including the established MACE worker.
    for p in Path(__file__).parent.glob('*.py'):
        dest=impl/p.name;shutil.copyfile(p,dest);pins[p.name]=record(dest)
    tasks=[];states=[];mtasks=[]
    import mace_omol as omol
    model=omol.model(verify(cfg['software']))
    for item in cfg['cases']:
        state=parent_state(item,cfg['topology']);rows,audit=expansion(state);case=item['case']
        audit_path=out/(case+'_preparation.json');write_new(audit_path,audit)
        states.append({'case':case,'preparation':record(audit_path),'source':item})
        for metal in ('Ca','La'):
            ep=item['endpoints'][metal];q=ep['charge']+audit['added_formal_charge'];tid=case+'__expanded__'+metal
            d=out/'tasks'/tid;d.mkdir(parents=True);xp=d/'core.xyz';write_xyz(xp,rows[metal]);ip=d/'endpoint.inp'
            ip.write_text('! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3\n'+f'* xyzfile {q} 1 core.xyz\n')
            check_atoms(rows[metal],q)
            task={'task_id':tid,'case':case,'metal':metal,'charge':q,'multiplicity':1,'input':record(ip),
                  'xyz':record(xp),'output_path':str(d/'endpoint.out'),'preparation':record(audit_path),'atom_count':len(rows[metal])}
            task['cache_key']=cache_key(task | {'policy':POLICY,'protocol_id':PROTOCOL,'orca':cfg['orca']});tasks.append(task)
            for variant,xpin,charge in [('core',ep['xyz'],ep['charge']),('expanded',record(xp),q)]:
                mtasks.append({'task_id':case+'__'+variant+'__'+metal,'case_id':case,'metal':metal,'metal_index':0,
                  'kind':'core','variant':variant,'energy_component':omol.COMPONENT,'energy_only':False,
                  'charge':charge,'spin_multiplicity':1,'xyz':xpin,'state':check_atoms(xyz(verify(xpin)),charge)})
    result={'schema_version':'alquemia.second_shell.v1','protocol_id':PROTOCOL,'agreement':cfg['agreement'],
      'configuration':record(config),'policy':POLICY,'states':states,'tasks':tasks,'implementation':pins,
      'orca':cfg['orca'],'execution_policy':cfg['execution_policy'],'execution_resources':{'mpi_ranks':16,'concurrent_tasks':4},
      'reference':None,'calibrated_decision':None,'baseline_changed':False,'evidence_use':'consumed_development',
      'new_DFT_endpoints':len(tasks),'compute_budget':None}
    write_new(out/'manifest.json',result)
    mm={'schema_version':'alquemia.mace_omol.v1','stage':'second_shell_fixed','protocol_id':'native_OMOL_second_shell_fixed_v1',
        'agreement':cfg['agreement'],'software':cfg['software'],'model':model,'implementation':pins,'tasks':mtasks,
        'source_DFT_manifest':record(out/'manifest.json'),'compute_budget':None}
    for t in mtasks:
        t['cache_key']=cache_key({'task':t,'model':model,'software':cfg['software'],'implementation':pins})
    write_new(out/'mace_manifest.json',mm)
    return validate(out/'manifest.json')


def validate(manifest):
    m=read_json(manifest);cfg=read_json(verify(m['configuration']));verify(m['agreement'])
    if m['policy']!=POLICY or m['protocol_id']!=PROTOCOL or len(m['tasks'])!=10:raise InvalidArtifact('frozen scope changed')
    for pin in m['implementation'].values():verify(pin)
    for s,item in zip(m['states'],cfg['cases']):
        source=parent_state(item,cfg['topology']);rows,audit=expansion(source);saved=read_json(verify(s['preparation']))
        if json.loads(json.dumps(audit))!=saved:raise InvalidArtifact('preparation regeneration differs')
        for metal in ('Ca','La'):
            t=next(t for t in m['tasks'] if t['case']==item['case'] and t['metal']==metal)
            actual=xyz(verify(t['xyz']))
            if actual!=rows[metal]:raise InvalidArtifact('exact expanded coordinates differ')
            check_atoms(actual,t['charge'])
    from affordable_workflow import dry_run
    return dry_run(manifest)


def collect(manifest, output):
    m=read_json(manifest);rows=[]
    for s in m['states']:
        item=s['source'];case=s['case'];eps={};failure=[]
        for metal in ('Ca','La'):
            t=next(t for t in m['tasks'] if t['case']==case and t['metal']==metal)
            try:
                op=Path(t['output_path']);eps[metal]=endpoint(record(op),record(str(op)+'.execution.json'),t['xyz'],t['input'])
            except (ValueError,OSError,KeyError) as exc:failure.append({'metal':metal,'reason':str(exc)})
        core=item['endpoints']['Ca']['energy_hartree']-item['endpoints']['La']['energy_hartree']
        r=eps['Ca']['energy_hartree']-eps['La']['energy_hartree'] if not failure else None
        prep=read_json(verify(s['preparation']))
        rows.append({'case':case,'group':item['group'],'status':'complete' if not failure else 'unavailable',
          'failures':failure,'endpoints':eps,'core_R_hartree':core,'expanded_R_hartree':r,
          'delta_R_kcal_mol':(r-core)*HA_TO_KCAL if r is not None else None,
          'added_formal_charge':prep['added_formal_charge'],'atoms_before':prep['original_atom_count'],'atoms_after':prep['new_atom_count'],
          'expanded_R_components_kcal_mol':{term:(eps['Ca']['components_hartree'][term]-eps['La']['components_hartree'][term])*HA_TO_KCAL
             if not failure and eps['Ca']['components_hartree'].get(term) is not None and eps['La']['components_hartree'].get(term) is not None else None
             for term in ('CPCM_dielectric','dispersion','gCP','SCF')}})
    lookup={r['case']:r for r in rows};comparisons=[]
    for a,b,label in [('4MAE','1H4I','PQQ functional association'),('1F6S','1GLG','condition-qualified affinity ordering'),('6IP9','1GLG','condition-qualified affinity ordering')]:
        x,y=lookup[a],lookup[b]
        comparisons.append({'La_like':a,'Ca_like':b,'evidence':label,'before_kcal_mol':(x['core_R_hartree']-y['core_R_hartree'])*HA_TO_KCAL,
         'after_kcal_mol':(x['expanded_R_hartree']-y['expanded_R_hartree'])*HA_TO_KCAL
           if x['expanded_R_hartree'] is not None and y['expanded_R_hartree'] is not None else None})
    result={'manifest':record(manifest),'status':'complete' if all(r['status']=='complete' for r in rows) else 'incomplete',
       'rows':rows,'comparisons':comparisons,'reference':None,'calibrated_decision':None,'baseline_changed':False,
       'interpretation':'new representation, charge/cavity changes retained; comparisons consumed, alpha replicas one biological group'}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    q=s.add_parser('config');q.add_argument('--root',required=True);q.add_argument('--agreement',required=True);q.add_argument('--output',required=True)
    for name in ('inventory','prepare'):
        q=s.add_parser(name);q.add_argument('--config',required=True);q.add_argument('--output',required=True)
    q=s.add_parser('validate');q.add_argument('--manifest',required=True)
    q=s.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('op')
    if op=='config':
        output=a.pop('output');r=default_config(**a);write_new(output,r)
    else:r=globals()[op](**a)
    print(json.dumps(r,indent=2))


if __name__=='__main__':main()
