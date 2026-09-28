"""Prepare-only exact-source C-bound hydrogen transfer; no electronic evaluation."""
import argparse,json
from pathlib import Path
import numpy as np
from affordable_common import read_json,record,verify,xyz,write_new,InvalidArtifact

def transfer_source(source,design,repair_map,review_record,output):
    review=read_json(verify(review_record))
    # The caller supplies a separately reviewed admitted source record, not old RESULT.json.
    if review.get('source_id')!=source or review.get('admitted') is not True:
        raise InvalidArtifact('Chemistry-aware source admission required')
    if not all(review.get('checks',{}).values()) or review.get('unresolved_inversions'):raise InvalidArtifact('Review checks unavailable or failed')
    recovered=read_json(verify(review['recovered_record']))
    config=read_json(verify(recovered['inputs']));export=read_json(verify(config['export']))
    atoms=read_json(verify(export['atoms']));bonds=read_json(verify(export['bonds']));saved=read_json(verify(recovered['saved_coordinates']))
    ids={a['id']:i for i,a in enumerate(atoms)}
    if saved['source_ids']!=[a['id'] for a in atoms]:raise InvalidArtifact('Repair coordinate ID order differs')
    coords=np.asarray(saved['xyz_A'],float)
    if coords.shape!=(len(atoms),3) or not np.isfinite(coords).all():raise InvalidArtifact('Repair coordinates invalid')
    adj={k:set() for k in ids}
    for a,b in bonds:adj[a].add(b);adj[b].add(a)
    free=set(recovered['free_H_ids']);rm=read_json(verify(repair_map));maps=rm['atom_graph']['source_to_qm'];targets={};unchanged_real_H=[];originals=design['sources'][source]['xyz']
    old={m:xyz(verify(p)) for m,p in originals.items()}
    for entry in maps:
        if 'source' not in entry:continue
        a=entry['source'];sid=f"{a['chain']}/{a['resnum']}/{a.get('insertion_code','')}/{a['atom']}";qi=entry['qm_index']
        if not np.allclose(old['La'][qi][1:],entry['xyz_A'],rtol=0,atol=1e-8):raise InvalidArtifact('Old compact source map differs')
        if sid not in ids:
            if a['element']=='H':unchanged_real_H.append(dict(index=qi,source_id=sid,reason='outside protein export'))
            continue
        ia=ids[sid]
        if atoms[ia]['element']!=a['element']:raise InvalidArtifact('Element mapping differs')
        if a['element']!='H':
            if not np.allclose(old['La'][qi][1:],coords[ia],rtol=0,atol=1e-8):raise InvalidArtifact('Heavy source frame changed')
            continue
        is_ch=len(adj[sid])==1 and atoms[ids[next(iter(adj[sid]))]]['element']=='C'
        if is_ch:
            if sid not in free:raise InvalidArtifact('C-H missing admitted free set')
            targets[qi]=dict(source_id=sid,xyz_A=coords[ia].tolist(),parent_id=next(iter(adj[sid])))
        else:unchanged_real_H.append(dict(index=qi,source_id=sid,reason='exchangeable or water H retained exactly'))
    if not targets:raise InvalidArtifact('No mapped C-bound H to transfer')
    out=Path(output)/source;out.mkdir(parents=True,exist_ok=False);newrefs={};checks={}
    for metal,ref in originals.items():
        lines=verify(ref).read_text().splitlines();new=lines.copy()
        for qi,v in targets.items():
            if old[metal][qi][0]!='H':raise InvalidArtifact('Transfer target is not H')
            new[qi+2]='H '+' '.join(f'{x:.10f}' for x in v['xyz_A'])
        path=out/(metal+'.xyz');path.write_text('\n'.join(new)+'\n');actual=xyz(path)
        fixed=[i for i in range(len(actual)) if i not in targets]
        checks[metal+'_fixed_lines_exact']=all(lines[i+2]==new[i+2] for i in fixed)
        checks[metal+'_inventory_exact']=[a[0] for a in old[metal]]==[a[0] for a in actual]
        checks[metal+'_transferred_coordinates_match']=all(np.allclose(actual[i][1:],v['xyz_A'],rtol=0,atol=6e-11) for i,v in targets.items())
        newrefs[metal]=record(path)
    aa={m:xyz(verify(p)) for m,p in newrefs.items()};checks['paired_coordinates_exact']=np.array_equal([a[1:] for a in aa['La']],[a[1:] for a in aa['Dy']]);checks['all_cap_lines_unchanged']=all(all(verify(originals[m]).read_text().splitlines()[c['qm_index']+2]==verify(newrefs[m]).read_text().splitlines()[c['qm_index']+2] for c in rm['atom_graph']['cut_bonds_and_caps']) for m in ('La','Dy'))
    if not all(checks.values()):raise InvalidArtifact('Transfer invariants failed')
    for i,v in targets.items():v['old_xyz_A']=old['La'][i][1:];v['displacement_A']=float(np.linalg.norm(np.array(v['xyz_A'])-v['old_xyz_A']))
    return dict(source_id=source,xyz=newrefs,original_xyz=originals,source_map=repair_map,reviewed_admission=review_record,repaired_H_indices=sorted(targets),transferred_atoms=[dict(qm_index=i,**v) for i,v in sorted(targets.items())],retained_real_H=unchanged_real_H,checks=checks,charge=-1,physical_states={'La':'LaIII 4f0 singlet','Dy':'DyIII 4f9 sextet'},effective_multiplicity=1,preparation='targeted_C_bound_H_only',energy_evaluations=0)

def main():
    p=argparse.ArgumentParser();p.add_argument('--request',required=True);p.add_argument('--output',required=True);a=p.parse_args();r=read_json(a.request);design=read_json(verify(r['design']));out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
    if {x['source_id'] for x in r['sources']}!={'Hans8DQ2','Hans8FNR','Mex8FNS'}:raise InvalidArtifact('All three sources required')
    rows=[transfer_source(x['source_id'],design,x['source_map'],x['reviewed_admission'],out) for x in r['sources']]
    write_new(out/'manifest.json',dict(protocol_id='nikasha_compact_CH_repair_preparation_v1',request=record(a.request),implementation=record(__file__),sources=rows,electronic_endpoints_prepared=6,electronic_evaluations=0,limitations='Old compact exchangeable H/waters/caps retained; not generic normalized preparation'))
if __name__=='__main__':main()
