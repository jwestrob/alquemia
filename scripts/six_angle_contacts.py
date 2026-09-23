"""Apply the existing omitted-protein contact flag to actual six-angle proposals."""
import argparse,json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from omitted_context_contacts import POLICY,source_state,within_three,metrics,atom_key
from mace_site_kinematics import Kinematics


def audit(collection,preparation,output):
    c=read_json(collection);lm=read_json(verify(c['manifest']));m=read_json(verify(lm['source_manifest']))
    topology=read_json(preparation)['config']['topology'];verify(topology);rows=[]
    for case in c['cases']:
        cid=case['case_id'];row={'case_id':cid,'status':'unavailable'}
        try:
            prior=next(x for x in m['cases'] if x['case_id']==cid)['prior_source'];old=read_json(verify(prior))
            src=next(x for x in old['cases'] if x['case_id']==cid)['union'];core=src['original_core']
            state,raw,protein,adjacency,ga=source_state(core['parent']['path'],str(verify(src['source']['source_structure'])),str(verify(topology)),json.dumps(core,sort_keys=True))
            g=state['graph'];prep=read_json(verify(src['representations']['context']['preparation']))
            if prep['source']!=state['source']:raise InvalidArtifact('context and source differ')
            mapping={atom_key(a['source']):a['qm_index'] for a in prep['mapping']['source_to_qm'] if 'source' in a and a['source']['element'] not in ('H','D')}
            pqq={k for k in g.atoms if g.meta[k]['canonical_resname']=='PQQ' and g.meta[k]['element'] not in ('H','D')}
            if not pqq or not pqq<=mapping.keys():raise InvalidArtifact('complete PQQ map unavailable')
            query=sorted(set(mapping)&protein);outside=sorted(protein-set(mapping))
            if not query or not outside:raise InvalidArtifact('empty protein contact set')
            ts={z:next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(cid,z)) for z in ('Ca','La')}
            maps={z:read_json(verify(t['mapping']))['context'] for z,t in ts.items()}
            if maps['Ca']!=maps['La']:raise InvalidArtifact('paired maps differ')
            kin=Kinematics(maps['Ca']);kin_ids={tuple(v):i for i,v in enumerate(maps['Ca']['physical_ids']) if isinstance(v,list)}
            for key,index in mapping.items():
                if key not in kin_ids or [index,'source',kin_ids[key]] not in maps['Ca']['core_links']:raise InvalidArtifact('source/context link mismatch')
            points={};pins={}
            for name in ('origin','adaptive_Ca','adaptive_La','six_Ca','six_La'):
                pin=case['matrix']['Ca'][name]['xyz'];atoms=xyz(verify(pin))
                if name=='origin':q=np.zeros(len(kin.modes))
                else:
                    z=name.split('_')[1];task=ts[z]
                    r=read_json(verify(task['prior_receipt'])) if name.startswith('adaptive') else read_json(Path(lm['source_manifest']['path']).parent/'proposals'/task['task_id']/'result.json')
                    q=np.asarray(r['proposal']['full_q'])
                actual=kin.evaluate(q)[1]
                if not np.allclose(actual,[a[1:] for a in atoms],atol=1e-12,rtol=0):raise InvalidArtifact('actual coordinate map differs')
                for key in pqq:
                    if not np.array_equal(actual[mapping[key]],kin.core[mapping[key]]):raise InvalidArtifact('PQQ moved')
                for z in ('Ca','La'):
                    if not np.allclose(actual,[a[1:] for a in xyz(verify(case['matrix'][z][name]['xyz']))],atol=1e-12,rtol=0):raise InvalidArtifact('paired candidate differs')
                points[name]=actual[[mapping[k] for k in query]];pins[name]=pin
            rawq=np.asarray([[raw.atoms[k].pos.x,raw.atoms[k].pos.y,raw.atoms[k].pos.z] for k in query])
            if not np.allclose(points['origin'],rawq,atol=1e-8,rtol=0):raise InvalidArtifact('q0 differs from raw observed protein')
            exterior=np.asarray([[raw.atoms[k].pos.x,raw.atoms[k].pos.y,raw.atoms[k].pos.z] for k in outside])
            excluded=[within_three(k,adjacency,g) for k in query];keep=np.asarray([[k not in excluded[i] for k in outside] for i in range(len(query))])
            distances={name:np.linalg.norm(v[:,None,:]-exterior[None,:,:],axis=2) for name,v in points.items()}
            contacts={name:metrics(v,distances['origin'],keep,query,outside,g) for name,v in distances.items()}
            row.update(status='available',graph=ga,source_preparation=src['representations']['context']['preparation'],mapping=ts['Ca']['mapping'],coordinates=pins,
                       query_heavy_count=len(query),outside_heavy_count=len(outside),contacts=contacts)
        except Exception as exc:row['reason']=str(exc)
        rows.append(row)
    summary={'sources':9,'available':sum(r['status']=='available' for r in rows),
             'new_proposals_with_new_severe_contact':sum(r['contacts'][q]['new_below_2_0_A']>0 for r in rows if r['status']=='available' for q in ('six_Ca','six_La'))}
    result={'collection':record(collection),'topology':topology,'topology_authority':record(preparation),'policy':POLICY,'implementation':record(__file__),
            'shared_diagnostic':record(Path(__file__).with_name('omitted_context_contacts.py')),'rows':rows,'summary':summary,'new_molecular_calls':0,'scoring_filter_applied':False}
    write_new(output,result);return {'output':record(output),'summary':summary}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('collection','preparation','output'):p.add_argument('--'+k,required=True,type=Path)
    print(json.dumps(audit(**vars(p.parse_args())),indent=2))
