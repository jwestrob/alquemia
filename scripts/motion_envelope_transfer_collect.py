"""Read-only recovery of the missing scalar-to-pool cell join metadata."""
import argparse,json
from pathlib import Path
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
import consistent_context as context
import motion_envelope_transfer as runner


def collect(manifest,output):
    runner.validate_pool(manifest);m=read_json(manifest);root=Path(manifest).parent
    cells={(t['case_id'],t['metal'],t['candidate']):t for t in m['tasks']}
    if len(cells)!=len(m['tasks']):raise InvalidArtifact('ambiguous pool cell identity')
    manifests={};joins=[]
    for pin in read_json(root/'solvent/INDEX.json')['shards']:
        path=verify(pin);lm=read_json(path)
        if lm['source']!=record(manifest):raise InvalidArtifact('scalar manifest is not from this exact pool')
        for t in lm['tasks']:
            target=cells[t['case_id'],t['metal'],t['candidate']]
            if not context.reusable_state(t,target):raise InvalidArtifact('physical scalar/cell identity differs')
            if t.get('cell_id',target['task_id'])!=target['task_id']:raise InvalidArtifact('existing explicit cell mapping differs')
            t['cell_id']=target['task_id'];joins.append({'scalar_task_id':t['task_id'],'cell_id':target['task_id'],'medium':t['medium'],'manifest':pin})
        manifests[str(path.resolve())]=lm
    _,p=runner.engine()
    original=p.read_json
    def joined_read(path):return manifests.get(str(Path(path).resolve()))or original(path)
    p.read_json=joined_read
    join_path=Path(output).with_name(Path(output).stem+'_JOIN.json')
    write_new(join_path,{'pool_manifest':record(manifest),'implementation':record(__file__),'joins':joins,
        'rule':'exact case/metal/candidate identity plus coordinate/charge/multiplicity; no energy alteration',
        'new_molecular_calls':0,'original_manifests_modified':False})
    return p.collect(manifest,output)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True);p.add_argument('--output',required=True)
    print(json.dumps(collect(**vars(p.parse_args())),indent=2))
