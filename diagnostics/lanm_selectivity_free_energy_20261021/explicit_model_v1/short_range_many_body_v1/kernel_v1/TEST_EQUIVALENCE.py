"""Additional chemical-exchange checks; no repeat of the 24-check suite."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import traceback
import TEST_KERNEL as t
from dataclasses import replace


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);args=parser.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();cpu=time.process_time();rows=[]
    x,g=t.fixture();p=t.parameters();sites=tuple(range(len(x)))
    baseline=t.k.correction(x,g,sites,p)
    def compare(actual,perm=None):
        ix=list(range(len(baseline['h']))) if perm is None else [3*i+j for i in perm for j in range(3)]
        return {'V':t.close(actual['V'],baseline['V']),
                'h':t.close(actual['h'],baseline['h'][ix]),
                'H':t.close(actual['H'],baseline['H'][ix][:,ix])}
    def add(name,call):
        try:rows.append({'name':name,'passed':True,'details':call()})
        except Exception as exc:rows.append({'name':name,'passed':False,'error':repr(exc),'traceback':traceback.format_exc()})
    def swap(a,b):
        order=list(range(len(x)));order[a],order[b]=order[b],order[a]
        # Metadata deliberately stays FIXED: this is physical exchange of
        # equivalent atom positions, not coordinated relabeling of frame lists.
        return compare(t.k.correction(x[order],g,sites,p),order)
    add('physical_water_H_exchange_fixed_metadata',lambda:swap(2,3))
    add('physical_acetate_O_exchange_fixed_metadata',lambda:swap(4,5))
    def reverse_frame():
        water=replace(g.ligands[0],frame_atoms_by_donor={1:(3,2)})
        return compare(t.k.correction(x,replace(g,ligands=(water,)+g.ligands[1:]),sites,p))
    add('equivalent_H_frame_order_no_privileged_atom',reverse_frame)
    add('intact_ligand_order_all_response_channels',lambda:compare(t.k.correction(x,replace(g,ligands=tuple(reversed(g.ligands))),sites,p)))
    result={'status':'passed' if all(r['passed'] for r in rows) else 'failed',
            'scope':'four additional manufactured chemical symmetry checks only',
            'checks':rows,'native_calls':0,'quantum_calls':0,'fits':0,'cached_operator_replays':0,
            'cost':{'test_wall_seconds':time.monotonic()-start,'test_cpu_seconds':time.process_time()-cpu,
                    'process_cpu_seconds_including_imports':time.process_time()},
            'runtime':{'torch_threads':t.torch.get_num_threads(),'interop_threads':t.torch.get_num_interop_threads()},
            'sources':{name:{'path':str(Path(__file__).with_name(name).resolve()),'sha256':hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()} for name in ('KERNEL.py','TEST_KERNEL.py','TEST_EQUIVALENCE.py','SPEC.md')}}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'passed':sum(r['passed'] for r in rows),'total':len(rows),'cost':result['cost']}))
    raise SystemExit(0 if result['status']=='passed' else 1)


if __name__=='__main__':main()
