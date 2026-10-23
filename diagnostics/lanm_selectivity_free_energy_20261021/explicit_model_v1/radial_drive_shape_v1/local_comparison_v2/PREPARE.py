#!/usr/bin/env python3
"""Metadata-only recovery preparation; retain admitted science and failed v1."""
import argparse,ast,difflib,json,time
from pathlib import Path
import RUN as R

def main(config,out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    c=json.loads(Path(config).read_text());old=R.HERE.parent/'local_comparison_v1';ready=json.loads((old/'READY.json').read_text());previous=R.read(ready['configuration_template'])
    delta={'implementation','adapter','math','solver','specification','run_id','output_directory','recovery_authority'}
    R.require(set(c)==set(previous)|{'recovery_authority'},'Only recovery provenance added')
    R.require(all(c[k]==v for k,v in previous.items() if k not in delta),'Scientific configuration exact')
    for k,name in [('implementation','RUN.py'),('adapter','ADAPTER.py'),('math','MATH.py'),('solver','SOLVE.py'),('specification','SPEC.md')]:R.require(c[k]==R.pin(R.HERE/name),'Exact final new source')
    R.require(c['run_id']=='radial_local_comparison_v2' and c['output_directory']==previous['output_directory'].replace('/local_comparison_v1/','/local_comparison_v2/') and not Path(c['output_directory']).exists(),'Fresh exact recovery output')
    checks=[]
    for name in ('MATH.py','SOLVE.py','SPEC.md'):
        R.require((old/name).read_bytes()==(R.HERE/name).read_bytes(),'Unchanged '+name);checks.append(name+' byte identical')
    before=(old/'ADAPTER.py').read_text();after=(R.HERE/'ADAPTER.py').read_text()
    R.require(after==before.replace('self.ledger.update(new_policy=','self.ledger.update(limits=copy.deepcopy(cfg[\'limits\']),new_policy='),'Sole adapter limits delta');checks.append('adapter sole detached-limits insertion')
    a,b=ast.parse((old/'RUN.py').read_text()),ast.parse((R.HERE/'RUN.py').read_text())
    funcs=lambda t:{n.name:n for n in t.body if isinstance(n,ast.FunctionDef)}
    aa,bb=funcs(a),funcs(b)
    R.require(set(bb)==set(aa)|{'persist_result'},'Sole new reporting helper')
    for name,node in aa.items():
        if name!='run':R.require(ast.dump(node)==ast.dump(bb[name]),'Unchanged runtime authority '+name)
    original_try=next(n for n in aa['run'].body if isinstance(n,ast.Try));new_try=next(n for n in bb['run'].body if isinstance(n,ast.Try))
    R.require(ast.dump(ast.Module(body=original_try.body,type_ignores=[]))==ast.dump(ast.Module(body=new_try.body,type_ignores=[])),'Scientific run body unchanged');checks.append('scientific run try body and all authority checks AST identical')
    authority=c['recovery_authority'];[R.verify(p) for p in authority.values()]
    prep=R.read(ready['preparation']);R.require(prep['status']=='prepared_radial_local_comparison_zero_calls' and prep['native_programs']==prep['solver_calls']==0,'Inherited completed zero-call preparation')
    ctrl=R.read(authority['runtime_controls']);R.require(ctrl['status']=='passed_radial_runtime_recovery_controls' and ctrl['run']==c['implementation'] and ctrl['adapter']==c['adapter'] and ctrl['checks_passed']==14 and ctrl['native_programs']==ctrl['solver_calls']==0,'Actual bounded recovery controls')
    for p in [ready['inputs'],ready['inventory'],ready['controls'],c['root_cache_bridge'],c['qualification_review'],c['qualification_root_review']]:R.verify(p)
    diff=''.join(''.join(difflib.unified_diff((old/n).read_text().splitlines(True),(R.HERE/n).read_text().splitlines(True),fromfile=str(old/n),tofile=str(R.HERE/n))) for n in ('ADAPTER.py','RUN.py'))
    (out/'SOURCE_DIFF.patch').write_text(diff)
    result={'status':'prepared_radial_local_comparison_zero_calls','preparation_implementation':R.pin(__file__),'configuration':R.pin(config),'configuration_fingerprint':R.fingerprint(c),'implementation':c['implementation'],'adapter':c['adapter'],'math':c['math'],'solver':c['solver'],'source_diff':R.pin(out/'SOURCE_DIFF.patch'),'checks':checks,'inherited_scientific_preparation':ready['preparation'],'inventory':ready['inventory'],'inherited_mathematical_controls':ready['controls'],'runtime_controls':authority['runtime_controls'],'recovery_authority':authority,'native_programs':0,'solver_calls':0,'scientific_vector_replay':False,'actual_Jacobian_constructed':False,'limits':R.LIMITS,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    R.save(out/'RESULT.json',result);print(json.dumps({'result':R.pin(out/'RESULT.json'),'status':result['status'],'cost':result['cost']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();main(a.config,a.out)
