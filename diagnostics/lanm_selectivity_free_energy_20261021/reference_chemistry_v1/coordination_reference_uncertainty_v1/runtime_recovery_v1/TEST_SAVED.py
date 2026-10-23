#!/usr/bin/env python3
"""Narrow versioned producer joins and staged release tests; no native calls."""
import argparse,copy,importlib.util,json,time
from pathlib import Path
s=importlib.util.spec_from_file_location('runtime_recovery',Path(__file__).with_name('REFERENCE.py'));m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--recovery-inventory',required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();tick,cpu=time.monotonic(),time.process_time();p=m.load(a.panel);checks=[]
    def check(x,name):m.b.check(x,name);checks.append(name)
    def reject(fn,name):
        try:fn()
        except (ValueError,KeyError,TypeError):checks.append(name);return
        raise ValueError('Expected rejection: '+name)
    a.out.mkdir(parents=True,exist_ok=False)
    for c in p['cells']:check(m.deck(p,c,4,3000)==m.b.verify(c['template_input']).read_text(),'unchanged deck '+c['cell_id'])
    sets=m.subsets(p);check([len(sets[x]) for x in ('controls_remaining','cluster_scout','cluster_remaining')]==[14,1,7],'exact14/1/7 counts')
    flat=sum(sets.values(),[]);check(len(set(flat))==22 and not set(flat)&set(m.old.CONTROL_SCOUTS),'no original scout or duplicate in continuation')
    inv=json.loads(Path(a.recovery_inventory).read_text());actual={}
    for cid,pin in inv.items():
        c=next(x for x in p['cells'] if x['cell_id']==cid);q=m.admit_recovered(m.b.verify(pin),c,p,a.panel);actual[cid]=q
        check(q['status']=='complete' and q['native_producer']==p['adapter'] and q['native_result']==q['original_result'],'actual preserved native producer '+cid)
    # Mock metadata only for prospective four-scout and cluster gates. Nothing
    # here is an actual fourth scout, qualification, cost approval or release.
    panel=m.b.pin(a.panel);source=m.b.pin(m.__file__);parser=m.b.pin(m.rec.__file__)
    original={cid:{'synthetic_original':cid} for cid in m.old.CONTROL_SCOUTS};recoveries={cid:{'synthetic_recovery':cid} for cid in m.old.CONTROL_SCOUTS}
    records={
        'prep':dict(status='passed_static_preflight',panel=panel),
        'runtime':dict(status='passed_runtime_handoff_preparation',implementation=source,parser=parser,panel=panel,all_26_original_decks_identical=True,new_native_calls=0),
        'review':dict(status='passed_runtime_handoff_review',implementation=source,parser=parser,panel=panel),
        'parser_review':dict(status='passed_saved_parser_recovery_review',implementation=parser,panel=panel),
        'scouts':dict(status='passed_control_scout_review',panel=panel,parser=parser,parser_recovery_review='parser_review',results=original,recoveries=recoveries),
        'qualification':dict(status='local_controls_passed',panel=panel,complete=18,declared=18,all_controls_passed=True,checks={str(i):True for i in range(35)},implementation=source),
        'cluster':dict(status='passed_cluster_scout_review',panel=panel,cost_and_memory_accepted=True,result='synthetic_cluster_result')}
    release=dict(status='approved_coordination_local_cc_subset',implementation=source,native_input_adapter=p['adapter'],parser=parser,panel=panel,
        cell_ids=[c['cell_id'] for c in p['cells']],maximum_native_calls=26,ranks=4,preparation='prep',runtime_preparation='runtime',independent_review='review',
        parser_recovery_review='parser_review',stage='controls_remaining',released_cell_ids=sets['controls_remaining'],control_scout_review='scouts',
        control_qualification='qualification',cluster_scout_review='cluster')
    oldread,oldverify,oldrec,oldadmit=m.read,m.b.verify,m.admit_recovered,m.admit_result
    m.read=lambda pin:copy.deepcopy(records[pin])
    m.b.verify=lambda pin:pin
    m.admit_recovered=lambda pin,c,p,path:dict(original_result=original[c['cell_id']])
    m.admit_result=lambda path,c,p,panel:dict(status='complete')
    try:
        check(m.release_contract(p,a.panel,release)==release,'mock complete four-scout gate admits exact14')
        for field,value,name in [('maximum_native_calls',27,'reject extra native budget'),('parser',{},'reject unreviewed parser'),('stage','controls_scout','reject repeated original scouts'),('released_cell_ids',sets['controls_remaining'][:-1],'reject changed remaining subset')]:
            bad=copy.deepcopy(release);bad[field]=value;reject(lambda:m.release_contract(p,a.panel,bad),name)
        saved=records['scouts']['recoveries'];records['scouts']['recoveries']={k:v for k,v in saved.items() if k!=m.old.CONTROL_SCOUTS[-1]}
        reject(lambda:m.release_contract(p,a.panel,release),'missing fourth recovery cannot release14');records['scouts']['recoveries']=saved
        records['scouts']['results']=dict(original);records['scouts']['results'][m.old.CONTROL_SCOUTS[0]]={'wrong':True}
        reject(lambda:m.release_contract(p,a.panel,release),'wrong original producer pin rejected');records['scouts']['results']=original
        scout=dict(release,stage='cluster_scout',released_cell_ids=sets['cluster_scout']);check(m.release_contract(p,a.panel,scout)==scout,'mock complete18 admits only biggest cluster')
        records['qualification']['checks']['0']=False;reject(lambda:m.release_contract(p,a.panel,scout),'failed approximation gate blocks cluster');records['qualification']['checks']['0']=True
        records['qualification']['complete']=17;reject(lambda:m.release_contract(p,a.panel,scout),'partial17 blocks cluster');records['qualification']['complete']=18
        remaining=dict(release,stage='cluster_remaining',released_cell_ids=sets['cluster_remaining']);check(m.release_contract(p,a.panel,remaining)==remaining,'mock admitted largest scout permits only7')
        records['cluster']['cost_and_memory_accepted']=False;reject(lambda:m.release_contract(p,a.panel,remaining),'unaccepted actual scout cost blocks7')
    finally:m.read,m.b.verify,m.admit_recovered,m.admit_result=oldread,oldverify,oldrec,oldadmit
    # The changed collector must still retain the full denominator on absence.
    missing={c['cell_id']:None for c in p['cells']};m.b.write(a.out/'MISSING26.json',missing)
    q=m.collect(a.panel,a.out/'MISSING26.json',a.out/'missing26','full')
    check(q['complete']==0 and q['declared']==len(q['rows'])==26 and not q['all_controls_passed'],'missing full26 remains unavailable')
    check(m.control_score is m.old.control_score,'unchanged35-control arithmetic reused directly')
    q=dict(status='passed_runtime_handoff_controls',implementation=source,parser=parser,test=m.b.pin(__file__),panel=panel,recovery_inventory=m.b.pin(a.recovery_inventory),
        actual_recovery_pins=inv,actual_recovery_count=len(actual),prospective_gates_mock_only=True,checks=checks,count=len(checks),new_native_calls=0,local_threads=1,
        wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu)
    pin=m.b.write(a.out/'RESULT.json',q);print(json.dumps(dict(result=pin,status=q['status'],count=len(checks))))
if __name__=='__main__':main()
