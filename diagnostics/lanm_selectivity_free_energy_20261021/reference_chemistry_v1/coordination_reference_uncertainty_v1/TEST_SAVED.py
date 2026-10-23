#!/usr/bin/env python3
"""Narrow saved parser, arithmetic, missing-denominator and release controls. No engines."""
import argparse,copy,importlib.util,json,math,time
from pathlib import Path
s=importlib.util.spec_from_file_location('test_local_cc',Path(__file__).with_name('REFERENCE.py'));m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();tick,cpu=time.monotonic(),time.process_time();p=m.load(a.panel);checks=[]
    def check(ok,name):m.b.check(ok,name);checks.append(name)
    def reject(fun,name):
        try:fun()
        except (ValueError,KeyError,TypeError):checks.append(name);return
        raise ValueError('Expected rejection: '+name)
    a.out.mkdir(parents=True,exist_ok=False)
    check(len(p['cells'])==26 and sum(c['stage']=='controls' for c in p['cells'])==18,'26/18 denominator')
    for c in p['cells']:
        text=m.deck(p,c,4,3000);check(text==m.b.verify(c['template_input']).read_text(),'exact prepared input '+c['cell_id'])
    cc='Iter E(tot) E(Corr) Delta-E Residual Time\n 1 -1.0 -.2 -.00001 5e-9 .2\n--- The CCSD iterations have converged ---\n'
    trace=m.iteration_rows(cc,r'Iter\s+E\(tot\)[^\n]*','cc')[0];check(trace['terminal_value']==5e-9 and trace['recognized'],'CC header explicit residual')
    split=cc.replace(' 1 -1.0',' -1.0');check(m.iteration_rows(split,r'Iter\s+E\(tot\)[^\n]*','cc')[0]['terminal_value']==5e-9,'PrintLevel3 separate iteration label')
    bad=cc.replace('Residual Time','Time Residual');check(not m.iteration_rows(bad,r'Iter\s+E\(tot\)[^\n]*','cc')[0]['recognized'],'reordered residual header rejected')
    bad=cc.replace('5e-9','********');check(m.iteration_rows(bad,r'Iter\s+E\(tot\)[^\n]*','cc')[0]['unsupported_numeric_row'],'starred CC residual rejected')
    tri='Iter E(Triple) Delta-E Time\n 3 -.012 7e-9 .2\nTRIPLES CORRECTION (T1)\n';tr=m.iteration_rows(tri,r'Iter\s+E\(Triple\)[^\n]*','triples')[0]
    check(tr['terminal_value']==7e-9 and tr['recognized'],'triples named Delta-E column')
    # Synthetic additional DLPNO metadata on an actual admitted canonical output:
    # validates parser mechanics only, never a molecular result or method qualification.
    c=next(c for c in p['cells'] if c['cell_id']=='control__water__Dy__complex');saved=p['controls'][c['canonical_id']];raw=m.b.verify(saved['output']).read_text()
    extra='\nLOCALIZATION HAS CONVERGED.\n LMP2 iterations CONVERGED\nParameters used in iterative triples\nETol ... 1.000e-08\n'+tri
    extra+='\n'.join(k+' = '+str(v) for k,v in m.TIERS['tight'].items())+'\n'
    extra+='AUXILIARY/C BASIS SET INFORMATION\n---------------------------------\n Group 1 Type Dy : 369s contracted to 369s\n Group 2 Type O : 10s contracted to 10s\n Group 3 Type H : 5s contracted to 5s\n'
    extra+=''.join(f'Atom {i}{s} basis set group => {g}\n' for i,s,g in [(0,'Dy',1),(1,'O',2),(2,'H',3),(3,'H',3)])
    extra+='---------------------------------\n# of basis functions in Aux-C ... 389\nNumber of AuxC basis functions ... 389\n'
    parsed=m.parse(raw+extra,c,0);check(parsed['status']=='complete','synthetic DLPNO parser positives on saved inventory')
    for old,new,key in [('LOCALIZATION HAS CONVERGED.','ORBITAL LOCALIZATION','localization_evidence'),('ETol ... 1.000e-08','ETol ... 1.000e-04','iterative_triples'),('TCutPNO = 1e-07','TCutPNO = 1e-06','local_TCutPNO'),('Group 1 Type Dy : 369s contracted to 369s','Group 1 Type Dy : 368s contracted to 368s','exact_metal_aux_C')]:
        check(not m.parse(raw+extra.replace(old,new),c,0)['checks'][key],'reject '+key)
    check(m.parse(raw+extra,c,1)['status']=='failed_admission','exit nonzero rejection')
    check(not m.parse(raw+extra+'\nFINAL SINGLE POINT ENERGY ********\n',c,0)['checks']['complete_energy_tokens'],'overflow final token rejected')
    model={}
    for c in p['cells']:
        if c['stage']!='controls':continue
        r=p['controls'][c['canonical_id']];model[c['cell_id']]=dict(energy_hartree=r['energy_hartree'],RHF_energies_hartree=[r['RHF_energy_hartree']],aux_C={'synthetic_equal':True})
    score=m.control_score(p,model);check(score['all_controls_passed'] and len(score['checks'])==35,'canonical self-control exact full35 checks')
    changed=copy.deepcopy(model);changed['control__water__Dy__complex']['energy_hartree']+=.5/m.b.HA;sc=m.control_score(p,changed);check(not sc['checks']['CP__water__Dy'] and not sc['checks']['paired_CP__water'],'half-kcal CP perturbation rejected')
    check(abs(sc['paired_exchange_local_canonical']['acetamide'][0]-sc['paired_exchange_local_canonical']['acetamide'][1]+.5)<1e-9,'water-to-donor sign')
    check(set(m.CONTROL_SCOUTS)<={c['cell_id'] for c in p['cells']} and len(m.CONTROL_SCOUTS)==4,'exact root four scouts')
    reject(lambda:m.release_contract(p,a.panel,dict(status='unreleased')),'unreleased native gate rejected')
    inv={c['cell_id']:None for c in p['cells']};m.b.write(a.out/'MISSING_INVENTORY.json',inv)
    result=m.collect(a.panel,a.out/'MISSING_INVENTORY.json',a.out/'missing26','full');check(result['complete']==0 and result['declared']==len(result['rows'])==26 and result['status']=='unavailable' and not result['all_controls_passed'],'missing full26 cannot pass')
    receipt=dict(status='passed_saved_controls',implementation=m.b.pin(m.__file__),test=m.b.pin(__file__),panel=m.b.pin(a.panel),checks=checks,count=len(checks),new_native_calls=0,local_threads=1,wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu,synthetic_parser_evidence_is_not_native_admission=True)
    m.b.write(a.out/'RESULT.json',receipt);print(json.dumps(dict(status=receipt['status'],checks=len(checks))))
if __name__=='__main__':main()
