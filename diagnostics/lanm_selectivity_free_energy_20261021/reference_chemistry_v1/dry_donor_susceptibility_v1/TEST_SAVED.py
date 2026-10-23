#!/usr/bin/env python3
"""Small saved controls; transformed fixtures are explicitly synthetic, never outcomes."""
import argparse,copy,importlib.util,json,math,re,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('dryfield',HERE/'REFERENCE.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def main():
 a=argparse.ArgumentParser();a.add_argument('--panel',required=True);a.add_argument('--out',type=Path,required=True);args=a.parse_args();args.out=args.out.resolve();tick,cpu=time.monotonic(),time.process_time();args.out.mkdir(parents=True,exist_ok=False)
 checks={}
 def check(name,condition):
  checks[name]=bool(condition);m.b.check(condition,name)
 def rejects(name,fn):
  try:fn()
  except Exception as exc:checks[name]={'passed':True,'error':repr(exc)};return
  raise AssertionError('Expected rejection: '+name)
 p=m.load(args.panel);check('actual_ten_prepared',m.validate(args.panel)['declared']==10)
 rows={}
 for metal in m.METALS:
  c=next(c for c in p['cells'] if c['metal']==metal and c['field_tag']=='zero');oc=p['centers'][metal]
  text=m.b.verify(oc['zero_EnGrad_output']).read_text();raw=m.b.verify(oc['zero_EnGrad_property']).read_text();q=p['old_EnGrad_context'][metal]['native']
  # Synthetic origin-zero fixture: same raw central energies and native metadata, explicitly transformed moment tokens.
  geom=m.prop.section(raw,'Geometry').split('&Units "Bohr"]',1)[1]
  rr=re.findall(r'^\s*([A-Za-z]+)(:?)\s+('+m.moment.NUMBER+r')\s+('+m.moment.NUMBER+r')\s+('+m.moment.NUMBER+r')\s*$',geom,re.M)
  pos=[[float(v) for v in r[2:]] for r in rr];nuc=[math.fsum(z*r[j] for z,r in zip(c['expected_effective_nuclear_charges'],pos)) for j in range(3)]
  total=oc['old_total_dipole_lab_au'];elec=[a-z for a,z in zip(total,nuc)]
  for name,label,values in [('dipoleTotal','Total Dipole Moment',total),('dipoleNucContrib','Nuclear contribution',nuc),('dipoleElecContrib','Electronic contribution',elec)]:
   before,after=raw.split('&'+name+' ',1);part,*rest=after.split('\n   &',1)
   part=re.sub(r'^(\s*)([012])\s+('+m.moment.NUMBER+r')\s*$',lambda hit:hit[1]+hit[2]+' '+format(values[int(hit[2])],'.17e'),part,flags=re.M)
   raw=before+'&'+name+' '+part+('\n   &'+rest[0] if rest else '')
   text=re.sub(r'^\s*'+re.escape(label)+r'\s*:.*$',label+' : '+' '.join(f'{v:.9f}' for v in values),text,flags=re.M)
  text=re.sub(r'^Position of electric origin\s*\.\.\..*$', 'Position of electric origin ... 0.000000 0.000000 0.000000',text,flags=re.M)
  folder=args.out/(metal+'_synthetic_lab0');folder.mkdir();op=folder/'synthetic.out';pp=folder/'synthetic.property.txt';op.write_text(text);pp.write_text(raw)
  props=m.field_properties(op,pp,c,q);check(metal+'_D4_on_property',abs(props['energy_hartree']-props['final_energy_hartree'])>.001)
  check(metal+'_raw_energy_unchanged',props['energy_hartree']==q['electronic_energy_hartree'])
  bad=raw.replace('&Charge [&Type "Integer"] 2','&Charge [&Type "Integer"] 3');bp=folder/'wrong_charge.property.txt';bp.write_text(bad);rejects(metal+'_wrong_charge',lambda:m.field_properties(op,bp,c,q))
  bp=folder/'wrong_final.property.txt';bp.write_text(re.sub(r'(&FinalEnergy\s+\[&Type "Double"\]\s+)('+m.moment.NUMBER+r')',lambda h:h[1]+str(float(h[2])+.001),raw));rejects(metal+'_wrong_final',lambda:m.field_properties(op,bp,c,q))
  bp=folder/'nonfinite.property.txt';bp.write_text(raw.replace(props['raw_tokens']['electronic'],'nan'));rejects(metal+'_nonfinite',lambda:m.field_properties(op,bp,c,q))
  bp=folder/'wrong_nuclear.property.txt';bp.write_text(raw.replace(props['raw_tokens']['dipoles']['nuclear'][0],str(nuc[0]+1)));rejects(metal+'_wrong_nuclear_origin',lambda:m.field_properties(op,bp,c,q))
  # Echo/parser control only; remove analytic-stage marker in a clearly synthetic output, never relabel an actual EnGrad result.
  echo='\n'.join('| %3d> %s'%(i,l) for i,l in enumerate(m.deck(p,c,4,3000).splitlines(),1))+'\n'
  head,tail=text.split('****END OF INPUT****',1);tt=re.sub(r'^\|\s*\d+>.*(?:\n|$)','',head,flags=re.M)+echo+'****END OF INPUT****'+tail;tt=tt.replace('ORCA SCF GRADIENT CALCULATION','SYNTHETIC REMOVED GRADIENT MARKER')
  parsed=m.parse(tt,c,0);check(metal+'_field_echo_and_native_parser',parsed['status']=='complete')
  check(metal+'_wrong_field_echo_reject',m.parse(tt.replace(m.field_line(c),' EField 0,0,1',1),c,0)['status']=='failed_admission')
  check(metal+'_native_failure_reject',m.parse(tt,c,1)['status']=='failed_admission')
  # Analytic quadratic synthetic field values, known susceptibility and conjugate moment.
  chi=25. if metal=='La' else 27.;mu=.3 if metal=='La' else .4;t=p['field_direction_lab']
  for cell in [v for v in p['cells'] if v['metal']==metal]:
   f=cell['field_au'];pr=copy.deepcopy(props);pr['energy_hartree']=-1-mu*f-.5*chi*f*f;pr['final_energy_hartree']=pr['energy_hartree']+props['D4_energy_hartree'];pr['dipole_au']=[(mu+chi*f)*a for a in t]
   rows[cell['cell_id']]={'status':'complete','properties':pr}
 result=m.score(rows,p,'full');check('known_quadratic_all16',result['numerical_checks_passed'] and len(result['checks'])==16)
 check('paired_sign',abs(result['Dy_minus_La']['estimates']['h']['chi_energy']-2)<1e-8)
 check('projection_not_lab_z',abs(result['per_metal']['La']['estimates']['h']['zero_moment']-.3)<1e-14)
 check('six_check_scout',len(m.score(rows,p,'scout')['checks'])==6)
 missing=copy.deepcopy(rows);missing.pop(p['cells'][0]['cell_id']);v=m.score(missing,p,'full');check('missing_zero_retains16_and_fails',len(v['checks'])==16 and not v['passed'] and v['per_metal']['La'] is None)
 bad=copy.deepcopy(rows);bad[p['cells'][1]['cell_id']]['properties']['energy_hartree']=float('nan');rejects('nonfinite_estimator',lambda:m.score(bad,p,'full'))
 bad=copy.deepcopy(rows);bad[p['cells'][1]['cell_id']]['properties']['D4_energy_hartree']+=1e-4;check('field_dependent_D4_fails',not m.score(bad,p,'full')['passed'])
 bad=copy.deepcopy(rows);bad[p['cells'][1]['cell_id']]['properties']['energy_hartree']+=1e-5;check('numerical_curvature_failure_retained',not m.score(bad,p,'full')['passed'])
 comp=m.compare_model(result,p);check('actual_model_center_join',set(comp['model_minus_quantum'])==set(m.METALS) and set(comp['paired_model_minus_quantum'])=={'h','half'})
 # Exact scoped release gate, actual preparation + explicitly synthetic review fixture.
 review=args.out/'synthetic_review.json';m.save(review,dict(status='passed_dry_donor_susceptibility_source_review',implementation=p['adapter'],panel=m.b.pin(args.panel)))
 release=dict(status='approved_dry_donor_susceptibility_reference',implementation=p['adapter'],panel=m.b.pin(args.panel),specification=p['specification'],cell_ids=[c['cell_id'] for c in p['cells']],maximum_native_calls=10,ranks=4,nodes=1,preparation=m.b.pin(Path(args.panel).with_name('PREFLIGHT.json')),independent_review=m.b.pin(review),stage='scout',released_cell_ids=p['scout_cell_ids'],execution_root=str(args.out/'synthetic_root'))
 check('La_five_release',m.release_contract(p,args.panel,release)==release)
 for key,value in [('released_cell_ids',p['remaining_cell_ids']),('maximum_native_calls',11),('nodes',2)]:
  bad=copy.deepcopy(release);bad[key]=value;rejects('release_'+key,lambda:m.release_contract(p,args.panel,bad))
 bad=copy.deepcopy(release);bad.update(stage='remaining',released_cell_ids=p['remaining_cell_ids']);rejects('Dy_requires_actual_La',lambda:m.release_contract(p,args.panel,bad))
 root=args.out/'synthetic_claims';ids=release['cell_ids']
 for cid in ids:m.claim(root,root/'tasks'/cid/'execution',cid,ids)
 rejects('ten_cap_duplicate',lambda:m.claim(root,root/'tasks'/ids[0]/'execution',ids[0],ids));check('ten_ledger',len(json.loads((root/'NATIVE_LEDGER.json').read_text()))==10)
 # No native calls; full missing denominator collector is an actual saved-only execution path.
 inv=args.out/'missing_inventory.json';m.save(inv,{c['cell_id']:str(args.out/'absent'/c['cell_id']/'RESULT.json') for c in p['cells']})
 check('missing_collector_nonzero',m.collect(args.panel,inv,args.out/'missing_collection','full')==1)
 absent=json.loads((args.out/'missing_collection/RESULT.json').read_text());check('missing_collector_all10_all16',absent['complete']==0 and len(absent['rows'])==10 and len(absent['response']['checks'])==16)
 out=dict(status='passed_saved_controls',implementation=m.b.pin(m.__file__),test_implementation=m.b.pin(__file__),panel=m.b.pin(args.panel),checks=checks,check_count=len(checks),new_native_calls=0,synthetic_fixtures_explicit=True,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
 m.save(args.out/'RESULT.json',out);print(json.dumps(dict(status=out['status'],checks=len(checks),cost=out['cost'])))
if __name__=='__main__':main()
