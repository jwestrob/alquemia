#!/usr/bin/env python3
"""Read-only joins for one frozen donor source tree; no build/native invocation."""
import hashlib,json,re,time
from pathlib import Path
ROOT=Path.cwd();D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1';W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1'
def pin(p):
 p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def verify(p):
 assert pin(p['path'])==p,p['path'];return Path(p['path'])
def routine(text,name):
 m=re.findall(r'^      subroutine '+re.escape(name)+r'\b.*?^      end\s*$',text,re.M|re.S);assert len(m)==1,name;return m[0]
def main():
 tick,cpu=time.monotonic(),time.process_time();mp=W/'donor_extension_native_v1/source_v1/MANIFEST.json';assert pin(mp)['sha256']=='45d45d65adc28d041298f8fe7e90b45e9a8425e7d4ff3f5636bbbc6f589e125a'
 m=json.loads(mp.read_text());base=json.loads(verify(m['baseline_source_manifest']).read_text());assert m['baseline_source_manifest']['sha256']=='82886b766d192b1fb0b89012dcd2bebe6315a2dcc83979b982441da0d464a3b4'
 checks={}
 def check(k,v):checks[k]=bool(v);assert v,k
 for label,v in [('candidate',m),('baseline',base)]:
  for p in v['source_outputs']:verify(p)
  verify(v['cmake'])
  check(label+'_615_source_and_CMake_pins',len(v['source_outputs'])==615)
 for k in ['design_sources','design','declaration','helper','implementation','patch','donor_inventory','water_native_actual_qualification']:verify(m[k])
 root=Path(m['source_root']);prior=Path(base['source_root']);new={str(Path(x['path']).relative_to(root)):x['sha256'] for x in m['source_outputs']};old={str(Path(x['path']).relative_to(prior)):x['sha256'] for x in base['source_outputs']}
 check('same_exact_source_inventory',set(new)==set(old));changed=sorted(k for k in new if new[k]!=old[k]);check('only_declared_three_files',changed==['source/induce.f','source/polar.f','source/testgrad.f'] and m['changed_files']==['induce.f','polar.f','testgrad.f'])
 before=(prior/'source/induce.f').read_text();after=(root/'source/induce.f').read_text();helper=verify(m['helper']).read_text()
 check('literal_helper_appended',after.endswith(helper))
 for name in ['dfield0a','ufield0a','induce0a','exportresponsematrix','localresponsegeometry','localresponsespectrum']:
  check('unchanged_'+name,routine(before,name)==routine(after,name))
 after_setup=routine(after,'localresponsesetup').replace('      call donorresponsesetup\n      if (lrdonor) return\n','')
 check('water_setup_only_early_dispatch',after_setup==routine(before,'localresponsesetup'))
 after_grad=routine(after,'localresponsegradient').replace('      if (lrdonor) then\n         call donorresponsegradient\n         return\n      end if\n','')
 check('water_gradient_only_donor_dispatch',after_grad==routine(before,'localresponsegradient'))
 after_report=routine(after,'localresponsereport').replace('      if (lrdonor) call donorresponsepermanent\n','')
 check('water_report_only_readonly_donor_moment',after_report==routine(before,'localresponsereport'))
 for name in ['moments.f','rotpole.f','torque.f','epolar.f','epolar1.f','epolar3.f','exfield.f','exfield1.f','exfield3.f','damping.f','kpolar.f','kmpole.f']:
  check('unchanged_'+name,new['source/'+name]==old['source/'+name])
 report=routine(after,'donorresponsepermanent');check('fixed_lab_permanent_formula', 'mu=mu+rpole(1,i)*(/x(i),y(i),z(i)/)+rpole(2:4,i)' in report and all(x not in report for x in ['xcm','xmid','uind','uinp']))
 top=routine(after,'induce');mom=(root/'source/moments.f').read_text();check('report_precedes_moments_induced_mutation',top.index('call induce0a')<top.index('call localresponsereport')<top.index('call exportresponsematrix') and mom.index('call rotpole')<mom.index('call induce')<mom.index('rpole(2:4,k) = rpole(2:4,k)'))
 check('no_intrinsic_alpha_damping_assignment',not re.search(r'^\s*(?:polarity|pdamp|thlval|rpole)\([^\n]*\)\s*=',helper,re.M|re.I))
 tg=(root/'source/testgrad.f').read_text();check('testgrad_reporting_after_full_gradient','call gradient (eval,derivs)\n            if (lrdonor)\n     &         call donorresponsetotalgradient (eval,derivs)' in tg)
 sources=json.loads(verify(m['design_sources']).read_text());inputs={}
 for metal,vals in sources['exact_central_inputs'].items():
  for p in vals.values():verify(p)
  inputs[metal]=vals
 check('exact_selected17621_inputs_unchanged',inputs['La']['key']['sha256']=='ec7e7644150e9382aa57c77f3341ee8d7b904b797b89828365e3b64ff15f9b08' and inputs['Dy']['key']['sha256']=='d1b823cbdc6ff210ea2b1ed3b21f8c24e5ba148740016b1be30a89104bd799ca' and inputs['La']['parameters']['sha256']==inputs['Dy']['parameters']['sha256']=='e1cc2bcc6d7fb319aba7abaaf6b4483f2ef5ba730727212b5c157a18a290f38b')
 check('source_scope_zero_build_native',m['new_native_calls']==m['builds']==0)
 out=W/'donor_extension_native_v1/independent_review_v1';out.mkdir(parents=True,exist_ok=False)
 r=dict(status='passed_frozen_donor_source_joins',implementation=pin(__file__),source_manifest=pin(mp),baseline_manifest=m['baseline_source_manifest'],checks=checks,check_count=len(checks),source_files_verified=1230,changed_files=changed,selected17621_inputs=inputs,new_native_calls=0,builds=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
 p=out/'SOURCE_CHECKS.json';p.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(dict(result=pin(p),check_count=len(checks),cost=r['cost']),indent=2))
if __name__=='__main__':main()
