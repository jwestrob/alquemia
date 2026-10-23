"""Saved two-box operator review: no native work and no full stencil replay."""
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time


def pin(path):
    path=Path(path).resolve();h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return dict(path=str(path),sha256=h.hexdigest())


def main(result_path,output):
    start,cpu=time.monotonic(),time.process_time();checks=[]
    def check(name,value):
        if not value:raise ValueError(name)
        checks.append(name)
    def verify(p):
        check('pin '+p['path'],pin(p['path'])==p);return Path(p['path'])
    def read(p):return json.loads(verify(p).read_text())
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    report=dict(status='started',reviewer_implementation=pin(__file__),result=pin(result_path),checks=checks)
    try:
        check('one thread',all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')))
        result=read(report['result'])
        check('terminal success',result['status']=='full_box_operator_initialization_complete' and result['error'] is None)
        cfg=read(result['configuration']);release=read(cfg['root_release']);source_review=read(cfg['root_source_review'])
        check('released exact configuration',release['status']=='released_two_full_box_operator_initializations' and release['configuration_sources']=={k:v for k,v in cfg.items() if k!='root_release'})
        check('reviewed implementation',result['implementation']==cfg['implementation']==source_review['implementation'] and source_review['status']=='passed_operator_initialization_source_review')
        verify(result['implementation']);ready=read(cfg['ready']);inputs=read(ready['inputs'])
        check('reviewed preparation',source_review['ready']==cfg['ready']==result['ready'])
        author=read(source_review['independent_source_review']);verify(author['controls'])
        check('independent source binds implementation',author['implementation']==result['implementation'])
        composition=read(cfg['composition']);qualification=read(cfg['qualification']);finite=read(inputs['finite_export_root_actual_review'])
        check('finite actual qualification authority',finite['execution']==cfg['qualification'] and finite['status']=='admitted_finite_initialized_PME_export_operator_qualification')
        check('library/source identity',result['library']==qualification['library']==inputs['read_only_library_authority'] and result['source_manifest']==qualification['source_manifest']==inputs['read_only_source_authority'])
        verify(result['library']);verify(result['source_manifest']);verify(result['executable'])
        check('candidate/composition identity',result['candidate']==composition['candidate'] and result['composition']==cfg['composition'])
        candidate=composition['candidate'];actual=read(candidate['actual_review'])
        check('candidate actual admission',actual['status']=='passed_component_fit_actual_review' and actual['selected_native_admission_passed'] is True and actual['selected']==candidate['selected'] and actual['result']==candidate['result'])
        check('exact program counts',result['native_calls']==2 and result['link_calls']==1 and release['native_call_cap']==2 and release['link_cap']==1)
        check('no physical operations/admission',all(result[k]==0 for k in ('induced_solves','energy_calls','gradient_calls')) and result['full_stability_admitted'] is False and result['physical_accuracy_qualified'] is False)
        check('serial login resource identity',result['threads']==result['maximum_native_workers']==1 and result['execution_mode']=='login' and result['scheduler_allocated_CPU_seconds']==0 and result['host']==cfg['hostname'])
        deps=read(cfg['link_dependencies'])
        check('qualified link recipe',deps['qualified_execution']==cfg['qualification'] and deps['qualified_argv']==qualification['link']['argv'])
        for p in deps['files']+qualification['modules']:verify(p)
        argv=list(qualification['link']['argv']);argv[7]=ready['harness']['path'];argv[-1]=result['executable']['path']
        check('actual isolated link',result['link']['argv']==argv)
        spec=importlib.util.spec_from_file_location('bound_derive',verify(ready['derive']));derive=importlib.util.module_from_spec(spec);spec.loader.exec_module(derive)
        cells={};operation_cost={}
        for name in ('link','native_La','native_Dy'):
            rec=result[name];verify(rec['stdout']);verify(rec['stderr'])
            check(name+' completed',rec['status']=='complete' and rec['returncode']==0)
            check(name+' child limits',rec['wall_limit_seconds']==120 and rec['address_space_limit_bytes']==8*1024**3 and 0<=rec['wall_seconds']<=121 and rec['child_CPU_seconds']>=0)
            operation_cost[name]={k:rec[k] for k in ('wall_seconds','child_CPU_seconds','PID')}
        keys=('c1','cinf','cF','U_R_nominal','evenness_perturbation_nominal','q_min','q_max','q_abs_max','q_zero','q_conjugate_difference_max')
        for original in composition['files']:
            metal=original['metal'];cell=result['cells'][metal];folder=Path(cfg['output_directory'])/metal
            check(metal+' complete',cell['status']=='complete')
            check(metal+' argv/cwd',result['native_'+metal]['argv']==[result['executable']['path'],'system.xyz',metal] and result['native_'+metal]['cwd']==str(folder))
            for name,p in cell['inputs'].items():check(metal+' input location '+name,verify(p)==folder/name)
            for field,name in (('xyz','system.xyz'),('parameters','amoebabio09.prm')):
                verify(original['original'][field]);verify(original['prepared'][field])
                check(metal+' unchanged '+field,original['original'][field]['sha256']==original['prepared'][field]['sha256']==cell['inputs'][name]['sha256'])
            expected_key,derivation=derive.derive(verify(original['prepared']['key']).read_text())
            check(metal+' exact key derivation',(folder/'system.key').read_text()==expected_key and cell['derivation']==derivation)
            receipt=verify(cell['initialization_receipt']).read_text().splitlines()
            check(metal+' exact raw receipt',[x.split() for x in receipt]==[['PME_OPERATOR_ONLY_INITIALIZATION_V1'],['METAL',metal],['SITES','11958','11958'],['COUNTS','reciprocal_initializations=1','induced_solves=0','energies=0','gradients=0'],['ENERGY_MODEL_ADMISSION','FALSE'],['FULL_STABILITY_ADMISSION','FALSE'],['END']])
            raw={};last=None
            with verify(cell['export']).open() as f:
                for line in f:
                    fields=line.split();last=fields
                    if fields and fields[0] in ('GRID','FFT','CHUNKS','Q_STATS','NORMS'):
                        check(metal+' unique '+fields[0],fields[0] not in raw);raw[fields[0]]=fields[1:]
            check(metal+' raw snapshot denominator',list(map(int,raw['GRID']))==[1,64,64,64,5,11958,11958,262144] and last==['END','11958','179370','262144','1494750'])
            check(metal+' raw native summary',list(map(float,raw['NORMS']))==cell['native_norms'] and list(map(float,raw['Q_STATS']))==cell['native_q_stats'])
            sites=cell['site_checks'];check(metal+' persisted full site/spline checks',sites['passed'] is True and sites['sites']==11958 and sites['spline_scalars']==358740 and sites['coordinate_alpha_comparison']=='exact parsed doubles')
            norm=cell['operator_norms'];comparison=cell['summary_comparison'];ops=3*11958*5**3+128;eps=sys.float_info.epsilon;allow=max(1e-12,32*(ops*eps)/(1-ops*eps))
            check(metal+' nominal-only scope',norm['arithmetic']=='ordinary_double_not_outward' and norm['PSD_admitted'] is False and norm['full_operator_stability_admitted'] is False and norm['q_negative_count']==0 and norm['q_zero']==0 and norm['M']==262144)
            check(metal+' comparison metadata',comparison['passed'] is True and comparison['accumulation_operations']==ops and comparison['norm_relative_allowance']==allow and comparison['q_stat_relative_allowance']==1e-12 and len(comparison['checks'])==10)
            native=cell['native_norms']+cell['native_q_stats']
            for i,key in enumerate(keys):
                row=comparison['checks'][i];error=abs(native[i]-norm[key]);bound=(allow if i<5 else 1e-12)*max(1,abs(native[i]),abs(norm[key]))
                check(metal+' summary arithmetic '+key,row==dict(quantity=key,native=native[i],replay=norm[key],error=error,bound=bound,passed=error<=bound) and row['passed'])
            scale=262144*min(norm['c1']*norm['cinf'],norm['cF'])
            check(metal+' scalar norm propagation',norm['U_R_nominal']==scale*norm['q_abs_max'] and norm['evenness_perturbation_nominal']==scale*norm['q_conjugate_difference_max']/2)
            cells[metal]=dict(operator_norms=norm,summary_comparison=comparison,site_checks=sites,export=cell['export'],initialization_receipt=cell['initialization_receipt'])
        report.update(status='passed_full_box_operator_initialization_actual_review',implementation=result['implementation'],configuration=result['configuration'],source_review=cfg['root_source_review'],root_release=cfg['root_release'],cells=cells,actual_cost={k:result[k] for k in ('wall_seconds','controller_CPU_seconds','child_CPU_seconds','peak_child_RSS_KiB','resource_limits','scheduler_allocated_CPU_seconds')},operation_cost=operation_cost,limitations=['Saved raw identities, native summary tokens and persisted independent replay arithmetic checked; full stencil summation not repeated.','Ordinary-double nominal operator evidence only; no application-error certificate, full stability, response solve, energy, gradient or physical accuracy admission.'])
    except Exception as e:report.update(status='failed_full_box_operator_initialization_actual_review',error=repr(e))
    report.update(checks_passed=len(checks),review_cost=dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,numerical_threads=1,native_calls=0,full_snapshot_replays=0))
    (out/'REVIEW.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(review=pin(out/'REVIEW.json'),status=report['status'],checks=len(checks),cost=report['review_cost'])))
    return 0 if report['status'].startswith('passed_') else 1


if __name__=='__main__':raise SystemExit(main(*sys.argv[1:]))
