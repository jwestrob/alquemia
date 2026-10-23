#!/usr/bin/env python3
"""Immutable prepared-only BIO overlays; no native execution interface."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shlex
import shutil
import time

HERE=Path(__file__).resolve().parent
ROLES={402:73,403:74,405:76,414:85}
METALS={'La':(500,900),'Dy':(501,901)}
CAPABILITIES=('BIO_water_402_403','metal_900_901','response_900_901_402','periodic_PME','dual_direct_fields','conservative_Cartesian_and_frame_forces','unchanged_BIO_nonmetal_terms','Na405_Cl414_self_consistent_polarization','native_metal_salt_mixing')

def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(isinstance(p,dict) and set(p)=={'path','sha256'},'Explicit file pin required');require(pin(p['path'])==p,'Changed pin '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def campaign_fingerprint(c):return digest({k:v for k,v in c.items() if k not in ('root_release','source_review','preparation')})

def physical(values,metal,domain):
    require(len(values)==len(domain['parameters'])==12,'Exact12 parameter domain')
    require([p['name'] for p in domain['parameters']]==['La_a','La_b','La_u','La_v','Dy_a','Dy_b','Dy_u','Dy_v','La_R','La_epsilon','Dy_R','Dy_epsilon'],'Parameter roles/order')
    require(all(type(v) in (int,float) and math.isfinite(v) and p['lower']<=v<=p['upper'] for v,p in zip(values,domain['parameters'])),'Finite in-box vector')
    require(domain['S']==math.log(4.),'Frozen stiffness scale')
    i=0 if metal=='La' else 4;j=8 if metal=='La' else 10;a,b,u,v=values[i:i+4]
    return [a,b,domain['S']*(u+v)/4,domain['S']*(u-v)/4,values[j],values[j+1],domain['fixed_native_POLPAIR'][metal]]

def namespace(prm):
    atoms={};classes=set()
    for line in prm.splitlines():
        if not line.split() or line.split()[0].lower()!='atom':continue
        t=shlex.split(line,comments=True);typ,cls=int(t[1]),int(t[2]);require(typ not in atoms,'Duplicate BIO atom type');atoms[typ]=cls;classes.add(cls)
    require(not ({900,901}&(set(atoms)|classes)),'Mapped metal type/class collision')
    require(all(atoms.get(t)==c for t,c in ROLES.items()),'Ambiguous BIO water/salt type/class roles')
    return {'type_class_roles':ROLES,'reserved_metal_namespace_clear':True}

def render(base,values,domain):
    expected={('vdwpr',73,900),('vdwpr',73,901),('polpair',402,900),('polpair',402,901)};seen=set();out=[];overlay=[]
    for line in base.splitlines(keepends=True):
        t=line.split('#',1)[0].split()
        if not t:out.append(line);continue
        kind=t[0].lower();require(kind not in ('polar-local-response','polar-donor-response','polar-pair-response','external-field','exfld-freq','response-matrix-export'),'Unexpected pre-existing response/field')
        if kind in ('vdwpr','polpair'):
            pair=tuple(sorted(map(int,t[1:3])));key=(kind,*pair)
            require(key in expected and key not in seen,'Ambiguous/unapproved or duplicate pair override');seen.add(key)
            metal='La' if pair[1]==900 else 'Dy';v=physical(values,metal,domain)
            require(len(t)==(5 if kind=='vdwpr' else 4),'Malformed pair row')
            new=f'vdwpr 73 {pair[1]} {v[4]:.17g} {v[5]:.17g}\n' if kind=='vdwpr' else f'polpair 402 {pair[1]} {v[6]:.17g}\n'
            out.append(new);overlay.append(new)
        else:out.append(line)
    require(seen==expected,'Exactly four declared baseline pair rows')
    for metal,(_,mt) in METALS.items():
        v=physical(values,metal,domain);line=f'POLAR-LOCAL-RESPONSE {mt} 402 '+' '.join(f'{x:.17g}' for x in v[:4])+'\n';out.append(line);overlay.append(line)
    return ''.join(out),''.join(overlay)

def candidate_contract(result,config,review,vector,binding,sources):
    require(result.get('status')=='bounded_component_fit_complete' and result.get('mode')=='run' and result.get('error') is None and result.get('candidate_frozen') is True,'Only completed frozen campaign')
    require(result['implementation']==sources['campaign_driver'] and result['configuration']==binding['configuration'],'Exact driver/config binding')
    require(config['implementation']==sources['campaign_driver'] and config['adapter']==sources['campaign_adapter'] and config['core']==sources['campaign_core'],'Exact campaign sources')
    require(config['parent_evaluator']==sources['vacuum_evaluator'] and config['parent_configuration']==sources['parent_configuration'],'Exact scientific parent')
    require(result['configuration_fingerprint']==campaign_fingerprint(config),'Campaign configuration hash')
    require(review.get('status')=='passed_component_fit_actual_review' and review.get('selected_native_admission_passed') is True,'Passed actual terminal candidate review')
    require(review['result']==binding['result'] and review['configuration']==binding['configuration'] and review['selected']==binding['selected'] and review['implementation']==sources['campaign_driver'],'Actual review exact joins')
    key=digest(vector['parameters']);require(result['selected_key']==review['selected_key']==key and result['selected']==binding['selected'],'Selected parameter/artifact hash')
    require(result['vectors'][key]['status']=='complete' and result['vectors'][key]['artifact']==binding['selected'] and result['vectors'][key]['parameters']==vector['parameters'],'Final selected index join')
    require(vector['status']=='complete' and len(vector['cells'])==46 and all(c['status']=='complete' for c in vector['cells'].values()),'Complete selected46, no partial vector')
    return list(vector['parameters'])

def validate_engine_binding(binding,composition):
    """Prospective separate gate only; never launches or changes prepared status."""
    require(isinstance(binding,dict),'No supported engine bound')
    for k in ('source_manifest','executable','review','result'):verify(binding[k])
    review=read(binding['review']);result=read(binding['result'])
    require(review.get('status')=='passed_periodic_BIO_response_static_composition_review','No admitted periodic BIO response capability')
    require(review['source_manifest']==binding['source_manifest'] and review['executable']==binding['executable'] and review['result']==binding['result'],'Exact engine review joins')
    require(review['composition']==composition and result['composition']==composition and result.get('numerical_checks_passed') is True,'Exact prepared-composition numerical admission')
    require(all(review.get('capabilities',{}).get(k) is True for k in CAPABILITIES),'Unsupported salt/custom-response engine composition')
    require(isinstance(review.get('periodic_settings'),dict) and review['periodic_settings'] and review.get('image_convention') and review.get('stability_convention'),'Explicit periodic/image/stability definitions')
    return {'supported_bound_composition':True,'native_execution_performed':False,'review':binding['review']}

def prepare(config_path,out):
    tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(config_path.read_text());sources=read(cfg['inputs'])
    require(cfg['implementation']==pin(__file__) and cfg['inputs']==pin(HERE/'INPUTS.json'),'Adapter/static input pins')
    require(cfg['output_status']=='prepared_only_not_executable' and cfg['engine_binding'] is None,'Input preparation cannot bind or launch an engine')
    audit=read(cfg['source_review']);require(audit['status']=='passed_periodic_model_composition_source_review' and audit['implementation']==cfg['implementation'] and audit['inputs']==cfg['inputs'],'Exact adapter source review')
    for p in sources.values():verify(p)
    binding=cfg['candidate'];require(isinstance(binding,dict) and set(binding)=={'result','configuration','selected','actual_review'},'Explicit future candidate pins required')
    result=read(binding['result']);require(result.get('status')=='bounded_component_fit_complete' and result.get('error') is None,'Reject partial/failed fit before selected artifact read')
    config=read(binding['configuration']);review=read(binding['actual_review']);require(review.get('status')=='passed_component_fit_actual_review' and review['result']==binding['result'],'Terminal review before selected artifact read')
    vector=read(binding['selected']);values=candidate_contract(result,config,review,vector,binding,sources)
    parent=read(sources['parent_configuration']);require(parent['implementation']==sources['vacuum_evaluator'] and parent['domain']==sources['domain'],'Immutable parent/domain join');domain=read(sources['domain'])
    effective={m:physical(values,m,domain) for m in METALS};mapping=read(sources['BIO_mapping']);handoff=read(sources['fixture_inventory']);require(handoff['status']=='periodic_fixture_metadata_handoff_only' and mapping['common_metal_mass_Da']==138.91,'Frozen fixture/mass convention')
    require(len(mapping['aqueous_endpoints'])==len(handoff['fixtures'])==2,'Exact two boxes');prepared=[]
    for row,fixture in zip(mapping['aqueous_endpoints'],handoff['fixtures']):
        require(row['id']==fixture['id'] and row['xyz']==fixture['geometry'] and row['key']==fixture['key'] and row['parameters']==fixture['parameters'],'Exact fixture inputs')
        prm=verify(row['parameters']).read_text();base=verify(row['key']).read_text();verify(row['xyz']);namespace(prm)
        require(row['parameters']==mapping['BIO_parameters_preserved_byte_for_byte'] or row['parameters']['sha256']==mapping['BIO_parameters_preserved_byte_for_byte']['sha256'],'BIO parameter bytes')
        require('POLAR-LOCAL-RESPONSE' not in base and 'EWALD' not in base.upper(),'Baseline mapping key only')
        key,overlay=render(base,values,domain);prepared.append((row,key,overlay))
    require(not out.exists(),'Immutable new output only');out.mkdir(parents=True)
    files=[]
    for row,key,overlay in prepared:
        folder=out/row['id'];folder.mkdir();shutil.copyfile(verify(row['xyz']),folder/'system.xyz');shutil.copyfile(verify(row['parameters']),folder/'amoebabio09.prm')
        (folder/'system.prepared.key').write_text(key);(folder/'metal_water.overlay.key').write_text(overlay)
        files.append({'id':row['id'],'metal':row['metal'],'original':row['files'],'prepared':{k:pin(folder/n) for k,n in [('xyz','system.xyz'),('parameters','amoebabio09.prm'),('key','system.prepared.key'),('overlay','metal_water.overlay.key')]}})
    manifest={'status':'prepared_only_not_executable','implementation':pin(__file__),'configuration':pin(config_path),'inputs':cfg['inputs'],'source_review':cfg['source_review'],'candidate':binding,'selected_parameter_key':digest(values),'selected_parameters':values,'effective_a_b_c_d_R_epsilon_POLPAIR':effective,'files':files,'BIO_intramolecular_and_salt_parameters_unchanged':True,'periodic_settings_bound':False,'engine_binding':None,'native_capabilities_required':list(CAPABILITIES),'common_metal_mass_Da':138.91,'source_Dy_mass_Da':162.5,'mass_convention_requires_dynamics_decision':True,'native_calls':0,'physical_accuracy_qualified':False,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n');return manifest

def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args();r=prepare(a.config,a.out);print(json.dumps({'status':r['status'],'manifest':pin(a.out/'MANIFEST.json'),'native_calls':0}))

if __name__=='__main__':main()
