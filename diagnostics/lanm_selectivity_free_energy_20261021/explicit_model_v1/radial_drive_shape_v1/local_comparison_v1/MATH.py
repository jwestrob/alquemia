"""Dimension-only views of admitted affine/physical arithmetic; no solver here."""
import ast,copy,math,types
from pathlib import Path
import numpy as np

def require(ok,msg):
    if not ok:raise ValueError(msg)

def clone(module,ns,name,replacements):
    fn=copy.deepcopy(next(n for n in ast.parse(Path(module.__file__).read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name));hits={k:0 for k in replacements}
    for node in ast.walk(fn):
        if isinstance(node,ast.Constant) and type(node.value) is int and node.value in replacements:
            old=node.value;node.value=replacements[old];hits[old]+=1
    require(all(hits.values()),'Expected dimension constants in '+name)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),str(Path(__file__)),'exec'),ns)
    return hits

def core14(original):
    ns=dict(vars(original))
    # Rebind unchanged pure helpers so checked()/measure() use the14-coordinate view.
    for name,fn in vars(original).items():
        if isinstance(fn,types.FunctionType) and fn.__globals__ is vars(original):ns[name]=types.FunctionType(fn.__code__,ns,name,fn.__defaults__,fn.__closure__)
    edits={name:clone(original,ns,name,repl) for name,repl in [('checked',{12:14}),('linearize',{24:28,12:14}),('certify_fraction',{13:15,12:14})]}
    def points(values,domain):
        require(len(values)==14 and values[12:]==[0.,0.],'Exact zero radial center only')
        old_domain=dict(domain,parameters=domain['parameters'][:12])
        result=[p+[0.,0.] for p in original.probe_points(values[:12],old_domain)]
        for j in (12,13):
            for sign in (-1,1):
                p=list(values);p[j]=sign*.01;result.append(p)
        return result
    ns['probe_points']=points
    conic=original.conic;cns=dict(vars(conic));edits['conic_assemble']=clone(conic,cns,'assemble',{12:14,13:15,856:860,1817:1821});ns['conic']=types.SimpleNamespace(**cns)
    result=types.SimpleNamespace(**ns);result.dimension_edits=edits;return result

def attribution14(original):
    ns=dict(vars(original));edits={}
    for name in ('affine','certificate','assemble'):
        edits[name]=clone(original,ns,name,{12:14,13:15} if name=='assemble' else {12:14})
    # certificate calls affine through these shared globals.
    result=types.SimpleNamespace(**ns);result.dimension_edits=edits;return result

def extend_domain(domain):
    result=copy.deepcopy(domain)
    result['parameters']+= [{'name':m+'_radial_drive','lower':-1.,'upper':1.,'initial':0.} for m in ('La','Dy')]
    return result

def view(vector,pin,bridge):
    require(len(vector['parameters'])==12 and vector['status']=='complete','Admitted original complete vector')
    result=dict(vector,parameters=vector['parameters']+[0.,0.])
    result['radial_zero_view']={'original_artifact':pin,'original_parameters':vector['parameters'],'bridge':bridge,'raw_cells_unchanged':True,'new_native_calls':0}
    return result

def problem(core,attr,center,probes,fixed,targets):
    p=attr.normalized_problem(core,center,probes,fixed);p['electrical_targets']=targets
    c=p['electrical_unweighted']['C_dipoles_Debye'];dy=c['names'].index('aquo9_central_Dy');la=c['names'].index('aquo9_central_La')
    reference=targets['C_Debye']['aquo9_central_Dy'][2]-targets['C_Debye']['aquo9_central_La'][2]
    r=math.fsum([float(c['r'][dy,2]),-float(c['r'][la,2]),-reference]);J=c['J'][dy,2]-c['J'][la,2]
    p['signed_W9_z']={'r_Debye':r,'J_Debye':J,'absolute_ceiling_Debye':abs(r),'reference_Dy_minus_La_z_Debye':reference}
    return p
