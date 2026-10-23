"""Finite saved Stage2 collection. No native, FFT, factor or solver call."""
import hashlib
import json
import math
import re
import sys
import time
from decimal import Decimal, localcontext
from pathlib import Path

N=11958
DIM=3*N
STEPS=(Decimal('0.0001'),Decimal('0.00005'))
GATE=Decimal('0.001')
EPS=Decimal.from_float(sys.float_info.epsilon)
BYPASSES=('BYPASS_OLD_ABSENT','BYPASS_CLONE_ABSENT','BYPASS_CLONE_ZERO')
DIRECTION_IDS=('metal_0','metal_1','metal_2','water_rotation','whole_box','strain_0','strain_1','strain_2')


def require(ok,message):
    if not ok:raise ValueError(message)


def pin(path):
    path=Path(path).resolve()
    with path.open('rb') as stream:sha=hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(path=str(path),sha256=sha)


def verify(item):
    require(isinstance(item,dict) and set(item)=={'path','sha256'} and pin(item['path'])==item,'Missing/changed pin')
    return Path(item['path'])


def read(item):return json.loads(verify(item).read_text())


def number(token):
    require(re.fullmatch(r'[+-]?\d+\.\d+(?:[EeDd][+-]?\d+)?',token) is not None,'Malformed native decimal')
    value=Decimal(token.replace('D','E').replace('d','E'))
    require(value.is_finite() and math.isfinite(float(value)),'Nonfinite native value')
    return value


def rounding(token):
    value=number(token)
    return Decimal(5).scaleb(value.as_tuple().exponent-1)


def arithmetic(dim,*values):return Decimal(128)*dim*EPS*max(Decimal(1),*(abs(x) for x in values))


def verdict(label,error,bound,**extra):
    require(error.is_finite() and bound.is_finite() and error>=0 and bound>=0,'Finite numerical check')
    return dict(label=label,error=float(error),bound=float(bound),passed=error<=bound,**extra)


def _parse_v1(path,n=N):
    """Single strict read of full output; no implicit displaced/last-solve choice."""
    with Path(path).open() as stream,localcontext() as ctx:
        ctx.prec=70
        def literal(value):require(stream.readline().strip()==value,'Native line: '+value)
        def row(prefix,width):
            fields=stream.readline().split();head=prefix.split()
            require(fields[:len(head)]==head and len(fields)==len(head)+width,'Native record: '+prefix)
            return fields[len(head):]
        literal('STAGE2_ONE_ANALYTIC_GRADIENT_V1')
        literal('UNITS E=kcal/mol G=positive_dE/dr_kcal/mol/angstrom VIRIAL=kcal/mol')
        energy=row('TOTAL_ENERGY',1)[0];ep=row('POLARIZATION_ENERGY',1)[0]
        electric,dielec=row('CONTRACTION_CONSTANTS electric dielec',2)
        require(number(dielec)>0 and number(electric)>0,'Positive contraction constants')
        number(energy);number(ep);g=[];vir=[]
        for i in range(1,n+1):
            fields=row('TOTAL_GRADIENT',4);require(fields[0]==str(i),'Gradient atom order')
            for token in fields[1:]:number(token)
            g.extend(fields[1:])
        for i in range(1,4):
            fields=row('TOTAL_VIRIAL',4);require(fields[0]==str(i),'Virial row order')
            for token in fields[1:]:number(token)
            vir.extend(fields[1:])
        role=row('ROLE',1)[0];state_id=row('STATE_ID',1)[0];binding=row('INPUT_BINDING',1)[0]
        require(role in ('ACTIVE',)+BYPASSES,'Declared gradient role')
        require(re.fullmatch('[a-f0-9]{64}',binding) is not None,'Input binding digest')
        result=dict(role=role,state_id=state_id,input_binding=binding,energy_token=energy,polarization_token=ep,
                    electric_token=electric,dielec_token=dielec,gradient_tokens=g,virial_tokens=vir,
                    E_G_virial_token_sha256=hashlib.sha256(('\n'.join([energy]+g+vir)+'\n').encode()).hexdigest(),
                    counters=None,PCG_pre_peek=None,shared_h=None,polarization_contraction=None)
        if role=='BYPASS_OLD_ABSENT':
            literal('INTERNAL_COUNTERS unavailable_unchanged_old_archive')
        else:
            counts=row('COUNTS induced_solve direct mutual',3)
            require(all(x.isdecimal() for x in counts),'Nonnegative counters')
            solve,direct,mutual=map(int,counts)
            require(solve==direct==1 and mutual>=1,'One direct/solve; actual mutual applications')
            literal('COUNTS extra_diagnostic_operator_applications=0 internal_FD=0')
            p=row('PCG_PRE_PEEK',5);require(p[0].isdecimal(),'PCG iteration count')
            iterations=int(p[0]);eps,epsd,epsp,old=map(number,p[1:])
            require(3<=iterations<100 and 0<=eps<Decimal('1e-8') and eps<=old and epsd>=0 and epsp>=0,
                    'Original native induction termination')
            result['counters']=dict(induced_solves=solve,direct_initializations=direct,mutual_applications=mutual,
                                    extra_diagnostic_operator_applications=0,internal_FD=0)
            result['PCG_pre_peek']=dict(iterations=iterations,RMS_Debye=float(eps),squared_channel_residuals=[float(epsd),float(epsp)],
                                       previous_RMS_Debye=float(old),tokens=p,post_peek_residual_not_recomputed=True)
            if role=='ACTIVE':
                reported=number(row('SHARED_H_READBACK_MAX_ERROR',1)[0]);require(reported>=0,'Nonnegative shared-h readback')
                hook_error=Decimal(0);scale=Decimal(1);contraction=Decimal(0);dual_max=Decimal(0)
                for i in range(1,n+1):
                    records=[]
                    for prefix,width in (('RAW_DIRECT',7),('CORRECTED_RHS',7),('FINAL_POST_PEEK',8),('CORRECTION',5)):
                        tokens=row(prefix,width);require(tokens[0]==str(i),'Dual/correction site order')
                        records.append(list(map(number,tokens[1:])))
                    raw,rhs,final,correction=records;alpha=final[0];u=final[1:4];v=final[4:];h=correction[1:]
                    require(alpha>=Decimal('1e-8'),'Positive actual intrinsic alpha')
                    scale=max(scale,*(abs(x) for x in raw+rhs+h))
                    hook_error=max(hook_error,*(abs(rhs[j]-raw[j]-h[j%3]) for j in range(6)))
                    dual_max=max(dual_max,*(abs(u[j]-v[j]) for j in range(3)))
                    contraction+=sum((u[j]*rhs[3+j] for j in range(3)),Decimal(0))
                bound=arithmetic(3*n,scale)
                result['shared_h']=verdict('same h in both direct RHS',max(hook_error,reported),bound,
                                            reconstructed_max_error=float(hook_error),reported_max_error=float(reported))
                prediction=-number(electric)/(2*number(dielec))*contraction
                result['polarization_contraction']=verdict('stationary dual polarization energy',abs(number(ep)-prediction),
                    Decimal('1e-8')+arithmetic(3*n,number(ep),prediction),native_kcal_mol=float(number(ep)),
                    contraction_kcal_mol=float(prediction),maximum_u_v_difference=float(dual_max),
                    convention='-electric/(2*dielec) sum_i u_i dot corrected_Ep_i; final post-peek u/v')
        literal('END gradient_calls=1 energy_only_calls=0 internal_FD=0')
        require(not stream.read().strip(),'Trailing native output')
        return result


def _parse_v2(path,n=N):
    """Single strict read of full output; no implicit displaced/last-solve choice."""
    with Path(path).open() as stream,localcontext() as ctx:
        ctx.prec=70
        def literal(value):require(stream.readline().strip()==value,'Native line: '+value)
        def row(prefix,width):
            fields=stream.readline().split();head=prefix.split()
            require(fields[:len(head)]==head and len(fields)==len(head)+width,'Native record: '+prefix)
            return fields[len(head):]
        literal('STAGE2_ONE_ANALYTIC_GRADIENT_V1')
        literal('UNITS E=kcal/mol G=positive_dE/dr_kcal/mol/angstrom VIRIAL=kcal/mol')
        energy=row('TOTAL_ENERGY',1)[0];ep=row('POLARIZATION_ENERGY',1)[0]
        electric,dielec=row('CONTRACTION_CONSTANTS electric dielec',2)
        require(number(dielec)>0 and number(electric)>0,'Positive contraction constants')
        number(energy);number(ep);g=[];vir=[]
        for i in range(1,n+1):
            fields=row('TOTAL_GRADIENT',4);require(fields[0]==str(i),'Gradient atom order')
            for token in fields[1:]:number(token)
            g.extend(fields[1:])
        for i in range(1,4):
            fields=row('TOTAL_VIRIAL',4);require(fields[0]==str(i),'Virial row order')
            for token in fields[1:]:number(token)
            vir.extend(fields[1:])
        literal('STATE_SCHEMA PERIODIC_PARITY_STATE_V2')
        radial={}
        for metal in (900,901):
            fields=row('RADIAL_COEFFICIENT',2)
            require(fields[0]==str(metal),'Both ordered radial slots')
            value=number(fields[1]);require(abs(value)<=1,'Radial domain')
            radial[str(metal)]=float(value)
        role=row('ROLE',1)[0];state_id=row('STATE_ID',1)[0];binding=row('INPUT_BINDING',1)[0]
        require(role in ('ACTIVE',)+BYPASSES,'Declared gradient role')
        require(re.fullmatch('[a-f0-9]{64}',binding) is not None,'Input binding digest')
        result=dict(role=role,state_id=state_id,input_binding=binding,energy_token=energy,polarization_token=ep,
                    electric_token=electric,dielec_token=dielec,gradient_tokens=g,virial_tokens=vir,
                    E_G_virial_token_sha256=hashlib.sha256(('\n'.join([energy]+g+vir)+'\n').encode()).hexdigest(),
                    counters=None,PCG_pre_peek=None,shared_h=None,polarization_contraction=None)
        if role=='BYPASS_OLD_ABSENT':
            literal('INTERNAL_COUNTERS unavailable_unchanged_old_archive')
        else:
            counts=row('COUNTS induced_solve direct mutual',3)
            require(all(x.isdecimal() for x in counts),'Nonnegative counters')
            solve,direct,mutual=map(int,counts)
            require(solve==direct==1 and mutual>=1,'One direct/solve; actual mutual applications')
            literal('COUNTS extra_diagnostic_operator_applications=0 internal_FD=0')
            p=row('PCG_PRE_PEEK',5);require(p[0].isdecimal(),'PCG iteration count')
            iterations=int(p[0]);eps,epsd,epsp,old=map(number,p[1:])
            require(3<=iterations<100 and 0<=eps<Decimal('1e-8') and eps<=old and epsd>=0 and epsp>=0,
                    'Original native induction termination')
            result['counters']=dict(induced_solves=solve,direct_initializations=direct,mutual_applications=mutual,
                                    extra_diagnostic_operator_applications=0,internal_FD=0)
            result['PCG_pre_peek']=dict(iterations=iterations,RMS_Debye=float(eps),squared_channel_residuals=[float(epsd),float(epsp)],
                                       previous_RMS_Debye=float(old),tokens=p,post_peek_residual_not_recomputed=True)
            if role=='ACTIVE':
                reported=number(row('SHARED_H_READBACK_MAX_ERROR',1)[0]);require(reported>=0,'Nonnegative shared-h readback')
                hook_error=Decimal(0);scale=Decimal(1);contraction=Decimal(0);dual_max=Decimal(0)
                for i in range(1,n+1):
                    records=[]
                    for prefix,width in (('RAW_DIRECT',7),('CORRECTED_RHS',7),('FINAL_POST_PEEK',8),('CORRECTION',5)):
                        tokens=row(prefix,width);require(tokens[0]==str(i),'Dual/correction site order')
                        records.append(list(map(number,tokens[1:])))
                    raw,rhs,final,correction=records;alpha=final[0];u=final[1:4];v=final[4:];h=correction[1:]
                    require(alpha>=Decimal('1e-8'),'Positive actual intrinsic alpha')
                    scale=max(scale,*(abs(x) for x in raw+rhs+h))
                    hook_error=max(hook_error,*(abs(rhs[j]-raw[j]-h[j%3]) for j in range(6)))
                    dual_max=max(dual_max,*(abs(u[j]-v[j]) for j in range(3)))
                    contraction+=sum((u[j]*rhs[3+j] for j in range(3)),Decimal(0))
                bound=arithmetic(3*n,scale)
                result['shared_h']=verdict('same h in both direct RHS',max(hook_error,reported),bound,
                                            reconstructed_max_error=float(hook_error),reported_max_error=float(reported))
                prediction=-number(electric)/(2*number(dielec))*contraction
                result['polarization_contraction']=verdict('stationary dual polarization energy',abs(number(ep)-prediction),
                    Decimal('1e-8')+arithmetic(3*n,number(ep),prediction),native_kcal_mol=float(number(ep)),
                    contraction_kcal_mol=float(prediction),maximum_u_v_difference=float(dual_max),
                    convention='-electric/(2*dielec) sum_i u_i dot corrected_Ep_i; final post-peek u/v')
        literal('END gradient_calls=1 energy_only_calls=0 internal_FD=0')
        require(not stream.read().strip(),'Trailing native output')
        result['radial']=radial
        result['state_schema']='PERIODIC_PARITY_STATE_V2'
        return result


def parse_gradient(path,role,manifest,n=N):
    # Role selects grammar prospectively. No auto-detection or fallback.
    require(role in ('ACTIVE',)+BYPASSES,'Bound gradient role')
    parsed=_parse_v1(path,n) if role=='BYPASS_OLD_ABSENT' else _parse_v2(path,n)
    require(parsed['role']==role,'Exact role readback')
    if role=='ACTIVE':require(parsed['radial']==manifest['model']['radial'],'Both actual radial coefficients')
    elif role!='BYPASS_OLD_ABSENT':require(parsed['radial']=={'900':0.,'901':0.},'Bypass radial zeros')
    return parsed
