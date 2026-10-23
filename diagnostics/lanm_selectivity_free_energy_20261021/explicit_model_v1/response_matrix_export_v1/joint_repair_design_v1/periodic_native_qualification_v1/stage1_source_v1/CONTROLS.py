"""Bounded source and manufactured scalar/sparse checks. Never evaluate a box."""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import tempfile
import time
import numpy as np
from scipy.sparse import csc_array,csr_array
import DERIVE_SETUP
import PREPARE
import REPLAY


def run(out):
    tick,cpu=time.monotonic(),time.process_time();here=Path(__file__).parent;checks=[]
    def check(label,condition):
        checks.append(dict(label=label,passed=bool(condition)))
        if not condition:raise ValueError(label)
    def rejected(label,call):
        try:call()
        except ValueError:check(label,True)
        else:check(label,False)
    for p in here.glob('*.py'):
        ast.parse(p.read_text());check('Python syntax '+p.name,True)
    provenance=json.loads((here/'SETUP_PROVENANCE.json').read_text())
    original=PREPARE.verify(provenance['parent_induce']).read_text()
    old,new=DERIVE_SETUP.derive(original)
    check('Exact derived setup',new==(here/'PARITY_SETUP.f').read_text())
    check('Exact saved setup authority',old==(here/'SETUP_AUTHORITY.f').read_text())
    inverse=new.replace('parityresponsesetup (handled)','periodicresponsesetup (handled)',1)
    inverse=inverse.replace('c     Standalone explicit-operator initialization only.\n'
                           'c     Normal setup and both fatal guards remain unchanged.\n',
                           'c     This source has no admission of a real periodic execution.\n'
                           '      call periodicresponsebackendguard\n')
    check('Only declared setup name/end delta',inverse==old)
    for guard in ('periodicresponsebackendguard','periodicresponsestability'):
        body=original.split('      subroutine '+guard+'\n',1)[1].split('      end\n',1)[0]
        check('Unchanged normal fatal guard '+guard,'      call fatal\n' in body)
    harness=(here/'PARITY.f90').read_text()
    calls=re.findall(r'^\s*call\s+(\w+)',harness,re.I|re.M)
    check('No implicit solver/energy/gradient call',not set(x.lower() for x in calls)&
          {'induce','induce0a','gradient','energy','analysis','epolar','epolar1','epolar3','localresponsesetup'})
    check('One declared direct initialization call site',calls.count('dfield0c')==1)
    check('One full action call site',calls.count('ufield0')==1)
    check('One separate reciprocal call site',calls.count('reciprocal_apply')==1)
    check('Exact finite mode/probe loops','do mode=1,3' in harness and 'do col=1,10' in harness)
    check('Native H computed before actions',harness.index('call parityresponsesetup')<harness.index('call ufield0'))
    check('Expected state before actions',harness.index('Native H differs')<harness.index('call ufield0'))
    check('No copied static H assignment',not re.search(r'lrbeta\s*\([^\n]*\)\s*=',harness,re.I))
    check('No free source line truncation',all(len(s)<=132 for s in harness.splitlines()))
    check('Fixed setup executable line width',all(len(s)<=72 for s in new.splitlines() if s[:1].lower()!='c'))
    # Explicit independent sparse action on a manufactured 2x2 symmetric matrix.
    lower=csc_array(np.array([[2.,0.],[.25,3.]]));v=np.array([1+2j,-3+.5j])
    check('Sparse lower complex action',np.array_equal(REPLAY.real_action(lower,v),np.array([[2.,.25],[.25,3.]])@v))
    # Non-even q must retain imaginary coupling, rather than using q-even or real-only FFT.
    B=csr_array(np.array([[1.,0.],[0.,1.],[.25,.5],[.75,-.25]]))
    s={'shape':(4,1,1),'qfac':np.array([0.,1.,2.,4.]).reshape(4,1,1)}
    z=np.array([1+2j,-.5+1j]);actual=REPLAY.complex_reciprocal(B,s,z)
    F=np.exp(-2j*np.pi*np.outer(np.arange(4),np.arange(4))/4)
    expected=B.T@(F.conj().T@(s['qfac'].ravel()*(F@(B@z))))
    check('Complex non-even reciprocal arithmetic',np.max(np.abs(actual-expected))<1e-13)
    even=dict(s,qfac=((s['qfac'].ravel()+s['qfac'].ravel()[[0,3,2,1]])/2).reshape(4,1,1))
    check('Odd q materially retained',np.max(np.abs(actual-REPLAY.complex_reciprocal(B,even,z)))>.1)
    trace=[];REPLAY.compare([1.],[1.],'exact',trace)
    check('Exact comparison passes',trace[-1]['passed'])
    rejected('Nonfinite action rejected',lambda:REPLAY.compare([np.nan],[1.],'nan',[]))
    rejected('Wrong shape rejected',lambda:REPLAY.compare([1,2],[1],'shape',[]))
    rejected('Outside unchanged parity gate rejected',lambda:REPLAY.compare([1.000001],[1.],'wrong',[]))
    sites=[dict(type=900,xyz=[0.,0.,0.],alpha=1.6),dict(type=402,xyz=[2.5,0.,0.],alpha=.8)]
    damping=np.array([1.6,.8])**(1/6);gamma=np.array([0.,.17])
    h,b=REPLAY.compact_response(sites,[0.,0.,0.,0.],gamma,damping)
    check('Zero h/H manufactured bypass',np.array_equal(h,np.zeros((2,3))) and np.array_equal(b,np.zeros(2)))
    h,b=REPLAY.compact_response(sites,[.1,0.,0.,0.],gamma,damping)
    check('Drive-only native-unit support',h[1,0]>0 and np.array_equal(h[0],np.zeros(3)) and np.array_equal(b,np.zeros(2)))
    h,b=REPLAY.compact_response(sites,[0.,0.,.1,0.],gamma,damping)
    check('Stiffness-only primary block',np.array_equal(h,np.zeros((2,3))) and b[0]>0 and b[1]==0)
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)/'raw';p.write_text('ROW 1\nEXTRA\n');rows=REPLAY.Rows(p)
        check('Strict row accepts exact schema',rows.row('ROW',1)==['1'])
        rejected('Trailing raw record rejected',rows.end)
        rows.stream.close()
        p.write_text('BROKEN\n');rows=REPLAY.Rows(p)
        rejected('Malformed raw row rejected',lambda:rows.row('ROW',1));rows.stream.close()
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    result=dict(status='passed_stage1_source_and_manufactured_controls',implementation=PREPARE.pin(__file__),
                source={p.name:PREPARE.pin(p) for p in here.iterdir() if p.is_file()},checks=checks,count=len(checks),
                actual_box_evaluations=0,native_calls=0,builds=0,factorizations=0,actual_operator_applications=0,
                cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(result['status'],len(checks))


if __name__=='__main__':
    import sys
    run(sys.argv[1])
