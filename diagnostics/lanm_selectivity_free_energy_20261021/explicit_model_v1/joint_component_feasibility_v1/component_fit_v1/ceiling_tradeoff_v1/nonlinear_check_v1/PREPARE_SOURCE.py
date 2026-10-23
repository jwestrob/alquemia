#!/usr/bin/env python3
"""Version the admitted finite pilot; scientific execution methods stay identical."""
import ast
import difflib
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
OLD=HERE.parents[2]/'nonlinear_check_v1/RUN.py'
EXPECTED='d9d058dad0bb695c4a32848b61e2106a55604a8d14465ebff65875242714eccd'

def main():
    before=OLD.read_text();assert hashlib.sha256(OLD.read_bytes()).hexdigest()==EXPECTED
    text=before
    a=text.index('def context(cfg):');b=text.index('\n\nclass NotAttempted',a)
    text=text[:a]+(HERE/'CONTEXT_SOURCE.txt').read_text().rstrip()+text[b:]
    a=text.index('def component_report(');b=text.index('\n\ndef main(',a)
    text=text[:a]+(HERE/'COMPONENT_SOURCE.txt').read_text().rstrip()+text[b:]
    text=text.replace("('root_release','source_review','preparation')","('root_release','source_review','preparation','tradeoff_actual_review')")
    text=text.replace('fixed_half','ceiling_tradeoff')
    text=text.replace("('parent_configuration','half_review','origin','declaration','source_review','preparation','candidate','limits','executor')",
                      "('parent_configuration','tradeoff_result','tradeoff_actual_review','followup_policy','origin','declaration','source_review','preparation','candidate','limits','executor')")
    old="            review=read(self.batch['source_review']);prep=read(self.batch['preparation']);release=read(self.root_release)"
    new=old+'''\n            actual=read(self.batch['tradeoff_actual_review']);verify(actual['implementation'])
            expected={'fraction':self.batch['candidate']['fraction'],'physical_key':self.batch['candidate']['key'],'parameters':self.batch['candidate']['parameters']}
            source=read(self.batch['tradeoff_result'])
            require(actual['status']=='passed_ceiling_tradeoff_actual_arithmetic_review' and actual['certificate_checks_passed'] is True
                    and actual['result']==self.batch['tradeoff_result'] and actual['prospective_policy']==self.batch['followup_policy']
                    and actual['prospective_trigger_fraction']==1. and actual['prospective_trigger_witness']==expected
                    and actual['solver_status']==source['solver_status'] and actual['solver_success']==source['solver_success'],
                    'Exact actual arithmetic/policy admission, independent of solver success')'''
    assert text.count(old)==1;text=text.replace(old,new)
    text=text.replace("origin=cfg['origin'],half_review=cfg['half_review'],candidate=cfg['candidate']",
                      "origin=cfg['origin'],tradeoff_actual_review=cfg.get('tradeoff_actual_review'),tradeoff_result=cfg['tradeoff_result'],candidate=cfg['candidate']")
    a=ast.parse(before);b=ast.parse(text)
    def functions(tree):return {n.name:ast.dump(n,include_attributes=False) for n in ast.walk(tree) if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    aa,bb=functions(a),functions(b)
    unchanged=('allowed','_stop','_invoke','_cell','evaluate','snapshot')
    assert all(aa[n]==bb[n] for n in unchanged)
    # Literal paired energy cycles stay identical; prediction/reporting alone changes.
    assert aa['one']==bb['one']
    out=HERE/'RUN.py';assert not out.exists();out.write_text(text)
    (HERE/'RUN.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),text.splitlines(True),fromfile=str(OLD),tofile=str(out))))
    report={'status':'prepared_versioned_pilot_source','old_source':{'path':str(OLD.resolve()),'sha256':EXPECTED},
            'source':{'path':str(out.resolve()),'sha256':hashlib.sha256(out.read_bytes()).hexdigest()},
            'identical_AST_methods':list(unchanged)+['one'],'native_calls':0}
    (HERE/'SOURCE_DELTA.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
