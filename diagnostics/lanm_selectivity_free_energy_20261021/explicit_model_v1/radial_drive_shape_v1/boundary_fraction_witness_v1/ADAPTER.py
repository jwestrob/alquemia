"""One exact vector, no imported cells; qualified radial/native/journal paths."""
import ast
import copy
import os
from pathlib import Path
import types


def require(ok, message):
    if not ok: raise ValueError(message)


def bounded_parent(old):
    tree = ast.parse(Path(old.__file__).read_text())
    fn = copy.deepcopy(next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'evaluator_class'))
    hits = {200: 0, 6496: 0}
    for node in ast.walk(fn):
        if isinstance(node, ast.Constant) and type(node.value) is int and node.value in hits:
            before = node.value; hits[before] += 1; node.value = {200: 1, 6496: 56}[before]
    require(hits == {200: 2, 6496: 3}, 'Exact parent budget-only substitutions')
    cell = next(n for n in ast.walk(fn) if isinstance(n, ast.FunctionDef) and n.name == '_cell'); changed = 0
    for node in ast.walk(cell):
        if isinstance(node, ast.Dict):
            for i, key in enumerate(node.keys):
                if isinstance(key, ast.Constant) and key.value == 'job_id':
                    require(isinstance(node.values[i], ast.Constant) and node.values[i].value is None, 'Original scheduler metadata')
                    node.values[i] = ast.parse("self.executor_identity['job_id']", mode='eval').body; changed += 1
    require(changed == 1, 'One truthful cell scheduler metadata substitution')
    ns = dict(vars(old)); exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])), __file__, 'exec'), ns)
    return types.SimpleNamespace(**ns)


def evaluator_class(old, radial, base, checks, release_gate):
    bounded = bounded_parent(old)
    facade = types.SimpleNamespace(**vars(base))
    facade.effective = lambda v, m, d: base.effective(v[:12], m, d)+[v[12 if m == 'La' else 13]]
    Parent = bounded.evaluator_class(facade)
    Radial = radial.evaluator_class(bounded, base, checks, release_gate)

    class Evaluator(Radial):
        def __init__(self, cfg, path, point, parent_config, allocation):
            require(not Path(path).exists(), 'Fresh one-point cache')
            self.batch = copy.deepcopy(cfg); self.mode = 'run'; self.batch_fingerprint = cfg['fingerprint']
            self.allowed = {point['key']:dict(parameters=point['parameters'], role='fixed_fraction_witness', fraction=.0625)}
            self.complete_vectors = {}; self.imported_vector_keys = set(); self.registry = {}; self.loaded_saved_vectors = {}
            self.saved_metal_index = {}; self.old_config = copy.deepcopy(parent_config)
            self.executor_identity = dict(allocation, controller_pid=os.getpid(), run_id=cfg['run_id'], implementation=cfg['implementation'], adapter=cfg['adapter'], parent_implementation=cfg['parent_evaluator'])
            runtime = copy.deepcopy(parent_config)
            runtime.update(source_manifest=cfg['native_source_manifest'], executables=cfg['executables'], build_receipt=cfg['build'], qualification_review=cfg['qualification_review'])
            base.Evaluator.__init__(self, runtime, path, root_release=cfg.get('root_release'))
            self.initial_vector += [0., 0.]
            self.family['radial_extension'] = dict(append_order=['La_radial_drive','Dy_radial_drive'], bounds=[[-1.,1.],[-1.,1.]], original12_family_preserved=True)
            self.program_attempts = []; self.cell_attempts = {}; self.started_cells = set(); self.imported_cell_keys = set()
            require(len(self.m['rows']) == 46 and sum(1+int(r['needs_gradient']) for r in self.m['rows'].values()) == 56, 'Exact46 cells/56 programs')
            require(base.digest(point['parameters']) == point['key'] and self.checked_values(point['parameters']) == point['parameters'], 'One exact literal point')
            for metal in ('La','Dy'):
                require(facade.effective(point['parameters'],metal,self.domain) != facade.effective(self.initial_vector,metal,self.domain), 'No implicit baseline fallback')
            self.ledger.update(limits=copy.deepcopy(cfg['limits']), new_policy=cfg['adapter'], configuration_fingerprint=self.batch_fingerprint,
                               native_returncodes_observed=0, initial_imported_unique_cells=0, available_saved_vectors=0)
            self.persistence = old.J.Journal(self.cache_dir/'EVENTS.jsonl', self.ledger, self.lock,
                                            dict(implementation=cfg['adapter'], configuration_fingerprint=self.batch_fingerprint, root_release=cfg.get('root_release')))
            self.persistence.materialize(self.ledger_path, expected=self.ledger)
            require(self.required_programs(point['parameters']) == 56, 'Full one-vector native denominator')

        def _gate(self):
            require(self.mode == 'run' and not self.stop_event.is_set(), 'No native after stop')
            require(self.ledger['native_program_invocations'] < 56, 'Hard56 invocation boundary')
            require(all(os.environ.get(k) == '1' for k in self.batch['threads']) and os.environ.get('SLURM_JOB_ID') == self.executor_identity['job_id']
                    and os.environ.get('SLURM_CPUS_PER_TASK') == '13' and os.environ.get('SLURM_NTASKS') == '1', 'Actual13CPU Slurm and one-thread native')
            if not self._gated: release_gate(self.batch); self._gated = True

        def evaluate(self, values):
            values = self.checked_values(values); key = base.digest(values)
            require(key in self.allowed and self.allowed[key]['parameters'] == values, 'Only exact prepared1/16 point')
            require(not self.ledger['vector_keys'] and self.ledger['native_program_invocations'] == 0 and self.required_programs(values) == 56, 'One request/no retry/full56 budget')
            answer = Parent.evaluate(self, values)
            require(self.ledger['native_program_invocations'] == 56 and len(answer['cells']) == 46, 'Complete46/56 witness')
            return answer

        def snapshot(self):
            with self.lock: self.persistence.materialize(self.ledger_path, expected=self.ledger)
            s = base.Evaluator.snapshot(self)
            s.update(max_parameter_vectors=1, max_native_programs=56, active_limits=self.batch['limits'], new_parameter_vectors=len(self.ledger['vector_keys']),
                     native_returncodes_observed=self.ledger['native_returncodes_observed'], imported_unique_cells=0, available_saved_vectors=0,
                     old_mutable_ledger_read=False, executor=copy.deepcopy(self.executor_identity))
            return s
    return Evaluator
