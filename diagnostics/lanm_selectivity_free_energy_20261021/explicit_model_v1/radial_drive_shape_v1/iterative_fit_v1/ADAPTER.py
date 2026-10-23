"""Finite radial14 campaign policy; inherited native/parser/journal implementations."""
import ast
import copy
import hashlib
import os
from pathlib import Path
import socket
import types


def require(ok, message):
    if not ok: raise ValueError(message)


def bounded_parent(old):
    tree = ast.parse(Path(old.__file__).read_text())
    fn = copy.deepcopy(next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'evaluator_class'))
    counts = {200: 0, 6496: 0}
    for n in ast.walk(fn):
        if isinstance(n, ast.Constant) and type(n.value) is int and n.value in counts:
            before = n.value; counts[before] += 1; n.value = {200: 228, 6496: 7280}[before]
    require(counts == {200: 2, 6496: 3}, 'Exact parent budget-only substitutions')
    cell=next(n for n in ast.walk(fn) if isinstance(n,ast.FunctionDef) and n.name=='_cell'); changed=0
    for node in ast.walk(cell):
        if isinstance(node,ast.Dict):
            for i,key in enumerate(node.keys):
                if isinstance(key,ast.Constant) and key.value=='job_id':
                    require(isinstance(node.values[i],ast.Constant) and node.values[i].value is None,'Inherited inactive scheduler metadata')
                    node.values[i]=ast.parse("self.executor_identity['job_id']",mode='eval').body; changed+=1
    require(changed==1,'One cell scheduler metadata substitution')
    ns = dict(vars(old)); exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])), __file__, 'exec'), ns)
    return types.SimpleNamespace(**ns)


def evaluator_class(old, radial, base, checks, release_gate):
    bounded = bounded_parent(old)
    facade = types.SimpleNamespace(**vars(base))
    facade.effective = lambda v, m, d: base.effective(v[:12], m, d)+[v[12 if m == 'La' else 13]]
    Parent = bounded.evaluator_class(facade)
    Radial = radial.evaluator_class(bounded, base, checks, release_gate)

    class Evaluator(Radial):
        def __init__(self, cfg, path, mode, inventory, parent_config, allocation):
            require(mode in ('prepare', 'run') and not Path(path).exists(), 'Fresh bounded campaign cache')
            self.batch = copy.deepcopy(cfg); self.mode = mode; self.batch_fingerprint = cfg['fingerprint']
            self.allowed = {}; self.complete_vectors = {}; self.loaded_saved_vectors = {}
            self.registry = {e['key']: e for e in inventory['entries']}; self.imported_vector_keys = set(self.registry)
            require(len(self.registry) == 31, 'Only exact31 admitted seed entries')
            self.old_config = copy.deepcopy(parent_config)
            self.executor_identity = dict(allocation, controller_pid=os.getpid(), run_id=cfg['run_id'], implementation=cfg['implementation'], adapter=cfg['adapter'], parent_implementation=cfg['parent_evaluator'])
            runtime = copy.deepcopy(parent_config)
            runtime.update(source_manifest=cfg['native_source_manifest'], executables=cfg['executables'], build_receipt=cfg['build'], qualification_review=cfg['qualification_review'])
            base.Evaluator.__init__(self, runtime, path, root_release=cfg.get('root_release'))
            self.initial_vector = list(self.registry[inventory['origin_key']]['parameters'])
            self.family['radial_extension'] = dict(append_order=['La_radial_drive', 'Dy_radial_drive'], bounds=[[-1., 1.], [-1., 1.]])
            self.program_attempts = []; self.cell_attempts = {}; self.started_cells = set(); self.imported_cell_keys = set(); self.cell_authority = {}
            self.saved_metal_index = {m: {} for m in ('La', 'Dy')}
            for key, entry in self.registry.items():
                for metal in self.saved_metal_index:
                    effective = base.digest(facade.effective(entry['parameters'], metal, self.domain))
                    self.saved_metal_index[metal].setdefault(effective, key)
            self.ledger.update(limits=copy.deepcopy(cfg['limits']), new_policy=cfg['adapter'], configuration_fingerprint=self.batch_fingerprint,
                               initial_imported_unique_cells=0, available_saved_vectors=31, native_returncodes_observed=0)
            self.persistence = old.J.Journal(self.cache_dir/'EVENTS.jsonl', self.ledger, self.lock,
                                            dict(implementation=cfg['adapter'], configuration_fingerprint=self.batch_fingerprint, root_release=cfg.get('root_release')))
            self.persistence.materialize(self.ledger_path, expected=self.ledger)

        def load_vector(self, key):
            entry = self.registry[key]; v = base.read(entry['artifact'])
            require(v['status'] == 'complete', 'Complete saved seed')
            if len(v['parameters']) == 12:
                require(entry['lineage'] == 'old_zero_view' and v['parameters']+[0., 0.] == entry['parameters'], 'Explicit zero bridge only')
                v = dict(v, parameters=entry['parameters'], radial_zero_view=dict(original_artifact=entry['artifact'], root_cache_bridge=self.batch['root_cache_bridge'], original_producer_retained=True, new_native_calls=0))
            require(v['parameters'] == entry['parameters'] and base.digest(v['parameters']) == key and set(v['cells']) == set(self.m['rows']), 'Exact saved31 vector')
            return v

        def _import_saved(self, key, value=None):
            with self.lock:
                if key in self.loaded_saved_vectors: return self.loaded_saved_vectors[key]
                v = self.load_vector(key) if value is None else copy.deepcopy(value)
                for cid, record in v['cells'].items():
                    require(record['status'] == 'complete' and record['cache_key'] == base.digest(record['cache_identity']), 'Immutable complete cell identity')
                    token, identity, rendered, meta = self.cache_identity(v['parameters'], cid)
                    row = self.m['rows'][cid]; actual = record['cache_identity']
                    old_key, old_meta = base.render(v['parameters'][:12], row, self.domain, self.h['resolver'])
                    old_binary = actual['source_manifest'] == self.old_config['source_manifest']
                    if old_binary:
                        require(v['parameters'][12 if row['metal'] == 'La' else 13] == 0., 'No nonzero cross-binary reuse')
                        expected = dict(cell_id=cid, effective_values=old_meta['physical_values_a_b_c_d_R_epsilon'], xyz=row['source_input']['xyz'], parameters=row['source_input']['parameters'],
                                        key_sha256=hashlib.sha256(old_key.encode()).hexdigest(), implementation=self.old_config['implementation'], source_manifest=self.old_config['source_manifest'],
                                        executables=self.old_config['executables'], helpers={k:self.old_config[k] for k in base.HELPERS})
                        require(actual == expected, 'Exact old scientific cell under admitted zero bridge')
                    else:
                        expected = dict(identity, radial_adapter=actual.get('radial_adapter'))
                        require(actual.get('radial_adapter') in self.batch['admitted_radial_adapters'] and actual == expected, 'Exact qualified radial scientific cell; policy metadata only differs')
                    for program, inp in record['program_inputs'].items():
                        text = base.render(v['parameters'][:12], row, self.domain, self.h['resolver'], row['needs_export'] if program == 'analyze' else False)[0] if old_binary else self.render(v['parameters'], cid, row['needs_export'] if program == 'analyze' else False)[0]
                        require(base.verify(inp['key']).read_text() == text, 'Exact inherited rendered native key')
                        for name in ('xyz', 'parameters'):
                            require(base.verify(inp[name]).read_bytes() == base.verify(row['source_input'][name]).read_bytes(), 'Literal inherited geometry/parameters')
                    for receipt in record['receipts']:
                        base.verify(receipt['output']); require(receipt['returncode'] == 0, 'Successful inherited receipt')
                    if token not in self.cache:
                        self.cache[token] = copy.deepcopy(record); self.cell_authority[token] = dict(artifact=self.registry[key]['artifact'], key=key, old_zero_bridge=old_binary)
                    self.imported_cell_keys.add(token)
                self.loaded_saved_vectors[key] = v
                return v

        def _saved_metal(self, values, metal):
            return self.saved_metal_index[metal].get(base.digest(facade.effective(values, metal, self.domain)))

        def authorize(self, values, role, step):
            values = self.checked_values(values); key = base.digest(values)
            require(role in ('probe', 'trial') and 1 <= step <= 8, 'Finite explicit driver role')
            self.allowed[key] = dict(parameters=values, role=role, outer_step=step)
            return key

        def required_programs(self, values):
            values = self.checked_values(values); total = 0
            for cid, row in self.m['rows'].items():
                token = self.cache_identity(values, cid)[0]
                if token not in self.cache and self._saved_metal(values, row['metal']) is None:
                    total += 1+int(row['needs_gradient'])
            return total

        def _cell(self, values, cid):
            token = self.cache_identity(values, cid)[0]
            saved = self._saved_metal(values, self.m['rows'][cid]['metal'])
            if token not in self.cache and saved is not None: self._import_saved(saved)
            require(token in self.cache or facade.effective(values, self.m['rows'][cid]['metal'], self.domain) != facade.effective(self.initial_vector, self.m['rows'][cid]['metal'], self.domain), 'No unqualified baseline fallback')
            answer = Parent._cell(self, values, cid)
            if token in self.imported_cell_keys:
                answer['immutable_campaign_reuse'] = dict(self.cell_authority[token], root_cache_bridge=self.batch['root_cache_bridge'], requested_vector=values,
                                                        requested_identity=self.cache_identity(values, cid)[1], original_producer_and_raw_cell_retained=True, new_native_calls=0)
            return answer

        def _gate(self):
            require(self.mode == 'run' and not self.stop_event.is_set(), 'No native in preparation/after stop')
            require(self.ledger['native_program_invocations'] < 7280, 'Hard7280 invocation ceiling')
            require(all(os.environ.get(k) == '1' for k in self.batch['threads']) and os.environ.get('SLURM_JOB_ID') == self.executor_identity['job_id'] and os.environ.get('SLURM_CPUS_PER_TASK') == '13' and os.environ.get('SLURM_NTASKS') == '1', 'Truthful13-CPU Slurm/12 workers/one-thread processes')
            if not self._gated: release_gate(self.batch); self._gated = True

        def evaluate(self, values):
            values = self.checked_values(values)
            require(base.digest(values) in self.allowed, 'Explicit authorized vector role')
            # Validate every claimed saved-metal identity before any worker can
            # consume the corresponding reduced whole-vector native budget.
            for metal in ('La', 'Dy'):
                saved = self._saved_metal(values, metal)
                if saved is not None: self._import_saved(saved)
            return Parent.evaluate(self, values)

        def snapshot(self):
            with self.lock: self.persistence.materialize(self.ledger_path, expected=self.ledger)
            s = base.Evaluator.snapshot(self)
            s.update(active_limits=self.batch['limits'], new_parameter_vectors=len(self.ledger['vector_keys']), native_returncodes_observed=self.ledger['native_returncodes_observed'],
                     imported_unique_cells=len(self.imported_cell_keys), available_saved_vectors=31, loaded_saved_vectors=len(self.loaded_saved_vectors), old_mutable_ledger_read=False, executor=copy.deepcopy(self.executor_identity))
            return s
    return Evaluator
