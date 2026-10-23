"""Two fixed56-program vectors; unchanged radial renderer/parser/journal/invoke."""
import ast,copy,types
from pathlib import Path

def factory(parent):
    tree=ast.parse(Path(parent.__file__).read_text());fn=copy.deepcopy(next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='evaluator_class'))
    counts={28:0,4:0}
    for n in ast.walk(fn):
        if isinstance(n,ast.Constant) and type(n.value) is int and n.value in counts:
            old=n.value;counts[old]+=1;n.value={28:56,4:2}[old]
        elif isinstance(n,ast.Constant) and n.value=='radial_probe':n.value='fixed_native_comparison'
    if counts!={28:4,4:1}:raise ValueError('Unexpected parent policy constants')
    ns=dict(vars(parent));exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),__file__,'exec'),ns)
    return ns['evaluator_class']
