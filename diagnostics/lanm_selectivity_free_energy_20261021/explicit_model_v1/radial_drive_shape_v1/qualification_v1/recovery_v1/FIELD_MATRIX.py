"""One scoped header predicate; inherited matrix arithmetic is unchanged."""
import ast,copy,math,types
from pathlib import Path

ZERO_GUARD="require(lines[4][0]=='EXTERNAL' and len(lines[4])==5 and lines[4][1]=='F' and np.array_equal(numeric(lines[4][2:]),np.zeros(3)),'Zero-field export required')"

def header(tokens,row,units,job,matrix):
    require=matrix.require
    require(job['stage']=='fields' and job['program']=='analyze' and job['cell'] in ('water_central_La','water_central_Dy') and job['sign'] in (-1,1) and job['field_au']==job['sign']*.001 and not job['fd'],'Exact declared field slot')
    require(job['id']=='__'.join(('fields',job['cell'],'analyze',str(job['sign']))),'Exact field identity')
    for k in ('xyz','key','parameters'):
        require(matrix.pin(matrix.verify(row['input'][k]))['sha256']==matrix.pin(matrix.verify(job['input'][k]))['sha256'],'Exact saved/prepared field inputs')
    key=matrix.verify(row['input']['key']).read_text()
    records=[s.split('#',1)[0].split() for s in key.splitlines()]
    fields=[r for r in records if r and r[0].upper()=='EXTERNAL-FIELD']
    expected_key=['EXTERNAL-FIELD','0','0',format(job['field_au']*units['MV_cm_per_field_au'],'.17g')]
    require(fields==[expected_key] and not any(r and r[0].upper()=='EXFLD-FREQ' for r in records),'Exact static external field key')
    require(units['elefield']==1439.96455,'Pinned native elefield conversion')
    expected=[float(x)/units['elefield'] for x in expected_key[1:]]
    require(len(tokens)==5 and tokens[:2]==['EXTERNAL','T'],'Nonzero field export grammar/flag')
    actual=matrix.numeric(tokens[2:]).tolist();bound=matrix.arithmetic(3*row['atoms'],expected,actual)
    require(actual[:2]==[0.,0.] and actual[2]*job['sign']>0 and max(abs(a-b) for a,b in zip(actual,expected))<=bound,'Native field conversion/readback')
    return {'native_external_tokens':tokens,'key_tokens':expected_key,'expected_texfld_e_per_A2':expected,'actual_texfld_e_per_A2':actual,'maximum_error':max(abs(a-b) for a,b in zip(actual,expected)),'arithmetic_bound':bound}

def scoped(matrix,job):
    # Copy one function AST and replace exactly one complete require statement.
    tree=ast.parse(Path(matrix.__file__).read_text());fn=copy.deepcopy(next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='parse_export'))
    expected=ast.dump(ast.parse(ZERO_GUARD).body[0],include_attributes=False);hits=[i for i,n in enumerate(fn.body) if ast.dump(n,include_attributes=False)==expected]
    matrix.require(len(hits)==1,'Exact inherited zero-header predicate')
    fn.body[hits[0]]=ast.parse('check_external(lines[4],row,units)').body[0]
    metadata=[]
    def check(tokens,row,units):metadata.append(header(tokens,row,units,job,matrix))
    ns=dict(vars(matrix));ns['check_external']=check
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),str(Path(__file__)),'exec'),ns)
    parsed=ns['parse_export']
    def parse(path,row,units):
        before=len(metadata);result=parsed(path,row,units);matrix.require(len(metadata)==before+1,'One external record per export');result['external_field']=metadata[-1];return result
    ns['parse_export']=parse
    qualify=types.FunctionType(matrix.qualify_exports.__code__,ns,matrix.qualify_exports.__name__,matrix.qualify_exports.__defaults__,matrix.qualify_exports.__closure__)
    proxy=types.SimpleNamespace(**vars(matrix));proxy.parse_export=parse;proxy.qualify_exports=qualify
    return proxy
