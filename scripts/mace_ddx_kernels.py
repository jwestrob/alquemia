"""Prepare an isolated thin binding to unchanged native ddX PCM operators."""
from pathlib import Path
import argparse
import difflib
import shutil
import tarfile
from affordable_common import InvalidArtifact,read_json,record,verify,write_new

FORTRAN='''
! Alquemia diagnostic adapter. All operator algorithms remain upstream.
subroutine alquemia_pcm_apply(c_ddx, operation, nbasis, nsph, x, y, c_error) bind(C)
    type(c_ptr), intent(in), value :: c_ddx, c_error
    integer(c_int), intent(in), value :: operation, nbasis, nsph
    real(c_double), intent(in) :: x(nbasis, nsph)
    real(c_double), intent(out) :: y(nbasis, nsph)
    type(ddx_setup_type), pointer :: model
    type(ddx_error_type), pointer :: error
    logical :: original_diagonal
    call c_f_pointer(c_ddx, model)
    call c_f_pointer(c_error, error)
    original_diagonal = model%constants%dodiag
    model%constants%dodiag = .true.
    select case(operation)
    case(1)
        call rinfx(model%params,model%constants,model%workspace,x,y,error)
    case(2)
        call repsx(model%params,model%constants,model%workspace,x,y,error)
    case(3)
        call lx(model%params,model%constants,model%workspace,x,y,error)
    case(4)
        call prec_repsx(model%params,model%constants,model%workspace,x,y,error)
    case(5)
        call ldm1x(model%params,model%constants,model%workspace,x,y,error)
    case default
        call update_error(error,'Unsupported diagnostic PCM operator')
    end select
    model%constants%dodiag = original_diagonal
end subroutine

subroutine alquemia_pcm_source(c_state, nbasis, nsph, phi) bind(C)
    type(c_ptr), intent(in), value :: c_state
    integer(c_int), intent(in), value :: nbasis, nsph
    real(c_double), intent(out) :: phi(nbasis, nsph)
    type(ddx_state_type), pointer :: state
    call c_f_pointer(c_state, state)
    phi = state%phi
end subroutine

function alquemia_pcm_norm(c_ddx, nbasis, nsph, x) result(value) bind(C)
    type(c_ptr), intent(in), value :: c_ddx
    integer(c_int), intent(in), value :: nbasis, nsph
    real(c_double), intent(in) :: x(nbasis, nsph)
    real(c_double) :: value
    type(ddx_setup_type), pointer :: model
    call c_f_pointer(c_ddx, model)
    value = hnorm(model%params%lmax,nbasis,nsph,x)
end function
'''
HEADER='''
/* Alquemia: thin diagnostic access, not a replacement PCM implementation. */
void alquemia_pcm_apply(void* model, int operation, int nbasis, int nsph,
                       const double* x, double* y, void* error);
void alquemia_pcm_source(void* state, int nbasis, int nsph, double* phi);
double alquemia_pcm_norm(void* model, int nbasis, int nsph, const double* x);
'''
CPP='''
        .def("pcm_source_vector", [](std::shared_ptr<State> self) {
              auto model = self->model();
              if (model->model() != "pcm") throw py::value_error("PCM required");
              array_f_t result({model->n_basis(), model->n_spheres()});
              alquemia_pcm_source(self->holder(),model->n_basis(),model->n_spheres(),result.mutable_data());
              return result;
             }, "Native signed source vector assembled by pcm_setup.")
        .def("pcm_operator", [](std::shared_ptr<State> self, int operation, array_f_t x) {
              auto model = self->model();
              if (model->model() != "pcm" || operation < 1 || operation > 5 ||
                  x.ndim() != 2 || x.shape(0) != model->n_basis() || x.shape(1) != model->n_spheres())
                  throw py::value_error("PCM operator or array dimensions differ");
              array_f_t result({model->n_basis(),model->n_spheres()});
              alquemia_pcm_apply(model->holder(),operation,model->n_basis(),model->n_spheres(),x.data(),result.mutable_data(),model->error());
              self->throw_if_error();
              return result;
             }, "operation"_a, "x"_a,
             "Native full operators: 1 R_infinity, 2 R_epsilon, 3 L; native inverse-block preconditioners: 4 R, 5 L.")
        .def("pcm_native_norm", [](std::shared_ptr<State> self, array_f_t x) {
              auto model = self->model();
              if (model->model() != "pcm" || x.ndim() != 2 || x.shape(0) != model->n_basis() || x.shape(1) != model->n_spheres())
                  throw py::value_error("PCM norm array dimensions differ");
              return alquemia_pcm_norm(model->holder(),model->n_basis(),model->n_spheres(),x.data());
             }, "x"_a, "Upstream H(-1/2) norm; no iterative solve.")
'''


def prepare(parent_manifest,plan,output):
    parent=read_json(parent_manifest)
    if parent['protocol']!='isolated_pyddx_0p9p0_portable_cmake_build_v2':raise InvalidArtifact('pinned parent build required')
    for pin in parent['pins']:verify(pin)
    archive=next(pin for pin in parent['pins'] if pin['path']==parent['source_archive'])
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);source=root/'source';source.mkdir()
    with tarfile.open(verify(archive)) as f:f.extractall(source,filter='data')
    tree=source/'pyddx-0.9.0';patches=[];source_pins={}
    for name in ('ddx_cinterface.f90','ddx.h','pyddx_classes.cpp'):
        p=tree/'src'/name;before=p.read_text();original=record(p)
        if name=='ddx_cinterface.f90':
            if not before.endswith('\nend\n'):raise InvalidArtifact('unexpected native module terminator')
            after=before[:-4]+FORTRAN+'\nend\n'
        elif name=='ddx.h':
            marker='#ifdef __cplusplus\n}\n#endif'
            if before.count(marker)!=1:raise InvalidArtifact('C linkage terminator differs')
            after=before.replace(marker,HEADER+'\n'+marker)
        else:
            marker='        .def_property_readonly("model", &State::model, "Model definition")'
            if before.count(marker)!=1:raise InvalidArtifact('State binding anchor differs')
            after=before.replace(marker,CPP+marker)
            marker='  py::class_<State, std::shared_ptr<State>>('
            if after.count(marker)!=1:raise InvalidArtifact('State declaration anchor differs')
            after=after.replace(marker,'  m.attr("alquemia_pcm_kernel_version") = "native_pcm_operators_v1";\n'+marker)
        p.write_text(after);source_pins[name]=dict(original_sha256=original['sha256'],patched=record(p))
        patches.extend(difflib.unified_diff(before.splitlines(keepends=True),after.splitlines(keepends=True),fromfile='upstream/'+name,tofile='adapter/'+name))
    (root/'native_operator_adapter.patch').write_text(''.join(patches));shutil.copyfile(__file__,root/Path(__file__).name)
    unchanged={}
    # Every numerical Fortran module except the appended C adapter remains unchanged.
    with tarfile.open(verify(archive)) as f:
        for p in sorted((tree/'src').glob('*.f90')):
            if p.name=='ddx_cinterface.f90':continue
            original=f.extractfile('pyddx-0.9.0/src/'+p.name).read()
            if original!=p.read_bytes():raise InvalidArtifact('native numerical source changed: '+p.name)
            unchanged[p.name]=record(p)
    patched=root/'pyddx-0.9.0-alquemia-operators.tar.gz'
    with tarfile.open(patched,'w:gz') as f:f.add(tree,arcname='pyddx-0.9.0')
    wrapper=Path(__file__).resolve().parents[1]/'diagnostics/mace_omol_20260917/run_ddx_build_v2.sbatch'
    m={**parent,'protocol':'isolated_pyddx_native_operator_adapter_build_v1','source_archive':str(patched),
       'parent_build_manifest':record(parent_manifest),'plan':record(plan),'builder':record(root/Path(__file__).name),
       'patch':record(root/'native_operator_adapter.patch'),'adapter_sources':source_pins,'unchanged_native_modules':unchanged,
       'wrapper':record(wrapper),'scientific_calculations':0}
    m['pins']=[*parent['pins'],record(patched),m['plan'],m['builder'],m['patch'],m['wrapper'],*[v['patched'] for v in source_pins.values()],*unchanged.values()]
    write_new(root/'manifest.json',m);return m


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--parent-build-manifest',required=True);p.add_argument('--plan',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();m=prepare(a.parent_build_manifest,a.plan,a.output)
    print({'protocol':m['protocol'],'native_modules_unchanged':len(m['unchanged_native_modules']),'scientific_calculations':0})
