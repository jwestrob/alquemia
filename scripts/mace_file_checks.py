"""Operation-local digest reuse; every lookup still checks file identity/change time.

Only callers decorated here opt in. No scientific result/parsed record is cached,
and no shared file or persistent cache is changed. Immutable artifact DAGs often
verify the same checkpoint hundreds of times within one prepare/report operation.
"""
from contextvars import ContextVar
from functools import wraps
from pathlib import Path
import os
import time
import affordable_common as common

_CACHE=ContextVar('mace_operation_file_digests',default=None)


def identity(path):
    # Opening also requests NFS close-to-open attribute revalidation.
    with path.open('rb') as handle:s=os.fstat(handle.fileno())
    return (s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns)


def cached_file_checks(function):
    @wraps(function)
    def run(*args,**kwargs):
        if _CACHE.get() is not None:return function(*args,**kwargs)
        original=common.digest
        token=_CACHE.set({})

        def checked_digest(path):
            cache=_CACHE.get()
            if cache is None:return original(path)
            p=Path(path);before=identity(p);key=(str(p.absolute()),before)
            # Some filesystems coarsen timestamps. Never reuse a digest for a
            # file modified within that window, including same-tick rewrites.
            stable=max(before[-2:]) < time.time_ns()-2_000_000_000
            value=cache.get(key) if stable else None
            if value is None:value=original(p)
            if identity(p)!=before:
                raise common.InvalidArtifact(f'artifact changed during verification: {p}')
            if stable:cache[key]=value
            return value

        common.digest=checked_digest
        try:return function(*args,**kwargs)
        finally:
            common.digest=original
            _CACHE.reset(token)
    return run
