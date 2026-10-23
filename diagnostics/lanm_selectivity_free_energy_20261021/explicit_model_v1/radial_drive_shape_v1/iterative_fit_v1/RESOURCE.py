"""Scheduling-only13-CPU allocation binding; no scientific operations."""
import os
import socket
import subprocess


def require(ok, message):
    if not ok: raise ValueError(message)


def parse(text, env):
    fields=dict(token.split('=',1) for token in text.split() if '=' in token)
    job=env.get('SLURM_JOB_ID','')
    require(job.isdigit() and fields.get('JobId')==job, 'Exact real Slurm job')
    require(fields.get('JobState')=='RUNNING' and fields.get('Account')=='banfield' and fields.get('Partition')=='standard' and fields.get('QOS')=='standard_expand', 'Running declared allocation')
    require(fields.get('NumNodes')=='1' and fields.get('NumTasks')=='1' and fields.get('NumCPUs')=='13' and fields.get('CPUs/Task')=='13' and fields.get('Requeue')=='0', 'One nonrequeued13-CPU task')
    require(env.get('SLURM_NTASKS')=='1' and env.get('SLURM_JOB_NUM_NODES',env.get('SLURM_NNODES'))=='1' and env.get('SLURM_CPUS_PER_TASK')=='13', 'Actual native/controller envelope')
    require(env.get('SLURM_RESTART_COUNT','0')=='0', 'No automatic restart')
    return dict(mode='slurm',job_id=job,hostname=socket.gethostname(),allocated_CPUs=13,nodes=1,tasks=1,native_workers=12,controller_threads=1,
                scheduler_fields=fields,default_memory_and_time=True,scheduler_allocated_CPU_seconds=None,terminal_accounting_pending=True)


def allocation(command):
    job=os.environ.get('SLURM_JOB_ID',''); require(job.isdigit(),'Real Slurm allocation required')
    result=subprocess.run([command,'show','job','-o',job],capture_output=True,text=True,check=True,timeout=30)
    require(not result.stderr.strip(),'Clean scheduler readback')
    return parse(result.stdout,os.environ)
