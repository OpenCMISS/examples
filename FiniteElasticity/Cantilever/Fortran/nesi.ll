#@ shell = /bin/sh
#@ job_name = nesi
#@ class = default
#@ group = nesi
#@ account_no = uoa
#@ notification = never
#@ wall_clock_limit = 1:00:00
#@ resources = ConsumableMemory(1024mb) ConsumableVirtualMemory(1024mb)
#@ job_type = MPICH
#@ total_tasks = 2
#@ blocking=unlimited
#@ initialdir = /home/tyu011/opencmiss/OpenCMISS/examples/FiniteElasticity/Cantilever/Fortran
#@ output = $(job_name).$(jobid).out
#@ error = $(job_name).$(jobid).out
#@ queue
 
mpirun  bin/x86_64-linux/openmpi/gnu_4.4/FortranExample-debug
