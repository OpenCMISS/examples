#@ shell = /bin/sh
#@ job_name = classicalfield_laplace_laplace_fortran
#@ class = default
#@ group = nesi
#@ account_no = uoa
#@ notification = never
#@ account_no = /nz/nesi
#@ wall_clock_limit = 1:00:00
#@ resources = ConsumableMemory(1024mb) ConsumableVirtualMemory(1024mb)
#@ job_type = MPICH
#@ total_tasks = 2
#@ blocking=unlimited
#@ notification = complete
#@ notify_user = ting.yu@auckland.ac.nz
#@ initialdir = /home/tyu011/opencmiss/OpenCMISS/examples/ClassicalField/Laplace/Laplace/Fortran
#@ output = /home/tyu011/opencmiss/OpenCMISS/logs/examples/ClassicalField/Laplace/Laplace/Fortran/joblogs/$(job_name).$(jobid).out
#@ error = /home/tyu011/opencmiss/OpenCMISS/logs/examples/ClassicalField/Laplace/Laplace/Fortran/joblogs/$(job_name).$(jobid).err
#@ queue
 
mpirun  bin/x86_64-linux/openmpi/gnu_4.4/FortranExample-debug
