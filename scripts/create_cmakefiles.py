#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import collections
import copy
from sets import Set
import os
from shutil import copyfile

cmake_content = """
# Add example executable
add_executable({executable_name} {src})

# Turn on Fortran preprocessing (#include directives)
if (MSVC)
    set(PROC_OPT "/fpp")
else()
    set(PROC_OPT "-cpp")
endif()

# in debug mode, add flag that make Floating point exceptions fatal
set(OPT "")
if (FALSE)
    if(CMAKE_BUILD_TYPE MATCHES DEBUG)
        set(OPT "-ffpe-trap=invalid,zero,underflow,overflow,denormal")		# GCC: make floating point exceptions fatal
    endif()
endif()

target_compile_options({executable_name} PRIVATE ${PROC_OPT} ${OPT} "-fbacktrace")

# Link to opencmiss
# -----------------
#
# This simply uses the opencmiss target created by OC_INIT() above.
#
# Alternatively, you can also directly invoke 'find_package(Iron <IRON_VERSION>)' to explicitly
# require that version. But then you have to deal with setting up the correct MPI include directives and library paths.
#
# If required, add any other required link libraries (cellml, petsc ..) here, too.
# For example, if you needed PetSC functionality, issue
#
#     find_package(PETSC <PETSC_VERSION> REQUIRED)
#     target_link_libraries(${EXAMPLE_BINARY} PRIVATE petsc)
#
# All the OpenCMISS dependencies provide a target you can link against corresponding to the (lowercase) component name.
target_link_libraries({executable_name} PRIVATE opencmiss)

# Add a simple test that runs the executable
add_test(NAME test_{executable_name} COMMAND {executable_name})
add_opencmiss_environment(test_{executable_name})

###################
# Developer notice!
#
# If you write Fortran code and use MPI, you need to use the following MPI directives:
#
# #ifndef NOMPIMOD
#   USE MPI
# #endif
# [...]
#   IMPLICIT NONE
# [...]
# #ifdef NOMPIMOD
#     #include "mpif.h"
# #endif
#
# Reasoning: In some cases like Windows/MPICH2 there sometimes is no mpi.mod file. In order to yet make
# the example work the build system adds the definition 'NOMPIMOD', which can be checked and acted to accordingly.
#
"""


for (dirpath, subdirectories, filenames) in os.walk(".."):
  print "directory {}".format(dirpath)
  src_found = False
  if "src" in subdirectories:
    directory_content = os.listdir(os.path.join(dirpath, "src"))
    if os.path.basename(dirpath) != "Fortran":
      #os.mkdir("Fortran")
      print 'rename "{}" to "{}"'.format(os.path.join(dirpath,"src"), os.path.join(dirpath,"Fortran","src"))
      try:
        os.mkdir(os.path.join(dirpath,"Fortran"))
      except:
        pass
      os.rename(os.path.join(dirpath,"src"), os.path.join(dirpath,"Fortran","src"))
      src_found = True
    
  if "Fortran" in subdirectories:
    src_found = True
      
  if src_found:
    print "process directory {}".format(dirpath)
    
    
    copyfile("resources/CMakeLists-1.txt", os.path.join(dirpath,"CMakeLists.txt"))
    copyfile("resources/OpenCMISS.cmake", os.path.join(dirpath,"OpenCMISS.cmake"))
    
    src_files = os.listdir(os.path.join(dirpath, "Fortran","src"))
    src_string = ""
    for filename in src_files:
      src_string += " src/"+filename
    executable_name = os.path.basename(dirpath)
    
    with open(os.path.join(dirpath,"Fortran","CMakeLists.txt"),'w') as file:
      file.write(cmake_content.format(src=src_string, executable_name=executable_name, PROC_OPT="{PROC_OPT}", OPT="{OPT}", EXAMPLE_BINARY="{EXAMPLE_BINARY}"))
      
    try:
      os.mkdir(os.path.join(dirpath,"build_debug"))
    except:
      pass
    copyfile("resources/build.sh", os.path.join(dirpath,"build_debug","build.sh"))
