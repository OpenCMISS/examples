#!/usr/bin/env python

import os
import sys
import shutil
import re

rootDir = os.getcwd()
template = os.path.join(rootDir, "CMakeLists.template.cmake")

re_findMPI = re.compile(r"(\s)*USE MPI\n((.|\n)*)IMPLICIT NONE\n", re.MULTILINE)

for dirName, subdirList, fileList in os.walk(rootDir):
    if dirName != rootDir:
        for file in fileList:
            #if file == 'Makefile':
            #    print "Creating CMakeLists.txt in " + dirName
            #    shutil.copy(template, os.path.join(dirName, "CMakeLists.txt"))
            #    continue
            if re.match(".*\.f(90)?$",file):
                thefile = os.path.join(dirName,file)
                mf = open(thefile, 'r')
                source = mf.read()
                mf.close()
                match = re_findMPI.search(source)
                print "Checking %s" % file
                if match:
                    repl = r'\n#ifndef NOMPIMOD\n\1 USE MPI\n#endif\n\2IMPLICIT NONE\n\n#ifdef NOMPIMOD\n#include "mpif.h"\n#endif\n\n'
                    source = re_findMPI.sub(repl,source,1)
                    mf = open(thefile, 'w')
                    mf.write(source)
                    mf.close()
                    print "Processed %s in %s" % (file, dirName) 
