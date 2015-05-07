__author__ = 'zwan145'

import sys, os
import numpy
from utils import *

## This script implements mesh convergence for the cantilever problem.
# It modifies the problem1_tricubic.py script for each mesh resolution
# and calls the python script for each mesh resolution.
# It then reads out the 2PK fibre stress for the 36 data point locations, and
# writes out the vectors themselves and the mesh topology in a text file.
###############################################################################
#"""
# Set up log file
so = se = open("Output.log", 'w', 0)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
os.dup2(so.fileno(), sys.stdout.fileno())
os.dup2(se.fileno(), sys.stderr.fileno())
#"""

# Mesh Refinements
circ_elems = [4,8,16,32]
trans_elems = [2,4,8,16]
# Read in single element script, change number of elements and the boundary conditions
try:
    fidr = open('LVInflation.py', 'r')
except IOError:
    print 'ERROR: Unable to open LVInflation.py'
    quit()

for i in range(0, len(trans_elems)):
    for j in range(0, len(circ_elems)):
        elems = [circ_elems[j], circ_elems[j], trans_elems[i]]
        print '++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
        print '    Solving for '+str(elems)+' elements\n'
        print '++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'

        # Rewrite python file
        filenamew = 'LVInflation_'+str(elems[0])+'-'+str(elems[1])+'-'+str(elems[2])+'.py'
        print filenamew

        RewritePythonScript('LVInflation.py', filenamew, elems)

        # Run new python script
        os.system('python '+filenamew)

        print '++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
        print '    Finished solving for '+str(elems)+' elements\n'
        print '++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'

        cmd = "python ConvertToVTK.py '"+str(elems[2])+"' '"+str(elems[1])+"' '"+str(elems[0])+"'\n"
        os.system(cmd)
