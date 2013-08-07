##################
#  Post Process  #
##################

import os
import subprocess
from subprocess import Popen, PIPE

def Post(nodes):
    # Set the reference values
    Ts = 0.001                   # Time     (s)
    Qs = 100.0e-6                # Flow     (m3/s)
    counter   = 1

    # Set the time parameters
    DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME      = 900.0
    DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT = 1.0

    # Node number for data extraction
    for node in nodes:
        # Create the Result file to store the flow and area versus time for each node
        createFile = open("Results/node_"+str(counter),'w+')
        # Loop through the output files
        for x in range(0,int(DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME/DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT)):
            outputFile = open("output/MainTime_"+str(x)+".part0.exnode")
            outputLINE = outputFile.readlines()
            for i in range(0,len(outputLINE)):
                if outputLINE[i].split() == ['Node:', str(node)]:
                    # Extract the variables from output files
                    Flow = float(''.join(outputLINE[i+4].split()))
                    Pressure = float(''.join(outputLINE[i+11].split()))
                    Time = x*Ts*DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT
                    Flow = Flow*Qs*1000000.0
                    # Write in the Result file
                    Result = open("Results/node_"+str(counter),'a+')
                    print >> Result,"%.4f"%Time,Flow,Pressure
        counter = counter+1
    return

nodes = input('Nodes: ')
Post(nodes)

print "."
print "."
print "."
print "Processing Completed!"

Popen(['gnuplot','Plot.p'])

