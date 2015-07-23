##################
#  Post Process  #
##################

import os
import subprocess
from subprocess import Popen,PIPE

def Post(nodes):
    counter = 1

    # Set the time parameters
    DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME      = 0.5
    DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT = 0.05

    # Node number for data extraction
    for node in nodes:
        # Create the Result file to store the flow and area versus time for each node
        createFile = open("Results/node_"+str(counter),'w+')
        # Loop through the output files
        for x in range(0,int(DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME/DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT),1):
            outputFile = open("output/Solid/Solid0000"+str(x)+".part0.exnode")
            outputLINE = outputFile.readlines()
            for i in range(0,len(outputLINE)):
                if outputLINE[i].split() == ['Node:', str(node)]:
                    # Extract the variables from output files
                    X1 = float(''.join(outputLINE[i+1].split()))
                    Y1 = float(''.join(outputLINE[i+2].split()))
                    X2 = float(''.join(outputLINE[i+3].split()))
                    Y2 = float(''.join(outputLINE[i+4].split()))
                    Time = x*DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT
                    # Write in the Result file
                    Result = open("Results/node_"+str(counter),'a+')
                    print >> Result,"%.4f"%Time,X1,Y1,X2,Y2
        counter = counter+1
    return

nodes = input('Nodes: ')
Post(nodes)

print "."
print "."
print "."
print "Processing Completed!"

Popen(['gnuplot','PlotNode.p'])

