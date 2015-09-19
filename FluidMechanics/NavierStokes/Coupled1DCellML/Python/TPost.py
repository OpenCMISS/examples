##################
#  Post Process  #
##################

import os
import subprocess
from subprocess import Popen, PIPE

def Post(times):
    counter = 1

    # Node number for data extraction
    for time in times:
        # Create the Result file to store the flow and area versus time for each node
        createFile = open("Results/time_"+str(counter),'w+')
        # Loop through the output files
        outputFile = open("output/MainTime_"+str(time)+".part0.exnode")
        outputLINE = outputFile.readlines()
        for node in range(1,len(outputLINE)):
            for i in range(0,len(outputLINE)):
                if outputLINE[i].split() == ['Node:', str(node)]:
                    # Extract the variables from output files
                    Flow = float(''.join(outputLINE[i+4].split()[0]))
                    Pressure = float(''.join(outputLINE[i+12].split()[0]))
                    Conc = float(''.join(outputLINE[i+14].split()[0]))
                    X = float(''.join(outputLINE[i+1].split()[0]))
                    Y = float(''.join(outputLINE[i+2].split()[0]))
                    Z = float(''.join(outputLINE[i+3].split()[0]))
                    L = (X**2.0+Y**2.0+Z**2.0)**0.5
                    # Write in the Result file
                    Result = open("Results/time_"+str(counter),'a+')
                    print >> Result,"%.4f"%L,Flow,Pressure,Conc
        counter=counter+1
    return

times = input('Times: ')
Post(times)

print "."
print "."
print "."
print "Processing Completed!"

Popen(['gnuplot','PlotTime.p'])

