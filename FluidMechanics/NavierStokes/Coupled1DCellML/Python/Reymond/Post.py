##################
#  Post Process  #
##################

import os
import numpy as np
import subprocess
import pylab
from subprocess import Popen, PIPE
import FluidExamples1DUtilities as Utilities

def Post(nodes):
    # Set the reference values
    Ts = 1.0#0.001                   # Time     (s)
    Qs = 100.0e-8                # Flow     (m3/s)
    nodeCounter   = 0

    # Set the time parameters
    solveTime      = 1345.#790.0
    timeIncrement = 0.05
    outputFrequency = 10

    # Set up numpy arrays
    numberOfTimesteps = (solveTime/timeIncrement)/outputFrequency
    timeResult = np.zeros((numberOfTimesteps))
    numberOfNodes = len(nodes)
    flowResult = np.zeros((numberOfNodes,numberOfTimesteps))
    pressureResult = np.zeros((numberOfNodes,numberOfTimesteps))

    # Read in node info
    inputNodeNumbers        = []
    bifurcationNodeNumbers  = []
    trifurcationNodeNumbers = []
    coupledNodeNumbers      = []
    arteryLabels            = []
    filename = 'input/Node.csv'
    totalNumberOfNodes = Utilities.GetNumberOfNodes(filename)
    nodeCoordinates = np.zeros([totalNumberOfNodes,4,3])
    Utilities.CsvNodeReader(filename,inputNodeNumbers,bifurcationNodeNumbers,trifurcationNodeNumbers,coupledNodeNumbers,nodeCoordinates,arteryLabels)        

    # Set up reference types
    refs = ['Reymond2009','Reymond2011']
    refTypes = ['Experimental','Model']

    # Node number for data extraction
    nodeLabels = []
    # TODO: should swap these loops around to speed things up for multiple nodes (currently re-reads files for each node)
    for node in nodes:
        # Get artery names
        arteryLabel = arteryLabels[node-1]
        print('Reading results for ' + arteryLabel)
        nodeLabels.append(arteryLabel)
        # Create the Result file to store the flow and area versus time for each node
        createFile = open("results/node_"+str(node),'w+')
        # Loop through the output files
        timestep = 0
        for x in range(0,int(solveTime/timeIncrement),outputFrequency):
            outputFile = open("output/MainTime_"+str(x)+".part0.exnode")
            outputLINE = outputFile.readlines()
            for i in range(0,len(outputLINE)):
                if outputLINE[i].split() == ['Node:', str(node)]:
                    # Extract the variables from output files
                    Flow = float(''.join(outputLINE[i+4].split()))
                    Pressure = float(''.join(outputLINE[i+12].split()))
                    Time = x*Ts*timeIncrement
                    Flow = Flow*Qs*1000000.0
                    # Write in the Result file
                    #Result = open("results/node_"+str(nodeCounter),'a+')
                    #Result = open("results/node_"+str(node),'a+')
                    #print >> Result,"%.4f"%Time,Flow,Pressure

                    # Store results
                    timeResult[timestep] = Time
                    flowResult[nodeCounter,timestep] = Flow
                    pressureResult[nodeCounter,timestep] = Pressure

            timestep += 1


        pylab.plot(timeResult,flowResult[nodeCounter,:],'r',alpha=0.5,label=arteryLabel)

        for ref in refs:
            for refType in refTypes:
                filename = './reference/'+ref+'/'+refType+'/node_'+str(node)
                if os.path.exists(filename):
                    refData = np.genfromtxt(filename,delimiter='\t')
                    refLabel = ref + ' ' + refType + ' ' + arteryLabel
                    print('Found reference '+refLabel)
                    pylab.plot(refData[:,0],refData[:,1],'--r',alpha=0.5,label=refLabel)

        # Plot this node
        pylab.show()


        
        nodeCounter = nodeCounter+1        

    return

nodes = input('Nodes: ')
Post(nodes)

print "."
print "."
print "."
print "Processing Completed!"

#Popen(['gnuplot','PlotNode.p'])

