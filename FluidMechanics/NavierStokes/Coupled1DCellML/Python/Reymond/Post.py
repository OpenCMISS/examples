##################
#  Post Process  #
##################

import os,sys
import numpy as np
import subprocess
import pylab
from matplotlib import pyplot as plt
from subprocess import Popen, PIPE
import FluidExamples1DUtilities as Utilities

def Post(nodes):
    # Set the reference values
    Ts = 1.0#0.001                   # Time     (s)
    Qs = 100.0e-8                # Flow     (m3/s)

    # Set the time parameters
    cycleTime = 790.0
    numberOfCycles = 4.0
    solveTime      = cycleTime*numberOfCycles
    timeIncrement = 0.2
    outputFrequency = 10
    
    #Choose type of plot(s)
    plotFlow = True
    plotPressure = False

    # Set up numpy arrays
    numberOfTimesteps = int((solveTime/timeIncrement)/outputFrequency)
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
    refs = ['Reymond2009'] # ,'Reymond2011'
    refTypes = ['Model','Experimental']
    colours = ['r','c','b','g','k']

    # Node number for data extraction
    nodeLabels = []

    #Read in all node data
    print("Reading output file data (this can take a while...)")
    flowNodeData = np.zeros((totalNumberOfNodes,numberOfTimesteps))
    pressureNodeData = np.zeros((totalNumberOfNodes,numberOfTimesteps))
    #for timestep in range(0,int(solveTime/timeIncrement),outputFrequency):
    for timestep in range(numberOfTimesteps):
        filename = "output/MainTime_"+str(timestep*outputFrequency)+".part0.exnode"
        #print(filename)

        outputFile = open(filename)
        lines = outputFile.readlines()
        for i in range(0,len(lines)):
            if 'Node:' in lines[i]:
                nodeNumber = [int(s) for s in lines[i].split() if s.isdigit()][0]
                #print(nodeNumber)
                # Extract the variables from output files
                Flow = float(lines[i+4].split()[0])
                Flow = Flow*Qs*1000000.0
                Pressure = float(lines[i+12].split()[0])

                # Store results
                flowNodeData[nodeNumber-1,timestep] = Flow
                pressureNodeData[nodeNumber-1,timestep] = Pressure
        outputFile.close()
    timeResult =  np.arange(0.0,solveTime,Ts*timeIncrement*outputFrequency)

    if plotFlow:
        nodeCounter   = 0
        for node in nodes:
            # Get artery names
            arteryLabel = arteryLabels[node-1]
            print('Generating plot for ' + arteryLabel)
            nodeLabels.append(arteryLabel)

            # Do a subplot if looking at several nodes at once
            if numberOfNodes > 1:
                plt.subplot(int((numberOfNodes)/2)+numberOfNodes%2,2,nodeCounter+1)

            # Plot this node
            plt.plot(timeResult,flowNodeData[node-1,:],'k-',label='OpenCMISS Model')
            plt.title(arteryLabel+' (Node '+str(node)+')')
            colourNum = 0
            for ref in refs:
                for refType in refTypes:
                    filename = './reference/'+ref+'/'+refType+'/node_'+str(node)
                    if os.path.exists(filename):
                        refData = np.genfromtxt(filename,delimiter='\t')
                        numberOfRefTimesteps = len(refData[:,0])
                        refLabel = ref + ' ' + refType# + ' ' + arteryLabel
                        print('Found reference from '+refLabel)
                        refDataCycles=refData
                        addTime = np.zeros(numberOfRefTimesteps*numberOfCycles)
                        for cycle in range(1,int(numberOfCycles)):
                            refDataCycles = np.vstack((refDataCycles,refData))
                            addTime[numberOfRefTimesteps*cycle:]+=cycleTime
                        refDataCycles[:,0]=np.add(refDataCycles[:,0],addTime)
                        plt.plot(refDataCycles[:,0],refDataCycles[:,1],colours[colourNum],alpha=1.0,label=refLabel)
                        colourNum+=1

            nodeCounter = nodeCounter+1        
        # Plot all nodes
        plt.legend(loc = (1.2, 0.5))
        plt.show()

    if plotPressure:
        nodeCounter   = 0
        for node in nodes:
            # Get artery names
            arteryLabel = arteryLabels[node-1]
            print('Generating plot for ' + arteryLabel)
            nodeLabels.append(arteryLabel)

            # Do a subplot if looking at several nodes at once
            if numberOfNodes > 1:
                plt.subplot(int((numberOfNodes)/2)+numberOfNodes%2,2,nodeCounter+1)

            # Plot this node
            plt.plot(timeResult,pressureNodeData[node-1,:],'b-')
            plt.title(arteryLabel+' (Node '+str(node)+')')
            nodeCounter = nodeCounter+1        
        # Plot all nodes
        plt.show()

    return

if len(sys.argv) > 1:
    nodes = int(sys.argv[1:])
else:
    nodes = [1,93,120,127,170]
Post(nodes)

print "Processing Completed!"

#Popen(['gnuplot','PlotNode.p'])

