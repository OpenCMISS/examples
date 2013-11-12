#!/usr/bin/env python

#> \file
#> \author David Ladd
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is openCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): 
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>


# Add Python bindings directory to PATH
import sys, os

import gzip
import numpy
import math
import re
from scipy import interpolate
from numpy import linalg,mean,sqrt
import pylab
import time
import matplotlib.pyplot as plt


class fieldInfo(object):
    'base class for info about fields'

    def __init__(self):

        self.group = ''
        self.numberOfFields = 0
        self.fieldNames = []
        self.numberOfFieldComponents = []

def tecplot_WriteRectilinearMesh(filename, X, Y, Z, vars):
    def pad(s, width):
        s2 = s
        while len(s2) < width:
            s2 = ' ' + s2
        if s2[0] != ' ':
            s2 = ' ' + s2
        if len(s2) > width:
            s2 = s2[:width]
        return s2
    def varline(vars, id, fw):
        s = ""
        for v in vars:
            s = s + pad(str(v[1][id]),fw)
        s = s + '\n'
        return s
 
    fw = 10 # field width
 
    f = open(filename, "wt")
 
    f.write('Variables="X","Y"')
    if len(Z) > 0:
        f.write(',"Z"')
    for v in vars:
        f.write(',"%s"' % v[0])
    f.write('\n\n')
 
    f.write('Zone I=' + pad(str(len(X)),6) + ',J=' + pad(str(len(Y)),6))
    if len(Z) > 0:
        f.write(',K=' + pad(str(len(Z)),6))
    f.write(', F=POINT\n')
 
    if len(Z) > 0:
        id = 0
        for k in xrange(len(Z)):
            for j in xrange(len(Y)):
                for i in xrange(len(X)):
                    f.write(pad(str(X[i]),fw) + pad(str(Y[j]),fw) + pad(str(Z[k]),fw))
                    f.write(varline(vars, id, fw))
                    id = id + 1
    else:
        id = 0
        for j in xrange(len(Y)):
            for i in xrange(len(X)):
                f.write(pad(str(X[i]),fw) + pad(str(Y[j]),fw))
                f.write(varline(vars, id, fw))
                id = id + 1
 
    f.close()

def vtkWriteAscii(filename, geometryData, velocityData, elementData):
    totalNumberOfNodes = geometryData.shape[0]
    totalNumberOfElements = elementData.shape[0]

    elementTransform = [0,2,8,6,18,20,26,24,1,5,7,3,19,23,25,21,9,11,17,15,12,14,10,16,4,22,13]

    f = open(filename, "wt")
    f.write('# vtk DataFile Version 2.0\n')
    f.write(filename + '\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')

    f.write('POINTS ' + str(totalNumberOfNodes) + ' float\n')
    for node in range(totalNumberOfNodes):
        f.write(str(geometryData[node,0]) + ' ' + str(geometryData[node,1]) + ' ' + str(geometryData[node,2]) + '\n')

    f.write('CELLS ' + str(totalNumberOfElements) + ' ' + str(27*totalNumberOfElements+totalNumberOfElements) +'\n')
    for element in range(totalNumberOfElements):
        line = '27 '
        for node in range(27):
            if node > 0:
                line += ' '
            line += str(int(elementData[element,elementTransform[node]]))
        line += '\n'
        f.write(line)

    f.write('CELL_TYPES  1\n29\n')

    f.write('POINT_DATA '+ str(totalNumberOfNodes) + '\n')
    f.write('VECTORS  Velocity float\n')
    for node in range(totalNumberOfNodes):
        f.write(str(velocityData[node,0]) + ' ' + str(velocityData[node,1]) + ' ' + str(velocityData[node,2]) + '\n')    

    f.close()


def findBetween( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


def readFirstHeader(f,info):

    #Read header info
#    print('reading HEADER')
    line=f.readline()
    s = re.findall(r'\d+', line)
    info.numberOfFields = int(s[0])
    for field in range(info.numberOfFields):
        line=f.readline()
        fieldName = findBetween(line, str(field + 1) + ') ', ',')
        info.fieldNames.append(fieldName)
        numberOfComponents = int(findBetween(line, '#Components=', '\n'))
        info.numberOfFieldComponents.append(numberOfComponents)
#        print('  number of components ' + str(numberOfComponents))
        for skip in range(numberOfComponents):
            line=f.readline()

def readExnodeHeader(f,numberOfFieldComponents):

    #Read header info
#    print('reading HEADER')
    line=f.readline()
    s = re.findall(r'\d+', line)
    numberOfFields = int(s[0])
#    print('  number of fields ' + str(numberOfFields))
    for field in range(numberOfFields):
        line=f.readline()
        fieldName = findBetween(line, str(field + 1) + ') ', ',')
        numberOfComponents = int(findBetween(line, '#Components=', '\n'))
        numberOfFieldComponents.append(numberOfComponents)
#        print('  number of components ' + str(numberOfComponents))
        for skip in range(numberOfComponents):
            line=f.readline()
#            print(line)


def readExnodeFile(filename,info,nodeData,totalNumberOfNodes):

    try:
        with open(filename):
            f = open(filename,"r")

            #Read header
            line=f.readline()
            info.group=findBetween(line, ' Group name: ', '\n')
            numberOfFieldComponents = []
            numberOfFields = 0
            readExnodeHeader(f,numberOfFieldComponents)
            numberOfFields = info.numberOfFields
            numberOfFieldComponents = info.numberOfFieldComponents

            #Read node data
            endOfFile = False
            while endOfFile == False:
                previousPosition = f.tell()
                line=f.readline()
                line = line.strip()
                if line:
                    if 'Node:' in line:
                        s = re.findall(r'\d+', line)
                        node = int(s[0])
                        for field in range(numberOfFields):
                            for component in range(numberOfFieldComponents[field]):
                                line=f.readline()
                                line = line.strip()
                                value = float(line)
                                if abs(value - 1.2345678806304932) < 1.0e-6:
                                    value =0.0
                                nodeData[node-1,field,component] = value

                    elif 'Fields' in line:
                        f.seek(previousPosition)
                        numberOfFieldComponents = []
                        numberOfFields = 0
                        readExnodeHeader(f,numberOfFieldComponents)
                        numberOfFields = len(numberOfFieldComponents)

                else:
                    endOfFile = True
                    f.close()
    except IOError:
       print ('Could not open file: ' + filename)



#=================================================================
# C o n t r o l   P a n e l
#=================================================================

ReynoldsNumbers = [5000]
meshResolution = [10,10]
compareSolutions = ['ghia.txt','erturk.txt','botella.txt']
path = "./output/"
vtuMesh = '/hpc/dlad004/velomap/extern/velomap_files/normal_aorta_1/aorta1IsoparamQuad.vtu'
meshName = 'Aorta'
resultFileRoot = meshName + '_t'
numberOfProcessors = 4
totalNumberOfNodes = 14494
fieldOfInterest = 2
componentsOfInterest = [i for i in range(3)]
elementFile = "./../../../input/Aorta/mesh/aorta1Isoparam.M"

#=================================================================
#=================================================================


nodeData = numpy.zeros([0,0,0,0])
numberOfRe = len(ReynoldsNumbers)

field = fieldInfo()        
filename = path + 'StaticSolution.part0.exnode'
try:
    with open(filename):
        firstFile = open(filename,"r")
        line=firstFile.readline()
        field.group=findBetween(line, ' Group name: ', '\n')
        readFirstHeader(firstFile,field)
        firstFile.close()
except IOError:
    print ('Could not open file: ' + filename)

nodeData = numpy.zeros([numberOfRe,totalNumberOfNodes,field.numberOfFields,max(field.numberOfFieldComponents)])

i = -1
for Re in ReynoldsNumbers:
    i+=1
    print('Reading data for Re ' + str(Re))
    for proc in range(numberOfProcessors):
        path = "./output/Re" + str(Re) + ;'Dim' + str(meshResolution[0]) + 'x' + str(meshResolution[1])'/'
        filename = path + 'LidDrivenCavity.part' + str(proc) +'.exnode'
        importNodeData = numpy.zeros([totalNumberOfNodes,field.numberOfFields,max(field.numberOfFieldComponents)])
        readExnodeFile(filename,field,importNodeData,totalNumberOfNodes)
        nodeData[i,:,:,:] += importNodeData[:,:,:]

nodeIds = []
nodeIds = [i for i in range(totalNumberOfNodes)]

for filename in compareSolutions:
    with open(filename,"r") as f:
        line = f.readline()
        dataInfo = line.strip().split(' ')
        compareAll = numpy.loadtxt(f)
        i = -1
        for Re in ReynoldsNumbers:
            i+=1
            j=0
            for refRe in dataInfo[1:]:
                j+=1
                if float(refRe) == Re:
                    compareData[i,:,:] = compareAll[j]
            
            

        for Re in dataInfo[1:]:
            if float(Re) in ReynoldsNumbers:
                compareData[Re,
                
                
        
        numpy.loadtxt(f)

    

plotData = True
if plotData:
    cycles = [i for i in range(numberOfFullCycles)]
    #pylab.plot(cycles,rmsCycle)
    pylab.bar(cycles,rmsCycle)
    pylab.xlabel('cycle')
    pylab.ylabel('RMS velocity component difference over cycle (cm/s)')
    pylab.title('Cycle Convergence ')
    pylab.grid(True)
    pylab.savefig('cycleConvergence')
    pylab.show()

    #times = [i*cmfeTimeIncrement for i in range(numberOfStepsInCycle)]
    times = [i*cmfeTimeIncrement for i in range(numberOfTimesteps)]
    pylab.plot(times,rmsPcv,'b-')
    pylab.plot(times,rmsPcv/meanVelocity,'g-')
    pylab.xlabel('Time')
    #pylab.ylabel('Percent RMS Difference CFD vs. PCV (%)')
    pylab.ylabel('RMS Difference CFD vs. PCV (cm/s)')
    pylab.title('RMS PCV vs. CFD')
    pylab.grid(True)
    pylab.savefig('cfdVsPcv')
    pylab.show()


    
