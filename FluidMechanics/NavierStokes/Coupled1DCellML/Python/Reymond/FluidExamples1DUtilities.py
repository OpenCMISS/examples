#> \file
#> \author David Ladd
#> \brief This is a utility script, with routines for reading data associated with 1D Navier-Stokes problems, analysing the mesh, and outputting results.
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
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): Soroush Safaei
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. If you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. If you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>
#> OpenCMISS/examples/FluidMechanics/NavierStokes/Coupled1DCellML/Python/Reymond/1DFluidUtilities.py
#>

import csv
import math
import numpy as np
from numpy import linalg  

def GetNumberOfNodes(filename):
    ''' 1D Navier-Stokes problems are currently described using csv files- this routine gets the number of nodes for array allocation.
    '''
    # Read the node file
    try:
        with open(filename,'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            rownum = 0
            for row in reader:
                if (rownum == 0):
                    # Read the header row
                    header = row
                else:
                    # Read the number of nodes
                    if (rownum == 1):
                        numberOfNodes = int(row[5])
                    else:
                        break
                rownum+=1
            return(numberOfNodes);
    except IOError:
        print ('Could not open Node csv file: ' + filename)

def CsvNodeReader(filename,inputNodeNumbers,bifurcationNodeNumbers,trifurcationNodeNumbers,coupledNodeNumbers,nodeCoordinates,arteryLabels):
    ''' 1D Navier-Stokes problems are currently described using csv files- this routine reads in nodal data and returns it to the user.
    '''
    numberOfInputNodes     = 0
    numberOfBifurcations   = 0
    numberOfTrifurcations  = 0
    numberOfTerminalNodes  = 0
    # Read the node file
    try:
        with open(filename,'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            rownum = -1
            for row in reader:
                if (rownum == -1):
                    # Read the header row
                    header = row
                else:
                    # Read the number of nodes
                    if (rownum == 0):
                        numberOfNodesSpace = int(row[5])
                        totalNumberOfNodes = numberOfNodesSpace*3
                        xValues = np.empty([numberOfNodesSpace,4])
                        yValues = np.empty([numberOfNodesSpace,4])
                        zValues = np.empty([numberOfNodesSpace,4])
                        xValues.fill(np.nan)
                        yValues.fill(np.nan)
                        zValues.fill(np.nan)
                    # Initialise the coordinates
                    arteryLabels.append(row[0])
                    xValues[rownum,0] = float(row[1])
                    yValues[rownum,0] = float(row[2])
                    zValues[rownum,0] = float(row[3])
                    # Read the input nodes
                    if (row[4] == 'input'):
                        inputNodeNumbers.append(rownum+1)
                        numberOfInputNodes = numberOfInputNodes+1
                    # Read the bifurcation nodes
                    elif (row[4] == 'bifurcation'):
                        numberOfBifurcations+=1
                        bifurcationNodeNumbers.append(rownum+1)
                        xValues[rownum,1] = float(row[1])
                        yValues[rownum,1] = float(row[2])
                        zValues[rownum,1] = float(row[3])
                        xValues[rownum,2] = float(row[1])
                        yValues[rownum,2] = float(row[2])
                        zValues[rownum,2] = float(row[3])
                    # Read the trifurcation nodes
                    elif (row[4] == 'trifurcation'):
                        numberOfTrifurcations+=1
                        trifurcationNodeNumbers.append(rownum+1)
                        xValues[rownum,1] = float(row[1])
                        yValues[rownum,1] = float(row[2])
                        zValues[rownum,1] = float(row[3])
                        xValues[rownum,2] = float(row[1])
                        yValues[rownum,2] = float(row[2])
                        zValues[rownum,2] = float(row[3])
                        xValues[rownum,3] = float(row[1])
                        yValues[rownum,3] = float(row[2])
                        zValues[rownum,3] = float(row[3])
                    # Read the terminal nodes
                    elif (row[4] == 'terminal'):
                        coupledNodeNumbers.append(rownum+1)
                        numberOfTerminalNodes = numberOfTerminalNodes+1
                # Next line
                rownum+=1      
            # Create coordinates numpy array - init to NaN
            nodeCoordinates[:,:,0] = xValues[:,:]
            nodeCoordinates[:,:,1] = yValues[:,:]
            nodeCoordinates[:,:,2] = zValues[:,:]
    except IOError:
        print ('Could not open Node csv file: ' + filename)


def CsvElementReader(filename,elementNodes,bifurcationElements,trifurcationElements,numberOfBifurcations,numberOfTrifurcations):
    ''' 1D Navier-Stokes problems are currently described using csv files- this routine reads in element data and returns it to the user.
    '''
    try:
        # Read the element file
        with open(filename,'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            rownum = 0
            i = 0
            k = 0
            for row in reader:
                if (rownum == 0):
                    # Read the header row
                    header = row
                else:
                    # Read the number of elements
                    if (rownum == 1):
                        totalNumberOfElements = int(row[11])
                        #elementNodes          = (totalNumberOfElements+1)*[3*[0]]
                    # Read the element nodes
                    #elementNodes[rownum] = [int(row[1]),int(row[2]),int(row[3])]
                    elementNodes.append([int(row[1]),int(row[2]),int(row[3])])
                    # Read the bifurcation elements
                    if (row[4]):
                        i+=1
                        bifurcationElements[i] = [int(row[4]),int(row[5]),int(row[6])]
                    # Read the trifurcation elements
                    elif (row[7]):
                        k+=1
                        trifurcationElements[k] = [int(row[7]),int(row[8]),int(row[9]),int(row[10])]
                # Next line
                rownum+=1
    except IOError:
        print ('Could not open Element csv file: ' + filename)

def CsvMaterialReader(filename,A0,E,H):
    ''' 1D Navier-Stokes problems are currently described using csv files- this routine reads in material data and returns it to the user.
    '''
    try:
        # Read the element file
        with open(filename,'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            rownum = 0
            for row in reader:
                if (rownum == 0):
                    # Read the header row
                    header = row
                else:
                    A0[rownum][0] = float(row[1])
                    E [rownum][0] = float(row[2])
                    H [rownum][0] = float(row[3])
                # Next line
                rownum+=1
    except IOError:
        print ('Could not open Material csv file: ' + filename)


def GetMaxStableTimestep(elementNodes,QMax,nodeCoordinates,H,E,A0,Rho):
    ''' Indicates the min/max timestep for the provided mesh
    '''
    maxTimestep = 0.0
    numberOfElements = len(elementNodes)-1
    numberOfNodes = len(nodeCoordinates)
    # Check the element length
    elementNumber = [0]*(numberOfElements+1)
    elementLength = [0]*(numberOfElements+1)
    numberOfNodes = len(nodeCoordinates[:,0,0])
    eig  = [0]*(numberOfNodes+1)
    for i in range(1,numberOfElements+1):
        Node1 = elementNodes[i][0]
        Node2 = elementNodes[i][1]
        Node3 = elementNodes[i][2]
        Length1 = linalg.norm(nodeCoordinates[Node1-1,0,:]-nodeCoordinates[Node2-1,0,:])
        Length2 = linalg.norm(nodeCoordinates[Node2-1,0,:]-nodeCoordinates[Node3-1,0,:])
        elementNumber[i] = i
        elementLength[i] = Length1 + Length2
        elementLength[0] = elementLength[i]
        # print "Element %1.0f" %elementNumber[i], 
        # print "Length: %1.1f" %elementLength[i],
        # print "Length1: %1.1f" %Length1,
        # print "Length2: %1.1f" %Length2
    maxElementLength = max(elementLength)
    minElementLength = min(elementLength)
    print("Max Element Length: %1.3f" % maxElementLength)
    print("Min Element Length: %1.3f" % minElementLength)
               
    # Check the timestep
    dt   = [0]*(numberOfNodes+1)
    for i in range(1,numberOfNodes+1):
        beta   = (3.0*math.sqrt(math.pi)*H[i,0]*E[i,0])/(4.0*A0[i,0])
        eig[i] = QMax/A0[i,0] + (A0[i,0]**0.25)*(math.sqrt(beta/(2.0*Rho)))
        dt[i]  = ((3.0**(0.5))/3.0)*minElementLength/eig[i]
        dt[0]  = dt[i]
    maxTimestep = min(dt)
    print("Max allowable timestep: %3.5f" % maxTimestep )
    return(maxTimestep);
