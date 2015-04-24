#!/usr/bin/env python

#> \file
#> \author Zhinuo Jenny Wang, 
#> \brief This is an example script to solve a finite elasticity equation using openCMISS calls in python.
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
#> Contributor(s): Sander Land
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

#> \example FiniteElasticity/UniAxialExtension/src/UniAxialExtensionExample.py
## Example script to solve a finite elasticity equation using openCMISS calls in python.
## \par Latest Builds:
## \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/UniAxialExtension/build-intel'>Linux Intel Build</a>
## \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/UniAxialExtension/build-gnu'>Linux GNU Build</a>
#<

# This python script solves one of 3 problems in the Cardiac Mechanics Benchmark Problems series organised by Sander Land in 2014-2015. 
 
# Problem 2: This python script solves the inflation of a truncated ellipsoide model 
# of the left ventricule. 
# The geometry and rule-based fibre field of the model is generated using a modified
# MATLAB script (originally written by Sander Land, modifications by Zhinuo Wang). 
# The constitutive model used is the Guccione transversely isotropic model, with 
# parameters as provided by the benchmark organiser. 
# W = C/2*exp(Q)
# where Q = b_f*E_11^2 + b_t*(E_22^2 + E_33^2 + E_23^2 + E_32^2) + b_fs*(E_12^2 + E_21^2 + E_13^2 + E_31^2)
# The material is fully incompressible.

import sys
import os
import numpy as np
import exfile
import math
from collections import OrderedDict

sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'], 'cm', 'bindings', 'python')))

# Initialise OpenCMISS
from lib import *
from numpy import array

### Set problem parameters3
numOfXi = 3
option = [1] # Trilinear
microstructure = 2 # Homogeneous fibre angles.
cellMLOption = [False]

### Set arbitrary user numbers which are unique to each object.
(coordinateSystemUserNumber,
 regionUserNumber,
 linearBasisUserNumber,
 cubicBasisUserNumber,
 generatedMeshUserNumber,
 meshUserNumber,
 decompositionUserNumber,
 geometricFieldUserNumber,
 equationsSetUserNumber,
 equationsSetFieldUserNumber,
 dependentFieldUserNumber,
 problemUserNumber,
 fibreFieldUserNumber,
 materialFieldUserNumber,
 deformedFieldUserNumber,
 strainFieldUserNumber,
 cellMLUserNumber,
 cellMLModelsFieldUserNumber,
 cellMLParametersFieldUserNumber,
 cellMLIntermediateFieldUserNumber) = range(1, 21)

# Set up geometric model
# longitudinal elements, circumferential elements, transmural elements
elems = [16,16,2]
cmd = '. ./matlab_batcher.sh benchmark_ellipse_linear ['+str(elems[0])+','+str(elems[1])+','+str(elems[2])+'];'
print cmd
if not os.path.exists("ellipse_benchmark_lin_"+str(elems[2])+"-"+str(elems[1])+"-"+str(elems[0])+".exnode"):
    os.system(cmd)

inputNodes = exfile.Exnode("ellipse_benchmark_lin_"+str(elems[2])+"-"+str(elems[1])+"-"+str(elems[0])+".exnode")
inputElems = exfile.Exelem("ellipse_benchmark_lin_"+str(elems[2])+"-"+str(elems[1])+"-"+str(elems[0])+".exelem")
num_apex_elem = elems[1]*elems[2]
apex_elems = []
for i in range(0, elems[2]):
    start = i * elems[1]*elems[0] + 1
    end = start + elems[1] - 1
    apex_elems.append([start, end])

# Set up region and CS
[numOfCompNodes, compNodeNum, CS, region] = BasicSetUp(regionUserNumber, coordinateSystemUserNumber)

# Set up tricubic Hermite basis functions
[linearBasis, colBasis] = BasisFunction(linearBasisUserNumber, numOfXi, option, collapsed=True)

# Set up mesh
mesh = CMISS.Mesh()
mesh.CreateStart(meshUserNumber, region, numOfXi)
mesh.NumberOfComponentsSet(1)
mesh.NumberOfElementsSet(inputElems.num_elements)

nodes = CMISS.Nodes()
nodes.CreateStart(region, inputNodes.num_nodes)
nodes.CreateFinish()

# Linear lagrange component
linearElem = CMISS.MeshElements()
linearElem.CreateStart(mesh, 1, linearBasis)
for elem in inputElems.elements:
    for i in range(1, elems[2]+1):
        if (elem.number >= (i-1)*elems[1]*elems[0]) & (elem.number <= i*elems[1]*elems[0]):
            if (elem.number >= apex_elems[i-1][0]) & (elem.number <= apex_elems[i-1][1]):
                linearElem.BasisSet(elem.number, colBasis)
                nodes = list(OrderedDict.fromkeys(elem.nodes))
                nodes = map(int, nodes)
                linearElem.NodesSet(elem.number, nodes)
            else:
                linearElem.NodesSet(elem.number, elem.nodes)

linearElem.CreateFinish()
mesh.CreateFinish()

# Set up decomposition for the mesh.
decomposition = DecompositionSetUp(decompositionUserNumber, mesh, numOfCompNodes)

# Set up geometric field.
geometricField = GeometricFieldSetUp(geometricFieldUserNumber, region, decomposition, option)

# Update the geometric field parameters manually
geometricField.ParameterSetUpdateStart(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)

def ExtractNodeCoords(nodes, field_name):
    coords = []
    for node_num in range(1, nodes.num_nodes+1):
        temp = []
        for component in [1,2,3]:
            component_name = ["x","y","z"][component-1]
            value = nodes.node_value(field_name, component_name, node_num, 1)
            temp.append(value)
        coords.append(temp)
    coords = tuple(coords)
    return coords

# Identify basal and endocardial nodes using mathematical formulation of the 
# truncated ellipsoid model.Intialise geometric field also.
basal_nodes = []
endocardial_nodes = []
all_nodes = ExtractNodeCoords(inputNodes, "coordinates")
for node_num in range(1, inputNodes.num_nodes+1):
    coord = []
    for component in [1,2,3]:
        value = all_nodes[node_num-1][component-1]
        geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1,
                                                node_num, component, value)
    coord = all_nodes[node_num-1]
    if coord[2] >= 5:
        basal_nodes.append(node_num)

    if coord[2] >= -17:
        u = math.acos(float(coord[2])/-17)
        u = -(math.pi - u)
        if coord[0] == 0:
            if coord[1] > 0.0:
                v = -math.pi/2
            else:
                v = math.pi/2
        elif coord[1] == 0:
            if coord[0] > 0.0:
                v = -math.pi
            else:
                v = 0
        else:
            v = math.atan(float(coord[1])/float(coord[0]))
        temp_x = abs(7*math.sin(u)*math.cos(v)) - abs(float(coord[0]))
        temp_y = abs(7*math.sin(u)*math.sin(v)) - abs(float(coord[1]))
        TOL = 1e-3
        if abs(temp_x) < TOL:
            if abs(temp_y) < TOL:
                endocardial_nodes.append(node_num)

unrefinedNodes = exfile.Exnode("unrefined_mesh.exnode")
unref_nodes = ExtractNodeCoords(unrefinedNodes, "coordinates")
no_base_unref_nodes = []
for node in unref_nodes:
    if node[2] < 5:
        no_base_unref_nodes.append(node)
no_base_unref_nodes = tuple(no_base_unref_nodes)

# Find nodes in current model which correspond to the unrefined nodes.
temp = [0,0,0]
eval_node_num = []
for i in range(0, inputNodes.num_nodes):
    if all_nodes[i] in no_base_unref_nodes:
        eval_node_num.append(i+1)
#print eval_node_num

geometricField.ParameterSetUpdateFinish(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)
    
# Export undeformed geometry.
GeometricFieldExport(region, "LVInflation_trilinear_undeformed_"+str(elems[2])+"-"+str(elems[1])+"-"+str(elems[0]))

# Set up fibre field.
fibreField = FibreFieldSetUp(fibreFieldUserNumber, region, decomposition, geometricField, option, microstructure,
                             inputNodes)

# Set up material field.
params = [10.0, 1.0, 1.0, 1.0]
materialField = MaterialFieldSetUp(materialFieldUserNumber, region, decomposition, geometricField, params, option,
                                   cellMLOption)

# Set up equations set
equationsSetField = CMISS.Field()  # Equations are also in a field
equationsSet = CMISS.EquationsSet()  # Initialise an equation set.
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField, CMISS.EquationsSetClasses.ELASTICITY,
                         CMISS.EquationsSetTypes.FINITE_ELASTICITY, CMISS.EquationsSetSubtypes.TRANSVERSE_ISOTROPIC_GUCCIONE,
                         equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()
print "----> Set up equations set <---\n"

# Set up material field in equations set.
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
equationsSet.MaterialsCreateFinish()

# Set up dependent field.
[dependentField, equationsSet] = DependentFieldSetUp(dependentFieldUserNumber, equationsSet, option, cellMLOption)

# Initialise dependent field.
DependentFieldInitialise(dependentField, geometricField, 0.0)

# Set up equations set
EquationsSetSetUp(equationsSet)

p = 0.0
for j in range(0,1):
    pressure_increments = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 1.0]
    tolerances = [1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5]
    iters = 1

    for i in range(0, len(pressure_increments)):
        increm = pressure_increments[i]
        p = p + increm
        tol = tolerances[i]
        print 'Applying pressure increment of: ', increm, ' using ', iters, ' iterations'
        print 'Current pressure is: ', p
        [problem, solverEquations] = ProblemSolverSetup(equationsSet, problemUserNumber, iters, tol, cellMLOption)
        BCEndoPressure(solverEquations, dependentField, endocardial_nodes, increm, basal_nodes, option)

        # Solve Problem
        problem.Solve()
        problem.Finalise()
        solverEquations.Finalise()

# Export results
materialField.Destroy()

# Export solutions
filename = "LVInflation_trilinear_"+str(elems[2])+"-"+str(elems[1])+"-"+str(elems[0])
ExportResults(dependentField, deformedFieldUserNumber, decomposition, region, filename, option)

# Evaluate displacement
deformed = exfile.Exnode("LVInflation_trilinear_"+str(elems[2])+"-"+str(elems[1])+"-"+str(elems[0])+".part0.exnode")
defNodes = ExtractNodeCoords(deformed, "DeformedGeometry")
defNodes = array(defNodes)
all_nodes = array(all_nodes)
disp = []
temp = [0,0,0]
for i in range(0, len(eval_node_num)):
    idx = eval_node_num[i]-1
    temp = defNodes[idx]- all_nodes[idx]
    disp.append(temp)
disp = array(disp)

rmse = math.sqrt(np.linalg.norm(disp))
# Write mse to a file
with open('rmse.txt', 'a') as f:
    f.write(str(elems[2])+" "+str(elems[1])+" "+str(elems[0])+" "+str(rmse)+"\n")
    f.close()

# To double check displacement,write out the resultant vectors in an exdata file.
with open("displacement"+str(elems[2])+"-"+str(elems[1])+"-"+str(elems[0])+".exdata", 'w') as fid:
    fid.write(' Group name: Displacement\n')
    fid.write(' #Fields=2\n')
    fid.write('  1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
    fid.write('   x.  Value index= 1, #Derivatives=0\n')
    fid.write('   y.  Value index= 2, #Derivatives=0\n')
    fid.write('   z.  Value index= 3, #Derivatives=0\n')
    fid.write('  1) displacement, coordinate, rectangular cartesian, #Components=3\n')
    fid.write('   x.  Value index= 1, #Derivatives=0\n')
    fid.write('   y.  Value index= 2, #Derivatives=0\n')
    fid.write('   z.  Value index= 3, #Derivatives=0\n')
    for i in range(0, len(disp)):
        idx = eval_node_num[i]
        fid.write(' Node:\t'+str(i+1)+'\n')
        fid.write(' '+str(all_nodes[idx-1][0])+' '+str(all_nodes[idx-1][1])+' '+str(all_nodes[idx-1][2])+'\n')
        fid.write(' '+str(disp[i][0])+' '+str(disp[i][1])+' '+str(disp[i][2])+'\n')


