#!/usr/bin/env python

#> \file
#> \author Adam Reeve
#> \brief This is an example script showing how to calculate strain at a Xi location.
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

#> \example FiniteElasticity/MeasureStrain/src/MeasureStrainExample.py
## Example script showing how to calculate strain at a Xi location and calculate
## a Gauss point based strain field.
#<

from __future__ import print_function

import numpy as np
from opencmiss import iron


# Problem parameters:
extensionRatio = 1.1

numberGlobalElements = [1, 1, 1]
dimensions = [1.0, 1.0, 1.0]
numberOfXi = 3

constitutiveRelation = iron.EquationsSetSubtypes.MOONEY_RIVLIN
constitutiveParameters = [2.0, 6.0]

# User numbers for identifying OpenCMISS-Iron objects:
coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
equationsSetUserNumber = 1
(geometricFieldUserNumber,
    materialFieldUserNumber,
    dependentFieldUserNumber,
    equationsSetFieldUserNumber,
    strainFieldUserNumber) = range(1, 6)


# Create a 3D rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = iron.Region()
region.CreateStart(regionUserNumber, iron.WorldRegion)
region.LabelSet("Region")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

# Define basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.NumberOfXiSet(numberOfXi)
basis.InterpolationXiSet([
        iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE] * numberOfXi)
basis.QuadratureNumberOfGaussXiSet([2] * numberOfXi)
basis.CreateFinish()

# Start the creation of a generated mesh in the region
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber, region)
generatedMesh.TypeSet(iron.GeneratedMeshTypes.REGULAR)
generatedMesh.BasisSet([basis])
generatedMesh.ExtentSet(dimensions)
generatedMesh.NumberOfElementsSet(numberGlobalElements)
mesh = iron.Mesh()
generatedMesh.CreateFinish(meshUserNumber, mesh)

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = 1
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")
geometricField.CreateFinish()

# Update the geometric field parameters from generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create the equations_set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                             iron.EquationsSetTypes.FINITE_ELASTICITY,
                             constitutiveRelation]
equationsSet.CreateStart(equationsSetUserNumber, region, geometricField,
    equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create default materials field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
equationsSet.MaterialsCreateFinish()

# Create default dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
equationsSet.DependentCreateFinish()

# Create a Guass point based field for calculated strain
strainField = iron.Field()
strainField.CreateStart(strainFieldUserNumber, region)
strainField.MeshDecompositionSet(decomposition)
strainField.TypeSet(iron.FieldTypes.GENERAL)
strainField.GeometricFieldSet(geometricField)
strainField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
strainField.VariableTypesSet([iron.FieldVariableTypes.U])
strainField.VariableLabelSet(iron.FieldVariableTypes.U, "Strain")
strainField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 6)
for component in range(1, 7):
    strainField.ComponentInterpolationSet(
            iron.FieldVariableTypes.U, component,
            iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
strainField.CreateFinish()

equationsSet.DerivedCreateStart(strainFieldUserNumber, strainField)
equationsSet.DerivedVariableSet(iron.EquationsSetDerivedTypes.STRAIN, iron.FieldVariableTypes.U)
equationsSet.DerivedCreateFinish()

# Set constitutive parameters
for component, parameter in enumerate(constitutiveParameters, 1):
    iron.Field.ComponentValuesInitialiseDP(
        materialField,iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES,
        component, parameter)

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
equationsSet.EquationsCreateFinish()

# Set deformed geometry
nodes = iron.Nodes()
region.NodesGet(nodes)
for node in range(1, nodes.NumberOfNodesGet() + 1):
    position = [geometricField.ParameterSetGetNode(
                iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                1, 1, node, component)
            for component in range(1, 4)]
    ratios = [extensionRatio,
            (1.0 / np.sqrt(extensionRatio)),
            (1.0 / np.sqrt(extensionRatio))]
    deformedPosition = [p * r for p, r in zip(position, ratios)]
    for component, value in enumerate(deformedPosition, 1):
        dependentField.ParameterSetUpdateNode(
                iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                1, 1, node, component, value)

# Calculate expected strain
FAxial = extensionRatio
FTransverse = (1.0 / np.sqrt(extensionRatio))
F = np.array([
    [FAxial,         0.0,         0.0],
    [   0.0, FTransverse,         0.0],
    [   0.0,         0.0, FTransverse]])
C = (F.T).dot(F)
E = 0.5 * (C - np.eye(3))

def matrixFromSymmetricComponents(components):
    return np.array([
        [components[0], components[1], components[2]],
        [components[1], components[3], components[4]],
        [components[2], components[4], components[5]],
        ])

print("Expected strain:")
print(E)

print("Calculated strain:")
elementNumber = 1
xiPosition = [0.5, 0.5, 0.5]
calculatedStrain = equationsSet.StrainInterpolateXi(elementNumber, xiPosition)
calculatedStrainTensor = matrixFromSymmetricComponents(calculatedStrain)
print(calculatedStrainTensor)

for i in range(3):
    for j in range(3):
        assert(abs(calculatedStrainTensor[i,j] - E[i,j]) < 1.0e-10)

# Calculate Gauss point based strain field:
equationsSet.DerivedVariableCalculate(iron.EquationsSetDerivedTypes.STRAIN)
# Check strain at a Gauss point:
elementNumber = 1
gaussPointNumber = 1
componentNumber = 4
gaussPointStrain = strainField.ParameterSetGetGaussPoint(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
        gaussPointNumber, elementNumber, componentNumber)
assert(abs(gaussPointStrain - E[1, 1]) < 1.0e-10)
