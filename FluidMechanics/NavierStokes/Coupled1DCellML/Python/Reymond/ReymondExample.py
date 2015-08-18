#> \file
#> \author David Ladd
#> \brief This is an example program to solve for flow using 1D transient Navier-Stokes 
#>  over an arterial tree with coupled 0D lumped models (RCR) defined in CellML. The geometry and
#>  boundary conditions are based on published data from Reymond et al. 2011: 'Validation of a patient-specific one-dimensional model of the systemic arterial tree'
#>  Results are compared against the data presented in this paper.
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
#> OpenCMISS/examples/FluidMechanics/NavierStokes/Coupled1DCellML/Python/Reymond/ReymondExample.py
#>

#================================================================================================================================
#  Initialise OpenCMISS and any other needed libraries
#================================================================================================================================
import numpy as np
import math,csv,time,sys,os,glob,shutil
import FluidExamples1DUtilities as Utilities1D

sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))
from opencmiss import CMISS

#================================================================================================================================
#  Set up field and system values 
#================================================================================================================================

(CoordinateSystemUserNumber,
 BasisUserNumber,
 RegionUserNumber,
 MeshUserNumber,
 DecompositionUserNumber,
 GeometricFieldUserNumber,
 DependentFieldUserNumber,
 MaterialsFieldUserNumber,  
 IndependentFieldUserNumber,
 EquationsSetUserNumberCharacteristic,
 EquationsSetUserNumberNavierStokes,
 EquationsSetFieldUserNumberCharacteristic,
 EquationsSetFieldUserNumberNavierStokes,
 ProblemUserNumber,
 CellMLUserNumber,
 CellMLModelsFieldUserNumber,
 CellMLStateFieldUserNumber,
 CellMLIntermediateFieldUserNumber,
 CellMLParametersFieldUserNumber,
 AnalyticFieldUserNumber) = range(1,21)

# Solver user numbers
SolverDAEUserNumber            = 1
SolverCharacteristicUserNumber = 2
SolverNavierStokesUserNumber   = 3

# Other system constants
numberOfDimensions     = 1  #(One-dimensional)
numberOfComponents     = 2  #(Flow & Area)

# Get the computational nodes info
numberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
computationalNodeNumber    = CMISS.ComputationalNodeNumberGet()

#================================================================================================================================
#  Problem Control Panel
#================================================================================================================================

# Set the flags
RCRBoundaries            = True   # Set to use coupled 0D Windkessel models (from CellML) at model outlet boundaries
nonReflecting            = False    # Set to use non-reflecting outlet boundaries
CheckTimestepStability   = False   # Set to do a basic check of the stability of the hyperbolic problem based on the timestep size
initialiseFromFile       = False   # Set to initialise values
ProgressDiagnostics      = True    # Set to diagnostics

if(nonReflecting and RCRBoundaries):
    sys.exit('Please set either RCR or non-reflecting boundaries- not both.')

#================================================================================================================================
#  Mesh Reading
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> Reading geometry from files... << == "

# Read nodes
inputNodeNumbers        = []
bifurcationNodeNumbers  = []
trifurcationNodeNumbers = []
coupledNodeNumbers      = []
arteryLabels            = []
filename = 'input/Node.csv'
numberOfNodes = Utilities1D.GetNumberOfNodes(filename)
nodeCoordinates = np.zeros([numberOfNodes,4,3])
Utilities1D.CsvNodeReader(filename,inputNodeNumbers,bifurcationNodeNumbers,trifurcationNodeNumbers,coupledNodeNumbers,
                        nodeCoordinates,arteryLabels)
numberOfInputNodes     = len(inputNodeNumbers)
numberOfBifurcations   = len(bifurcationNodeNumbers)
numberOfTrifurcations  = len(trifurcationNodeNumbers)
numberOfTerminalNodes  = len(coupledNodeNumbers)

# Read elements
elementNodes = []
elementNodes.append([0,0,0])
bifurcationElements = (numberOfBifurcations+1)*[3*[0]]
trifurcationElements = (numberOfTrifurcations+1)*[4*[0]]
Utilities1D.CsvElementReader('input/Element.csv',elementNodes,bifurcationElements,trifurcationElements,numberOfBifurcations,numberOfTrifurcations)
numberOfElements = len(elementNodes)-1
        
if (ProgressDiagnostics):
    print " Number of nodes: " + str(numberOfNodes)
    print " Number of elements: " + str(numberOfElements)
    print " Input at nodes: " + str(inputNodeNumbers)
    print " Bifurcations at nodes: " + str(bifurcationNodeNumbers)
    print " Trifurcations at nodes: " + str(trifurcationNodeNumbers)
    print " Terminal at nodes: " + str(coupledNodeNumbers)
    print " == >> Finished reading geometry... << == "

#================================================================================================================================
#  Initial Data & Default Values
#================================================================================================================================

# Set the material parameters
Rho  = 1050.0     # Density     (kg/m3)
Mu   = 0.004      # Viscosity   (Pa.s)
G0   = 0.0        # Gravitational acceleration (m/s2)
Pext = 0.0        # External pressure (Pa)
Alpha = 1.0       # Flow profile type

# Material parameter scaling factors
Ls = 1000.0              # Length   (m -> mm)
Ts = 1000.0              # Time     (s -> ms)
Ms = 1000.0              # Mass     (kg -> g)
Qs    = (Ls**3.0)/Ts     # Flow             (m3/s)  
As    = Ls**2.0          # Area             (m2)
Hs    = Ls               # vessel thickness (m)
Es    = Ms/(Ls*Ts**2.0)  # Elasticity Pa    (kg/(ms2) --> g/(mm.ms^2)
Rhos  = Ms/(Ls**3.0)     # Density          (kg/m3)
Mus   = Ms/(Ls*Ts)       # Viscosity        (kg/(ms))
Ps    = Ms/(Ls*Ts**2.0)  # Pressure         (kg/(ms2))
Gs    = Ls/(Ts**2.0)     # Acceleration    (m/s2)

# Initialise the node-based parameters
A0   = np.zeros((numberOfNodes+1,4))  # Area        (m2)
H    = np.zeros((numberOfNodes+1,4))  # Thickness   (m)
E    = np.zeros((numberOfNodes+1,4))  # Elasticity  (Pa)
# Read the MATERIAL csv file
Utilities1D.CsvMaterialReader('input/Material.csv',A0,E,H)

# Apply scale factors        
Rho = Rho*Rhos
Mu  = Mu*Mus
P   = Pext*Ps
A0  = A0*As
E   = E*Es
H   = H*Hs
G0  = G0*Gs

Q  = np.zeros((numberOfNodes+1,4))
A  = np.zeros((numberOfNodes+1,4))
dQ = np.zeros((numberOfNodes+1,4))
dA = np.zeros((numberOfNodes+1,4))

for bifIdx in range(1,numberOfBifurcations+1):
    nodeIdx = bifurcationNodeNumbers[bifIdx-1]
    for versionIdx in range(1,3):
        A0[nodeIdx][versionIdx] = A0[elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]
        E [nodeIdx][versionIdx] = E [elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]
        H [nodeIdx][versionIdx] = H [elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]
for trifIdx in range(1,numberOfTrifurcations+1):
    nodeIdx = trifurcationNodeNumbers[trifIdx-1]
    for versionIdx in range(1,4):
        A0[nodeIdx][versionIdx] = A0[elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]
        E [nodeIdx][versionIdx] = E [elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]
        H [nodeIdx][versionIdx] = H [elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]

# Start with Q=0, A=A0 state
A = A0

# Or initialise from init file
if (initialiseFromFile):
    init = np.zeros([numberOfNodes+1,4,4])
    init = np.load('./input/init.npy')
    Q[1:numberOfNodes+1,:] = init[:,0,:]
    A[1:numberOfNodes+1,:] = init[:,1,:]
    dQ[1:numberOfNodes+1,:] = init[:,2,:]
    dA[1:numberOfNodes+1,:] = init[:,3,:]

# Set the output parameters
# (NONE/PROGRESS/TIMING/SOLVER/MATRIX)
dynamicSolverNavierStokesOutputType    = CMISS.SolverOutputTypes.NONE
nonlinearSolverNavierStokesOutputType  = CMISS.SolverOutputTypes.NONE
nonlinearSolverCharacteristicsOutputType = CMISS.SolverOutputTypes.NONE
linearSolverCharacteristicOutputType    = CMISS.SolverOutputTypes.NONE
linearSolverNavierStokesOutputType     = CMISS.SolverOutputTypes.NONE
# (NONE/TIMING/SOLVER/MATRIX)
cmissSolverOutputType = CMISS.SolverOutputTypes.NONE
dynamicSolverNavierStokesOutputFrequency = 10

# Set the time parameters
numberOfPeriods = 4.0
timePeriod      = 790.
timeIncrement   = 0.2
startTime       = 0.0
stopTime  = numberOfPeriods*timePeriod
dynamicSolverNavierStokesTheta = [1.0]

# Set the solver parameters
relativeToleranceNonlinearNavierStokes   = 1.0E-05  # default: 1.0E-05
absoluteToleranceNonlinearNavierStokes   = 1.0E-08  # default: 1.0E-10
solutionToleranceNonlinearNavierStokes   = 1.0E-05  # default: 1.0E-05
relativeToleranceLinearNavierStokes      = 1.0E-05  # default: 1.0E-05
absoluteToleranceLinearNavierStokes      = 1.0E-08  # default: 1.0E-10
relativeToleranceNonlinearCharacteristic = 1.0E-05  # default: 1.0E-05
absoluteToleranceNonlinearCharacteristic = 1.0E-08  # default: 1.0E-10
solutionToleranceNonlinearCharacteristic = 1.0E-05  # default: 1.0E-05
relativeToleranceLinearCharacteristic    = 1.0E-05  # default: 1.0E-05
absoluteToleranceLinearCharacteristic    = 1.0E-08  # default: 1.0E-10

DIVERGENCE_TOLERANCE = 1.0E+10  # default: 1.0E+05
MAXIMUM_ITERATIONS   = 100000   # default: 100000
RESTART_VALUE        = 3000     # default: 30

# N-S/C coupling tolerance
couplingTolerance1D = 1.0E+10
# 1D-0D coupling tolerance
couplingTolerance1D0D = 0.001

# Navier-Stokes solver
if(RCRBoundaries):
    EquationsSetSubtype = CMISS.EquationsSetSubtypes.COUPLED1D0D_NAVIER_STOKES
    # Characteristic solver
    EquationsSetCharacteristicSubtype = CMISS.EquationsSetSubtypes.CHARACTERISTIC
    ProblemSubtype = CMISS.ProblemSubTypes.COUPLED1D0D_NAVIER_STOKES
else:
    EquationsSetSubtype = CMISS.EquationsSetSubtypes.TRANSIENT1D_NAVIER_STOKES
    # Characteristic solver
    EquationsSetCharacteristicSubtype = CMISS.EquationsSetSubtypes.CHARACTERISTIC
    ProblemSubtype = CMISS.ProblemSubTypes.TRANSIENT1D_NAVIER_STOKES

#================================================================================================================================
#  Coordinate System
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> COORDINATE SYSTEM << == "

# Start the creation of RC coordinate system
CoordinateSystem = CMISS.CoordinateSystem()
CoordinateSystem.CreateStart(CoordinateSystemUserNumber)
CoordinateSystem.DimensionSet(3)
CoordinateSystem.CreateFinish()

#================================================================================================================================
#  Region
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> REGION << == "

# Start the creation of  region
Region = CMISS.Region()
Region.CreateStart(RegionUserNumber,CMISS.WorldRegion)
Region.label = "ArterialSystem"
Region.coordinateSystem = CoordinateSystem
Region.CreateFinish()

#================================================================================================================================
#  Bases
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> BASIS << == "

# Start the creation of  bases
basisXiGauss = 3
Basis = CMISS.Basis()
Basis.CreateStart(BasisUserNumber)
Basis.type = CMISS.BasisTypes.LAGRANGE_HERMITE_TP
Basis.numberOfXi = numberOfDimensions
Basis.interpolationXi = [CMISS.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]
Basis.quadratureNumberOfGaussXi = [basisXiGauss]
Basis.CreateFinish()

#================================================================================================================================
#  Nodes
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> NODES << == "

# Start the creation of mesh nodes
Nodes = CMISS.Nodes()
Nodes.CreateStart(Region,numberOfNodes)
Nodes.CreateFinish()

#================================================================================================================================
#  Mesh
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> MESH << == "

# Start the creation of  mesh
Mesh = CMISS.Mesh()
Mesh.CreateStart(MeshUserNumber,Region,numberOfDimensions)
Mesh.NumberOfElementsSet(numberOfElements)
meshNumberOfComponents = 1
# Specify the mesh components
Mesh.NumberOfComponentsSet(meshNumberOfComponents)
# Specify the mesh components
MeshElements = CMISS.MeshElements()
meshComponentNumber = 1

# Specify the  mesh component
MeshElements.CreateStart(Mesh,meshComponentNumber,Basis)
for elemIdx in range(1,numberOfElements+1):
    MeshElements.NodesSet(elemIdx,elementNodes[elemIdx])
for bifIdx in range(1,numberOfBifurcations+1):
    MeshElements.LocalElementNodeVersionSet(int(bifurcationElements[bifIdx][0]),1,1,3)
    MeshElements.LocalElementNodeVersionSet(int(bifurcationElements[bifIdx][1]),2,1,1) 
    MeshElements.LocalElementNodeVersionSet(int(bifurcationElements[bifIdx][2]),3,1,1) 
for trifIdx in range(1,numberOfTrifurcations+1):
    MeshElements.LocalElementNodeVersionSet(int(trifurcationElements[trifIdx][0]),1,1,3)
    MeshElements.LocalElementNodeVersionSet(int(trifurcationElements[trifIdx][1]),2,1,1) 
    MeshElements.LocalElementNodeVersionSet(int(trifurcationElements[trifIdx][2]),3,1,1) 
    MeshElements.LocalElementNodeVersionSet(int(trifurcationElements[trifIdx][3]),4,1,1) 
MeshElements.CreateFinish()

# Finish the creation of the mesh
Mesh.CreateFinish()

#================================================================================================================================
#  Decomposition
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> MESH DECOMPOSITION << == "

# Start the creation of  mesh decomposition
Decomposition = CMISS.Decomposition()
Decomposition.CreateStart(DecompositionUserNumber,Mesh)
Decomposition.TypeSet(CMISS.DecompositionTypes.CALCULATED)
Decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
Decomposition.CreateFinish()

#================================================================================================================================
#  Geometric Field
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> GEOMETRIC FIELD << == "

# Start the creation of  geometric field
GeometricField = CMISS.Field()
GeometricField.CreateStart(GeometricFieldUserNumber,Region)
GeometricField.NumberOfVariablesSet(1)
GeometricField.VariableLabelSet(CMISS.FieldVariableTypes.U,'Coordinates')
GeometricField.TypeSet = CMISS.FieldTypes.GEOMETRIC
GeometricField.meshDecomposition = Decomposition
GeometricField.ScalingTypeSet = CMISS.FieldScalingTypes.NONE
# Set the mesh component to be used by the geometric field components
for componentNumber in range(1,CoordinateSystem.dimension+1):
    GeometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,componentNumber,
                                             meshComponentNumber)
GeometricField.CreateFinish()

# Set the geometric field values
for node in range(numberOfNodes):
    nodeNumber = node+1
    nodeDomain = Decomposition.NodeDomainGet(nodeNumber,meshComponentNumber)
    if (nodeDomain == computationalNodeNumber):
        for version in range(4):
            versionNumber = version + 1
            # If the version is undefined for this node (not a bi/trifurcation), continue to next node
            if (np.isnan(nodeCoordinates[node,version,0])):
                break
            else:
                for component in range(3):
                    componentNumber = component+1
                    GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
                                                            versionNumber,1,nodeNumber,componentNumber,nodeCoordinates[node,version,component])

# Finish the parameter update
GeometricField.ParameterSetUpdateStart(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)
GeometricField.ParameterSetUpdateFinish(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)     

# Export Geometry
Fields = CMISS.Fields()
Fields.CreateRegion(Region)
Fields.NodesExport("Geometry","FORTRAN")
Fields.ElementsExport("Geometry","FORTRAN")
Fields.Finalise()

#================================================================================================================================
#  Equations Sets
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> EQUATIONS SET << == "

# Create the equations set for CHARACTERISTIC
EquationsSetCharacteristic = CMISS.EquationsSet()
EquationsSetFieldCharacteristic = CMISS.Field()
# Set the equations set to be a static nonlinear problem
EquationsSetCharacteristic.CreateStart(EquationsSetUserNumberCharacteristic,Region,GeometricField,
                                       CMISS.EquationsSetClasses.FLUID_MECHANICS,CMISS.EquationsSetTypes.CHARACTERISTIC_EQUATION,
                                       EquationsSetCharacteristicSubtype,EquationsSetFieldUserNumberCharacteristic,EquationsSetFieldCharacteristic)
EquationsSetCharacteristic.CreateFinish()

# Create the equations set for NAVIER-STOKES
EquationsSetNavierStokes = CMISS.EquationsSet()
EquationsSetFieldNavierStokes = CMISS.Field()
# Set the equations set to be a dynamic nonlinear problem
EquationsSetNavierStokes.CreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField,
    CMISS.EquationsSetClasses.FLUID_MECHANICS,CMISS.EquationsSetTypes.NAVIER_STOKES_EQUATION,
     EquationsSetSubtype,EquationsSetFieldUserNumberNavierStokes,EquationsSetFieldNavierStokes)
EquationsSetNavierStokes.CreateFinish()

#================================================================================================================================
#  Dependent Field
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> DEPENDENT FIELD << == "

# CHARACTERISTIC
# Create the equations set dependent field variables
DependentFieldNavierStokes = CMISS.Field()
EquationsSetCharacteristic.DependentCreateStart(DependentFieldUserNumber,DependentFieldNavierStokes)
DependentFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.U,'General')
DependentFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.DELUDELN,'Derivatives')
DependentFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.V,'Characteristics')
if (RCRBoundaries):
    DependentFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.U1,'CellML Q and P')
DependentFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.U2,'Pressure')
# Set the mesh component to be used by the field components.
# Flow & Area
DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,1,meshComponentNumber)
DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,2,meshComponentNumber)
# Derivatives
DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN,1,meshComponentNumber)
DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN,2,meshComponentNumber)
# Riemann
DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.V,1,meshComponentNumber)
DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.V,2,meshComponentNumber)
# qCellML & pCellml
if (RCRBoundaries):
    DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U1,1,meshComponentNumber)
    DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U1,2,meshComponentNumber)
# Pressure
DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U2,1,meshComponentNumber)
DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U2,2,meshComponentNumber)

EquationsSetCharacteristic.DependentCreateFinish()

#------------------

# NAVIER-STOKES
EquationsSetNavierStokes.DependentCreateStart(DependentFieldUserNumber,DependentFieldNavierStokes)
EquationsSetNavierStokes.DependentCreateFinish()

DependentFieldNavierStokes.ParameterSetCreate(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.PREVIOUS_VALUES)

# Initialise the dependent field variables
for nodeIdx in range (1,numberOfNodes+1):
    nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumber)
    if (nodeDomain == computationalNodeNumber):
        if (nodeIdx in trifurcationNodeNumbers):
            versions = [1,2,3,4]
        elif (nodeIdx in bifurcationNodeNumbers):
            versions = [1,2,3]
        else:
            versions = [1]
        for versionIdx in versions:
            # U variables
            DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
                                                                versionIdx,1,nodeIdx,1,Q[nodeIdx][versionIdx-1])
            DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
                                                                versionIdx,1,nodeIdx,2,A[nodeIdx][versionIdx-1])
            DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.PREVIOUS_VALUES,
                                                                versionIdx,1,nodeIdx,1,Q[nodeIdx][versionIdx-1])
            DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.PREVIOUS_VALUES,
                                                                versionIdx,1,nodeIdx,2,A[nodeIdx][versionIdx-1])
            # delUdelN variables
            DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.DELUDELN,CMISS.FieldParameterSetTypes.VALUES,
                                                                versionIdx,1,nodeIdx,1,dQ[nodeIdx][versionIdx-1])
            DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.DELUDELN,CMISS.FieldParameterSetTypes.VALUES,
                                                                versionIdx,1,nodeIdx,2,dA[nodeIdx][versionIdx-1])

# revert default version to 1
versionIdx = 1
             
# Finish the parameter update
DependentFieldNavierStokes.ParameterSetUpdateStart(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)
DependentFieldNavierStokes.ParameterSetUpdateFinish(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)   

#================================================================================================================================
#  Materials Field
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> MATERIALS FIELD << == "

# CHARACTERISTIC
# Create the equations set materials field variables 
MaterialsFieldNavierStokes = CMISS.Field()
EquationsSetCharacteristic.MaterialsCreateStart(MaterialsFieldUserNumber,MaterialsFieldNavierStokes)
MaterialsFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.U,'MaterialsConstants')
MaterialsFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.V,'MaterialsVariables')
# Set the mesh component to be used by the field components.
for componentNumber in range(1,4):
    MaterialsFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.V,componentNumber,meshComponentNumber)
EquationsSetCharacteristic.MaterialsCreateFinish()

#------------------

# NAVIER-STOKES
EquationsSetNavierStokes.MaterialsCreateStart(MaterialsFieldUserNumber,MaterialsFieldNavierStokes)
EquationsSetNavierStokes.MaterialsCreateFinish()

# Set the materials field constants
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,
                                                       CMISS.FieldParameterSetTypes.VALUES,1,Mu)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,
                                                       CMISS.FieldParameterSetTypes.VALUES,2,Rho)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,
                                                       CMISS.FieldParameterSetTypes.VALUES,3,Alpha)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,
                                                       CMISS.FieldParameterSetTypes.VALUES,4,Pext)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,
                                                       CMISS.FieldParameterSetTypes.VALUES,5,Ls)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,
                                                       CMISS.FieldParameterSetTypes.VALUES,6,Ts)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,
                                                       CMISS.FieldParameterSetTypes.VALUES,7,Ms)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,
                                                       CMISS.FieldParameterSetTypes.VALUES,8,G0)

# Initialise the materials field variables (A0,E,H)
bifIdx = 0
trifIdx = 0
for nodeIdx in range(1,numberOfNodes+1,1):
    nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumber)
    if (nodeDomain == computationalNodeNumber):
        if (nodeIdx in trifurcationNodeNumbers):
            versions = [1,2,3,4]
        elif (nodeIdx in bifurcationNodeNumbers):
            versions = [1,2,3]
        else:
            versions = [1]
        for versionIdx in versions:
            MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,CMISS.FieldParameterSetTypes.VALUES,
                                                                versionIdx,1,nodeIdx,1,A0[nodeIdx][versionIdx-1])
            MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,CMISS.FieldParameterSetTypes.VALUES,
                                                                versionIdx,1,nodeIdx,2,E[nodeIdx][versionIdx-1])
            MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,CMISS.FieldParameterSetTypes.VALUES,
                                                                versionIdx,1,nodeIdx,3,H[nodeIdx][versionIdx-1])

# Finish the parameter update
MaterialsFieldNavierStokes.ParameterSetUpdateStart(CMISS.FieldVariableTypes.V,CMISS.FieldParameterSetTypes.VALUES)
MaterialsFieldNavierStokes.ParameterSetUpdateFinish(CMISS.FieldVariableTypes.V,CMISS.FieldParameterSetTypes.VALUES)

#================================================================================================================================
# Independent Field
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> INDEPENDENT FIELD << == "

# CHARACTERISTIC
# Create the equations set independent field variables  
IndependentFieldNavierStokes = CMISS.Field()
EquationsSetCharacteristic.IndependentCreateStart(IndependentFieldUserNumber,IndependentFieldNavierStokes)
IndependentFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.U,'Normal Wave Direction')
# Set the mesh component to be used by the field components.
IndependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,1,meshComponentNumber)
IndependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,2,meshComponentNumber)
EquationsSetCharacteristic.IndependentCreateFinish()

#------------------

# NAVIER-STOKES
EquationsSetNavierStokes.IndependentCreateStart(IndependentFieldUserNumber,IndependentFieldNavierStokes)
EquationsSetNavierStokes.IndependentCreateFinish()

# Set the normal wave direction for bifurcation
for bifIdx in range (1,numberOfBifurcations+1):
    nodeIdx = bifurcationNodeNumbers[bifIdx-1]
    nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumber)
    if (nodeDomain == computationalNodeNumber):
        # Incoming(parent)
        IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
                                                              1,1,nodeIdx,1,1.0)
        # Outgoing(branches)
        IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
                                                              2,1,nodeIdx,2,-1.0)
        IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
                                                              3,1,nodeIdx,2,-1.0)
# Set the normal wave direction for trifurcation
for trifIdx in range (1,numberOfTrifurcations+1):
    nodeIdx = trifurcationNodeNumbers[trifIdx-1]
    nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumber)
    if (nodeDomain == computationalNodeNumber):
        # Incoming(parent)
        IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
                                                              1,1,nodeIdx,1,1.0)
        # Outgoing(branches)
        IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
                                                              2,1,nodeIdx,2,-1.0)
        IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
                                                              3,1,nodeIdx,2,-1.0)
        IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
                                                              4,1,nodeIdx,2,-1.0)

# Set the normal wave direction for terminal
if (RCRBoundaries or nonReflecting):
    for terminalIdx in range (1,numberOfTerminalNodes+1):
        nodeIdx = coupledNodeNumbers[terminalIdx-1]
        nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumber)
        if (nodeDomain == computationalNodeNumber):
            # Incoming (parent) - outgoing component to come from 0D
            versionIdx = 1
            IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
                                                                  versionIdx,1,nodeIdx,1,1.0)

# Finish the parameter update
IndependentFieldNavierStokes.ParameterSetUpdateStart(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)
IndependentFieldNavierStokes.ParameterSetUpdateFinish(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)

#================================================================================================================================
# Analytic Field
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> ANALYTIC FIELD << == "

AnalyticFieldNavierStokes = CMISS.Field()
EquationsSetNavierStokes.AnalyticCreateStart(CMISS.NavierStokesAnalyticFunctionTypes.SPLINT_FROM_FILE,AnalyticFieldUserNumber,
                                             AnalyticFieldNavierStokes)
AnalyticFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.U,'Input Flow')
EquationsSetNavierStokes.AnalyticCreateFinish()

#DOC-START cellml define field maps
#================================================================================================================================
#  RCR CellML Model Maps
#================================================================================================================================


if (RCRBoundaries):

    #----------------------------------------------------------------------------------------------------------------------------
    # Description
    #----------------------------------------------------------------------------------------------------------------------------
    # A CellML OD model is used to provide the impedance from the downstream vascular bed beyond the termination
    # point of the 1D model. This is iteratively coupled with the the 1D solver. In the case of a simple resistance
    # model, P=RQ, which is analogous to Ohm's law: V=IR. A variable map copies the guess for the FlowRate, Q at 
    # the boundary from the OpenCMISS Dependent Field to the CellML equation, which then returns presssure, P.
    # The initial guess value for Q is taken from the previous time step or is 0 for t=0. In OpenCMISS this P value is 
    # then used to compute a new Area value based on the P-A relationship and the Riemann variable W_2, which gives a
    # new value for Q until the values for Q and P converge within tolerance of the previous value.
    #----------------------------------------------------------------------------------------------------------------------------

    if (ProgressDiagnostics):
        print " == >> RCR CELLML MODEL << == "
        
    qCellMLComponent = 1
    pCellMLComponent = 2

    # Create the CellML environment
    CellML = CMISS.CellML()
    CellML.CreateStart(CellMLUserNumber,Region)
    # Number of CellML models
    CellMLModelIndex = [0]*(numberOfTerminalNodes+1)

    # Windkessel Model
    for terminalIdx in range (1,numberOfTerminalNodes+1):
        nodeIdx = coupledNodeNumbers[terminalIdx-1]
        nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumber)
        print('reading model: ' + "./input/CellMLModels/outlet/"+str(terminalIdx)+"/ModelRCR.cellml")
        if (nodeDomain == computationalNodeNumber):
            CellMLModelIndex[terminalIdx] = CellML.ModelImport("./input/CellMLModels/outlet/"+str(terminalIdx)+"/ModelRCR.cellml")
            # known (to OpenCMISS) variables
            CellML.VariableSetAsKnown(CellMLModelIndex[terminalIdx],"Circuit/Qin")
            # to get from the CellML side 
            CellML.VariableSetAsWanted(CellMLModelIndex[terminalIdx],"Circuit/Pout")
    CellML.CreateFinish()

    # Start the creation of CellML <--> OpenCMISS field maps
    CellML.FieldMapsCreateStart()
    
    # ModelIndex
    for terminalIdx in range (1,numberOfTerminalNodes+1):
        nodeIdx = coupledNodeNumbers[terminalIdx-1]
        nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumber)
        if (nodeDomain == computationalNodeNumber):
            # Now we can set up the field variable component <--> CellML model variable mappings.
            # Map the OpenCMISS boundary flow rate values --> CellML
            # Q is component 1 of the DependentField
            CellML.CreateFieldToCellMLMap(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,1,
                                          CMISS.FieldParameterSetTypes.VALUES,CellMLModelIndex[terminalIdx],"Circuit/Qin",CMISS.FieldParameterSetTypes.VALUES)
            # Map the returned pressure values from CellML --> CMISS
            # pCellML is component 1 of the Dependent field U1 variable
            CellML.CreateCellMLToFieldMap(CellMLModelIndex[terminalIdx],"Circuit/Pout",CMISS.FieldParameterSetTypes.VALUES,
                                          DependentFieldNavierStokes,CMISS.FieldVariableTypes.U1,pCellMLComponent,CMISS.FieldParameterSetTypes.VALUES)

    # Finish the creation of CellML <--> OpenCMISS field maps
    CellML.FieldMapsCreateFinish()

    CellMLModelsField = CMISS.Field()
    CellML.ModelsFieldCreateStart(CellMLModelsFieldUserNumber,CellMLModelsField)
    CellML.ModelsFieldCreateFinish()
    
    # Set the models field at boundary nodes
    for terminalIdx in range (1,numberOfTerminalNodes+1):
        nodeIdx = coupledNodeNumbers[terminalIdx-1]
        nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumber)
        if (nodeDomain == computationalNodeNumber):
            #print("Terminal node: " + str(nodeIdx))
            versionIdx = 1
            CellMLModelsField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
                                                     versionIdx,1,nodeIdx,1,CellMLModelIndex[terminalIdx])

    CellMLStateField = CMISS.Field()
    CellML.StateFieldCreateStart(CellMLStateFieldUserNumber,CellMLStateField)
    CellML.StateFieldCreateFinish()

    CellMLParametersField = CMISS.Field()
    CellML.ParametersFieldCreateStart(CellMLParametersFieldUserNumber,CellMLParametersField)
    CellML.ParametersFieldCreateFinish()

    CellMLIntermediateField = CMISS.Field()
    CellML.IntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber,CellMLIntermediateField)
    CellML.IntermediateFieldCreateFinish()

    # Finish the parameter update
    DependentFieldNavierStokes.ParameterSetUpdateStart(CMISS.FieldVariableTypes.U1,CMISS.FieldParameterSetTypes.VALUES)
    DependentFieldNavierStokes.ParameterSetUpdateFinish(CMISS.FieldVariableTypes.U1,CMISS.FieldParameterSetTypes.VALUES)
# DOC-END cellml define field maps

#================================================================================================================================
#  Equations
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> EQUATIONS << == "

# 2nd Equations Set - CHARACTERISTIC
EquationsCharacteristic = CMISS.Equations()
EquationsSetCharacteristic.EquationsCreateStart(EquationsCharacteristic)
EquationsCharacteristic.sparsityType = CMISS.EquationsSparsityTypes.SPARSE
# (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
EquationsCharacteristic.outputType = CMISS.EquationsOutputTypes.NONE
EquationsSetCharacteristic.EquationsCreateFinish()

#------------------

# 3rd Equations Set - NAVIER-STOKES
EquationsNavierStokes = CMISS.Equations()
EquationsSetNavierStokes.EquationsCreateStart(EquationsNavierStokes)
EquationsNavierStokes.sparsityType = CMISS.EquationsSparsityTypes.FULL
EquationsNavierStokes.lumpingType = CMISS.EquationsLumpingTypes.UNLUMPED
# (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
EquationsNavierStokes.outputType = CMISS.EquationsOutputTypes.NONE
EquationsSetNavierStokes.EquationsCreateFinish()

#================================================================================================================================
#  Problems
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> PROBLEM << == "

# Start the creation of a problem.
Problem = CMISS.Problem()
Problem.CreateStart(ProblemUserNumber)
Problem.SpecificationSet(CMISS.ProblemClasses.FLUID_MECHANICS,
                         CMISS.ProblemTypes.NAVIER_STOKES_EQUATION,ProblemSubtype)    
Problem.CreateFinish()

#================================================================================================================================
#  Control Loops
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> PROBLEM CONTROL LOOP << == "
    
'''
   Solver Control Loops

                   L1                                 L2                        L3

1D0D
------

                                                      | 1) 0D Simple subloop   | 1) 0D/CellML DAE Solver
                                                      |                              
    Time Loop, L0  | 1) 1D-0D Iterative Coupling, L1  | 2) 1D NS/C coupling:   | 1) Characteristic Nonlinear Solver
                   |    Convergence Loop (while loop) |    (while loop)        | 2) 1DNavierStokes Transient Solver
                   |
                   | 2) (optional) Simple subloop     | 


'''

# Order of solvers within their respective subloops
SolverCharacteristicUserNumber = 1
SolverNavierStokesUserNumber   = 2
SolverCellmlUserNumber         = 1
if (RCRBoundaries):
   Iterative1d0dControlLoopNumber   = 1
   Simple0DControlLoopNumber        = 1
   Iterative1dControlLoopNumber     = 2
else:
   Iterative1dControlLoopNumber     = 1

# Start the creation of the problem control loop
TimeLoop = CMISS.ControlLoop()
Problem.ControlLoopCreateStart()
Problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE],TimeLoop)
TimeLoop.LabelSet('Time Loop')
TimeLoop.TimesSet(startTime,stopTime,timeIncrement)
TimeLoop.TimeOutputSet(dynamicSolverNavierStokesOutputFrequency)

# Set tolerances for iterative convergence loops
if (RCRBoundaries):
    Iterative1DCouplingLoop = CMISS.ControlLoop()
    Problem.ControlLoopGet([Iterative1d0dControlLoopNumber,Iterative1dControlLoopNumber,
                            CMISS.ControlLoopIdentifiers.NODE],Iterative1DCouplingLoop)
    Iterative1DCouplingLoop.AbsoluteToleranceSet(couplingTolerance1D)
    Iterative1D0DCouplingLoop = CMISS.ControlLoop()
    Problem.ControlLoopGet([Iterative1d0dControlLoopNumber,CMISS.ControlLoopIdentifiers.NODE],
                           Iterative1D0DCouplingLoop)
    Iterative1D0DCouplingLoop.AbsoluteToleranceSet(couplingTolerance1D0D)
else:
    Iterative1DCouplingLoop = CMISS.ControlLoop()
    Problem.ControlLoopGet([Iterative1dControlLoopNumber,CMISS.ControlLoopIdentifiers.NODE],
                           Iterative1DCouplingLoop)
    Iterative1DCouplingLoop.AbsoluteToleranceSet(couplingTolerance1D)

Problem.ControlLoopCreateFinish()

#================================================================================================================================
#  Solvers
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> SOLVERS << == "

# Start the creation of the problem solvers    
DynamicSolverNavierStokes     = CMISS.Solver()
NonlinearSolverNavierStokes   = CMISS.Solver()
LinearSolverNavierStokes      = CMISS.Solver()
NonlinearSolverCharacteristic = CMISS.Solver()
LinearSolverCharacteristic    = CMISS.Solver()

Problem.SolversCreateStart()

# 1st Solver, Simple 0D subloop - CellML
if (RCRBoundaries):
    CellMLSolver = CMISS.Solver()
    Problem.SolverGet([Iterative1d0dControlLoopNumber,Simple0DControlLoopNumber,
                       CMISS.ControlLoopIdentifiers.NODE],SolverDAEUserNumber,CellMLSolver)
    CellMLSolver.OutputTypeSet(cmissSolverOutputType)

# 1st Solver, Iterative 1D subloop - CHARACTERISTIC
if (RCRBoundaries):
    Problem.SolverGet([Iterative1d0dControlLoopNumber,Iterative1dControlLoopNumber,
                       CMISS.ControlLoopIdentifiers.NODE],SolverCharacteristicUserNumber,NonlinearSolverCharacteristic)
else:
    Problem.SolverGet([Iterative1dControlLoopNumber,CMISS.ControlLoopIdentifiers.NODE],
                      SolverCharacteristicUserNumber,NonlinearSolverCharacteristic)
# Set the nonlinear Jacobian type
NonlinearSolverCharacteristic.NewtonJacobianCalculationTypeSet(CMISS.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
NonlinearSolverCharacteristic.OutputTypeSet(nonlinearSolverCharacteristicsOutputType)
# Set the solver settings
NonlinearSolverCharacteristic.NewtonAbsoluteToleranceSet(absoluteToleranceNonlinearCharacteristic)
NonlinearSolverCharacteristic.NewtonSolutionToleranceSet(solutionToleranceNonlinearCharacteristic)
NonlinearSolverCharacteristic.NewtonRelativeToleranceSet(relativeToleranceNonlinearCharacteristic)
# Get the nonlinear linear solver
NonlinearSolverCharacteristic.NewtonLinearSolverGet(LinearSolverCharacteristic)
LinearSolverCharacteristic.OutputTypeSet(linearSolverCharacteristicOutputType)
# Set the solver settings
LinearSolverCharacteristic.LinearTypeSet(CMISS.LinearSolverTypes.ITERATIVE)
LinearSolverCharacteristic.LinearIterativeMaximumIterationsSet(MAXIMUM_ITERATIONS)
LinearSolverCharacteristic.LinearIterativeDivergenceToleranceSet(DIVERGENCE_TOLERANCE)
LinearSolverCharacteristic.LinearIterativeRelativeToleranceSet(relativeToleranceLinearCharacteristic)
LinearSolverCharacteristic.LinearIterativeAbsoluteToleranceSet(absoluteToleranceLinearCharacteristic)
LinearSolverCharacteristic.LinearIterativeGMRESRestartSet(RESTART_VALUE)

#------------------

# 2nd Solver, Iterative 1D subloop - NAVIER-STOKES
if (RCRBoundaries):
    Problem.SolverGet([Iterative1d0dControlLoopNumber,Iterative1dControlLoopNumber,
                       CMISS.ControlLoopIdentifiers.NODE],SolverNavierStokesUserNumber,DynamicSolverNavierStokes)
else:
    Problem.SolverGet([Iterative1dControlLoopNumber,CMISS.ControlLoopIdentifiers.NODE],
                      SolverNavierStokesUserNumber,DynamicSolverNavierStokes)
DynamicSolverNavierStokes.OutputTypeSet(dynamicSolverNavierStokesOutputType)
DynamicSolverNavierStokes.DynamicThetaSet(dynamicSolverNavierStokesTheta)
# Get the dynamic nonlinear solver
DynamicSolverNavierStokes.DynamicNonlinearSolverGet(NonlinearSolverNavierStokes)
# Set the nonlinear Jacobian type
NonlinearSolverNavierStokes.NewtonJacobianCalculationTypeSet(CMISS.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
NonlinearSolverNavierStokes.OutputTypeSet(nonlinearSolverNavierStokesOutputType)

# Set the solver settings
NonlinearSolverNavierStokes.NewtonAbsoluteToleranceSet(absoluteToleranceNonlinearNavierStokes)
NonlinearSolverNavierStokes.NewtonSolutionToleranceSet(solutionToleranceNonlinearNavierStokes)
NonlinearSolverNavierStokes.NewtonRelativeToleranceSet(relativeToleranceNonlinearNavierStokes)
# Get the dynamic nonlinear linear solver
NonlinearSolverNavierStokes.NewtonLinearSolverGet(LinearSolverNavierStokes)
LinearSolverNavierStokes.OutputTypeSet(linearSolverNavierStokesOutputType)
# Set the solver settings
LinearSolverNavierStokes.LinearTypeSet(CMISS.LinearSolverTypes.ITERATIVE)
LinearSolverNavierStokes.LinearIterativeMaximumIterationsSet(MAXIMUM_ITERATIONS)
LinearSolverNavierStokes.LinearIterativeDivergenceToleranceSet(DIVERGENCE_TOLERANCE)
LinearSolverNavierStokes.LinearIterativeRelativeToleranceSet(relativeToleranceLinearNavierStokes)
LinearSolverNavierStokes.LinearIterativeAbsoluteToleranceSet(absoluteToleranceLinearNavierStokes)
LinearSolverNavierStokes.LinearIterativeGMRESRestartSet(RESTART_VALUE)
    
# Finish the creation of the problem solver
Problem.SolversCreateFinish()

#================================================================================================================================
#  Solver Equations
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> SOLVER EQUATIONS << == "

# Start the creation of the problem solver equations
NonlinearSolverCharacteristic = CMISS.Solver()
SolverEquationsCharacteristic = CMISS.SolverEquations()
DynamicSolverNavierStokes     = CMISS.Solver()
SolverEquationsNavierStokes   = CMISS.SolverEquations()

Problem.SolverEquationsCreateStart()

# CellML Solver
if (RCRBoundaries):
    CellMLSolver = CMISS.Solver()
    CellMLEquations = CMISS.CellMLEquations()
    Problem.CellMLEquationsCreateStart()
    Problem.SolverGet([Iterative1d0dControlLoopNumber,Simple0DControlLoopNumber,
                       CMISS.ControlLoopIdentifiers.NODE],SolverDAEUserNumber,CellMLSolver)
    CellMLSolver.CellMLEquationsGet(CellMLEquations)
    # Add in the equations set
    CellMLEquations.CellMLAdd(CellML)    
    Problem.CellMLEquationsCreateFinish()

#------------------

# CHARACTERISTIC solver
if (RCRBoundaries):
    Problem.SolverGet([Iterative1d0dControlLoopNumber,Iterative1dControlLoopNumber,
                       CMISS.ControlLoopIdentifiers.NODE],SolverCharacteristicUserNumber,NonlinearSolverCharacteristic)
else:
    Problem.SolverGet([Iterative1dControlLoopNumber,CMISS.ControlLoopIdentifiers.NODE],
                      SolverCharacteristicUserNumber,NonlinearSolverCharacteristic)
NonlinearSolverCharacteristic.SolverEquationsGet(SolverEquationsCharacteristic)
SolverEquationsCharacteristic.sparsityType = CMISS.SolverEquationsSparsityTypes.SPARSE
# Add in the equations set
EquationsSetCharacteristic = SolverEquationsCharacteristic.EquationsSetAdd(EquationsSetCharacteristic)

#  NAVIER-STOKES solver
if (RCRBoundaries):
    Problem.SolverGet([Iterative1d0dControlLoopNumber,Iterative1dControlLoopNumber,
                       CMISS.ControlLoopIdentifiers.NODE],SolverNavierStokesUserNumber,DynamicSolverNavierStokes)
else:
    Problem.SolverGet([Iterative1dControlLoopNumber,CMISS.ControlLoopIdentifiers.NODE],
                      SolverNavierStokesUserNumber,DynamicSolverNavierStokes)
DynamicSolverNavierStokes.SolverEquationsGet(SolverEquationsNavierStokes)
SolverEquationsNavierStokes.sparsityType = CMISS.SolverEquationsSparsityTypes.SPARSE
# Add in the equations set
EquationsSetNavierStokes = SolverEquationsNavierStokes.EquationsSetAdd(EquationsSetNavierStokes)

# Finish the creation of the problem solver equations
Problem.SolverEquationsCreateFinish()
    
#================================================================================================================================
#  Boundary Conditions
#================================================================================================================================

if (ProgressDiagnostics):
    print " == >> BOUNDARY CONDITIONS << == "

# CHARACTERISTIC
BoundaryConditionsCharacteristic = CMISS.BoundaryConditions()
SolverEquationsCharacteristic.BoundaryConditionsCreateStart(BoundaryConditionsCharacteristic)

# Area-outlet
versionIdx = 1
for terminalIdx in range (1,numberOfTerminalNodes+1):
    nodeNumber = coupledNodeNumbers[terminalIdx-1]
    nodeDomain = Decomposition.NodeDomainGet(nodeNumber,meshComponentNumber)
    if (nodeDomain == computationalNodeNumber):
        if (nonReflecting):
            BoundaryConditionsCharacteristic.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
                                                     versionIdx,1,nodeNumber,2,CMISS.BoundaryConditionsTypes.FIXED_NONREFLECTING,A[nodeNumber][0])
        elif (RCRBoundaries):
            BoundaryConditionsCharacteristic.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
                                                     versionIdx,1,nodeNumber,2,CMISS.BoundaryConditionsTypes.FIXED_CELLML,A[nodeNumber][0])
        else:
            BoundaryConditionsCharacteristic.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
                                                     versionIdx,1,nodeNumber,2,CMISS.BoundaryConditionsTypes.FIXED_OUTLET,A[nodeNumber][0])

SolverEquationsCharacteristic.BoundaryConditionsCreateFinish()

#------------------

# NAVIER-STOKES
BoundaryConditionsNavierStokes = CMISS.BoundaryConditions()
SolverEquationsNavierStokes.BoundaryConditionsCreateStart(BoundaryConditionsNavierStokes)

# Inlet (Flow)
versionIdx = 1
for inputIdx in range (1,numberOfInputNodes+1):
    nodeNumber = inputNodeNumbers[inputIdx-1]
    nodeDomain = Decomposition.NodeDomainGet(nodeNumber,meshComponentNumber)
    if (nodeDomain == computationalNodeNumber):
        BoundaryConditionsNavierStokes.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
                                               versionIdx,1,nodeNumber,1,CMISS.BoundaryConditionsTypes.FIXED_FITTED,Q[inputIdx][0])
# Area-outlet
versionIdx = 1
for terminalIdx in range (1,numberOfTerminalNodes+1):
    nodeNumber = coupledNodeNumbers[terminalIdx-1]
    nodeDomain = Decomposition.NodeDomainGet(nodeNumber,meshComponentNumber)
    if (nodeDomain == computationalNodeNumber):
        if (nonReflecting):
            BoundaryConditionsNavierStokes.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
                                                   versionIdx,1,nodeNumber,2,CMISS.BoundaryConditionsTypes.FIXED_NONREFLECTING,A[nodeNumber][0])
        elif (RCRBoundaries):
            BoundaryConditionsNavierStokes.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
                                                   versionIdx,1,nodeNumber,2,CMISS.BoundaryConditionsTypes.FIXED_CELLML,A[nodeNumber][0])
        else:
            BoundaryConditionsNavierStokes.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
                                                   versionIdx,1,nodeNumber,2,CMISS.BoundaryConditionsTypes.FIXED_OUTLET,A[nodeNumber][0])

# Finish the creation of boundary conditions
SolverEquationsNavierStokes.BoundaryConditionsCreateFinish()

if (CheckTimestepStability):
    QMax = 430.0
    maxTimestep = Utilities1D.GetMaxStableTimestep(elementNodes,QMax,nodeCoordinates,H,E,A0,Rho)
    if (timeIncrement > maxTimestep):
        sys.exit('Timestep size '+str(timeIncrement)+' above maximum allowable size of '+str(maxTimestep)+'. Please reduce step size and re-run')

#================================================================================================================================
#  Run Solvers
#================================================================================================================================
         
# Solve the problem
print "Solving problem..."
start = time.time()
Problem.Solve()
end = time.time()
elapsed = end - start
print "Total Number of Elements = %d " %numberOfElements
print "Calculation Time = %3.4f" %elapsed
print "Problem solved!"

# Remove CellML tmp files
for filename in glob.glob("./tmp.cellml2code*"):
    shutil.rmtree(filename)

#================================================================================================================================
#  Finish Program
#================================================================================================================================
