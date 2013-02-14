#> \file
#> \author Soroush Safaei
#> \brief This is an example program to solve 1D Transient Navier-Stokes over a branch
#>  with coupled 0D lumped models (resistance, RCR) defined in CellML.
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
#> Contributor(s): David Ladd
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
#> OpenCMISS/example/FluidMechanics/NavierStokes/Coupled1DCellML/src/Coupled1DCellMLExample.f90
#<

#================================================================================================================================
#  Guidance
#================================================================================================================================

    #----------------------------------------------------------------------------------------------------------------------------
    # Don't forget to determine these values before running the example:
    #
    # cellmlFlag
    # References
    # TimeStep
    #
    #----------------------------------------------------------------------------------------------------------------------------

#================================================================================================================================
#  Start Program
#================================================================================================================================

# Program Variables
EquationsSetFieldUserNumberNavierStokes = 1337
EquationsSetFieldUserNumberCharacteristic = 1338

DomainUserNumber = 1
CoordinateSystemUserNumber = 1
RegionUserNumber = 2
MeshUserNumber = 3
DecompositionUserNumber = 4
GeometricFieldUserNumber = 5
DependentFieldUserNumber = 6
MaterialsFieldUserNumber = 7
IndependentFieldUserNumber = 8
EquationsSetUserNumberNavierStokes = 9
EquationsSetUserNumberCharacteristic = 10
ProblemUserNumber = 11
CellMLUserNumber = 12
CellMLModelsFieldUserNumber = 13
CellMLStateFieldUserNumber = 14
CellMLIntermediateFieldUserNumber = 15
CellMLParametersFieldUserNumber = 16
MaterialsFieldUserNumberCellML = 17

SolverCharacteristicUserNumber = 1
SolverNavierStokesUserNumber = 2

MaterialsFieldUserNumberMu = 1
MaterialsFieldUserNumberRho = 2
MaterialsFieldUserNumberK = 3
MaterialsFieldUserNumberAs = 4
MaterialsFieldUserNumberRe = 5
MaterialsFieldUserNumberFr = 6
MaterialsFieldUserNumberSt = 7
MaterialsFieldUserNumberA0 = 8
MaterialsFieldUserNumberE = 9
MaterialsFieldUserNumberH0 = 10

#================================================================================================================================
#  Initialise OpenCMISS
#================================================================================================================================

from opencmiss import CMISS
# Add Python bindings directory to PATH
import pdb
import sys,os
import fem_topology
sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

#CMISS.DiagnosticsSetOn(CMISS.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",[""])
#CMISS.ErrorHandlingModeSet(CMISS.ErrorHandlingModes.TRAP_ERROR)
#CMISS.OutputSetOn("Testing")

NumberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
ComputationalNodeNumber = CMISS.ComputationalNodeNumberGet()

#================================================================================================================================
#  Problem Control Panel
#================================================================================================================================

# Reading geometry files
IPNODE_FILE = open("Input/Geometry.ipnode").read()
IPELEM_FILE = open("Input/Geometry.ipelem").read()
IPNODE_LINE = IPNODE_FILE.split()
IPELEM_LINE = IPELEM_FILE.split()

# Set geometry parameters
NumberOfNodesSpace = int(IPNODE_LINE[14].split()[-1]) 
TotalNumberOfElements = int(IPELEM_LINE[14].split()[-1])

# Set branching nodes 
NumberOfBifurcations = 0
NumberOfTrifurcations = 0
bifurcationNodeNumber = []
bifurcationElementNumber = []
trifurcationNodeNumber = []
trifurcationElementNumber = []
bifurcationNodeNumber.append('null')
bifurcationElementNumber.append('null')
trifurcationNodeNumber.append('null')
trifurcationElementNumber.append('null')

# Set flags
cellmlFlag = False
windkesselFlag = False

# Set geometry parameters
NumberOfDimensions = 1
NumberOfNodesFlow = NumberOfNodesSpace
NumberOfNodesArea = NumberOfNodesSpace
TotalNumberOfNodes = NumberOfNodesSpace*3

# Set basis parameters
BasisUserNumberSpace = 1
BasisUserNumberFlow = 2
BasisUserNumberArea = 3
BasisUserNumberMaterial = 4
BasisXiInterpolationSpace = 2
BasisXiInterpolationFlow = 2
BasisXiInterpolationArea = 2

# Set interpolation parameters
BasisXiGaussSpace = 3
BasisXiGaussFlow = 3
BasisXiGaussArea = 3

#================================================================================================================================
#  Reading Mesh
#================================================================================================================================

# User numbers
NumberOfXi = 3
NumberOfComponents = 2
MeshNumberOfComponents = 1
MeshTotalNumberOfElements = 1
BifurcationNumberOfVersions = 3
TrifurcationNumberOfVersions = 4
GeometricFieldNumberOfVariables = 1
GeometricFieldNumberOfComponents = NumberOfXi

# Regions
WorldRegion = fem_topology.femInitialize()
WorldRegion.RegionsCreateStart(RegionUserNumber)
WorldRegion.RegionsCreateFinish(RegionUserNumber)
REGION = WorldRegion.RegionsRegionGet(RegionUserNumber)

# Basis
REGION.BASES.BasesCreateStart(BasisUserNumberSpace)
REGION.BASES.BasisTypeSet(BasisUserNumberSpace,"Quadratic_Lagrange")
REGION.BASES.BasisNumberOfXiCoordinatesSet(BasisUserNumberSpace,NumberOfXi)
REGION.BASES.BasesCreateFinish(BasisUserNumberSpace)
REGION.ReadMesh(BasisUserNumberSpace,MeshUserNumber,
    GeometricFieldUserNumber,"CMISS","Input/Geometry","Input/Geometry")

# Variables    
FieldVariable = 1
ElementNodes = [0,0,0]*(TotalNumberOfElements)
xValues = [0]*(NumberOfNodesSpace)
yValues = [0]*(NumberOfNodesSpace)
zValues = [0]*(NumberOfNodesSpace)

# Extract geometry data
for i in range(1,NumberOfNodesSpace+1):
    xValues[i-1] = REGION.FIELDS.FieldParameterSetNodeValuesGet(GeometricFieldUserNumber,FieldVariable,i,1)
    yValues[i-1] = REGION.FIELDS.FieldParameterSetNodeValuesGet(GeometricFieldUserNumber,FieldVariable,i,2)
    zValues[i-1] = REGION.FIELDS.FieldParameterSetNodeValuesGet(GeometricFieldUserNumber,FieldVariable,i,3)

# Number of bifurcations
for i in range(1,NumberOfNodesSpace+1):
    if (xValues[i-1][1] != 0):
        if (xValues[i-1][3] == 0):
            NumberOfBifurcations = NumberOfBifurcations+1
            bifurcationNodeNumber.append(i)

# Bifurcations elements
for i in range(1,TotalNumberOfElements+1):
    ElementNodes[i-1] = REGION.MESHES.MeshElementsNodesGet(MeshUserNumber,MeshNumberOfComponents,i)
    for j in range(1,NumberOfBifurcations+1):
        if (ElementNodes[i-1][2] == bifurcationNodeNumber[j]):
            bifurcationElementNumber.append(i)

# Number of trifurcations
for i in range(1,NumberOfNodesSpace+1):
    if (xValues[i-1][1] != 0):
        if (xValues[i-1][3] != 0):
            NumberOfTrifurcations = NumberOfTrifurcations+1
            trifurcationNodeNumber.append(i)

# Trifurcation elements
for i in range(1,TotalNumberOfElements+1):
    ElementNodes[i-1] = REGION.MESHES.MeshElementsNodesGet(MeshUserNumber,MeshNumberOfComponents,i)
    for j in range(1,NumberOfTrifurcations+1):
        if (ElementNodes[i-1][2] == trifurcationNodeNumber[j]):
            trifurcationElementNumber.append(i)
            
#================================================================================================================================
#  Set Initial & Default Values
#================================================================================================================================

# Set material parameters
Pi = 3.141593                               # Pi
MU_PARAM = 0.0035                           # Mu (pa.s)
RHO_PARAM = 1050.0                          # Rho (kg/m3)
A0_PARAM = [0]*(NumberOfNodesSpace+1)       # Area (m2)
E_PARAM = [0]*(NumberOfNodesSpace+1)        # Elasticity (pa)
H0_PARAM = [0]*(NumberOfNodesSpace+1)       # Thickness (m)

# Set reference values
Qs = 100.0e-6                               # Flow (m3/s)
As = 100.0e-6                               # Area (m2)
Xs = 1.0                                    # Distance (m)
Ts = 1.0                                    # Time (s)
K = 4.0/3.0                                 # Parabolic flow section
Re = 8.0*Pi*MU_PARAM*Xs/(RHO_PARAM*Qs)      # Reynolds number
Fr = (As**2.5/Qs**2.0)/(2.0*RHO_PARAM)      # Froude number
St = (As*Xs)/(Ts*Qs)                        # Strouhal number

# Reading Initial file
IPNODE_FILE = open("Input/Material").read()
IPNODE_LINE = IPNODE_FILE.split()

# Set Reference Area, Elasticity & Thickness
for i in range(1,NumberOfNodesSpace+1):
    A0_PARAM[i] = float(IPNODE_LINE[2+3*(i-1)].split()[-1]) 
for i in range(1,NumberOfNodesSpace+1):
    E_PARAM[i] = float(IPNODE_LINE[(2+3*(NumberOfNodesSpace-1))+3*i].split()[-1]) 
for i in range(1,NumberOfNodesSpace+1):
    H0_PARAM[i] = float(IPNODE_LINE[(2*(2+3*(NumberOfNodesSpace-1)))+1+3*i].split()[-1]) 

# Set initial conditions
Q = [0]*(NumberOfNodesSpace+1)
A = [0]*(NumberOfNodesSpace+1)
Q[0] = 1.0
for i in range(1,NumberOfNodesSpace+1):
    Q[i] = 0.0
    A[i] = A0_PARAM[i]/As

# Reading Terminal file
IPNODE_FILE = open("Input/Terminal").read()
IPNODE_LINE = IPNODE_FILE.split()

# Set terminal nodes
NumberOfTerminalNodes = int(IPNODE_LINE[1].split()[-1])

# Set terminal conditions
coupledNodeNumber = [0]*(NumberOfTerminalNodes+1)
resistanceProximal = [0]*(NumberOfTerminalNodes+1)
resistanceDistal = [0]*(NumberOfTerminalNodes+1)
capacitance = [0]*(NumberOfTerminalNodes+1)
Ae = [0]*(NumberOfTerminalNodes+1)

# Set lumped model parameters
for i in range(1,NumberOfTerminalNodes+1):
    coupledNodeNumber[i] = int(IPNODE_LINE[1+i].split()[-1])
for i in range(1,NumberOfTerminalNodes+1):
    Ae[i] = A0_PARAM[coupledNodeNumber[i]]/As
for i in range(1,NumberOfTerminalNodes+1):
    resistanceProximal[i] = float(IPNODE_LINE[1+NumberOfTerminalNodes+3*i].split()[-1])
for i in range(1,NumberOfTerminalNodes+1):
    resistanceDistal[i] = float(IPNODE_LINE[(1+NumberOfTerminalNodes+3*NumberOfTerminalNodes)+3*i].split()[-1])
for i in range(1,NumberOfTerminalNodes+1):
    capacitance[i] = float(IPNODE_LINE[(1+NumberOfTerminalNodes+6*NumberOfTerminalNodes)+3*i].split()[-1])

# Set output parameters
# (NONE/PROGRESS/TIMING/SOLVER/MATRIX)
DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE = CMISS.SolverOutputTypes.NONE
NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE = CMISS.SolverOutputTypes.NONE
LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE = CMISS.SolverOutputTypes.NONE
# (NONE/TIMING/SOLVER/MATRIX)
EQUATIONS_NAVIER_STOKES_OUTPUT = CMISS.SolverOutputTypes.NONE

# Set time parameters
DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME = 0.0
DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME = 1.0
DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT = 0.01
DYNAMIC_SOLVER_NAVIER_STOKES_THETA = [0.5]

# Set result output parameters
DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY = 1

# Set solver parameters
LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG = False
RELATIVE_TOLERANCE = 1.0E-10  # default: 1.0E-05
ABSOLUTE_TOLERANCE = 1.0E-10  # default: 1.0E-10
DIVERGENCE_TOLERANCE = 1.0E20 # default: 1.0E5
MAXIMUM_ITERATIONS = 100000   # default: 100000
RESTART_VALUE = 3000          # default: 30
LINESEARCH_ALPHA = 1.0

# Set CellML flag
if (cellmlFlag):
    # New equations set type to store p values in the Equations Set Field
    EquationsSetSubtype = CMISS.EquationsSetSubtypes.Coupled1D0D_NAVIER_STOKES
    # Characteristic (nodal/characteristic) solver remains the same
    EquationsSetCharacteristicSubtype = CMISS.EquationsSetSubtypes.Coupled1D0D_CHARACTERISTIC
    # New problem type to execute the 1D-0D coupling subloop at each timestep
    ProblemSubtype = CMISS.ProblemSubTypes.Coupled1D0D_NAVIER_STOKES
else: 
    EquationsSetSubtype = CMISS.EquationsSetSubtypes.OneDTRANSIENT_NAVIER_STOKES
    EquationsSetCharacteristicSubtype = CMISS.EquationsSetSubtypes.STATIC_CHARACTERISTIC
    ProblemSubtype = CMISS.ProblemSubTypes.OneDTRANSIENT_NAVIER_STOKES

#================================================================================================================================
#  Coordinate System
#================================================================================================================================

# Start the creation of a new RC coordinate system
CoordinateSystem = CMISS.CoordinateSystem()
CoordinateSystem.CreateStart(CoordinateSystemUserNumber)
# Set the coordinate system dimension
CoordinateSystem.dimension = 3
# Finish the creation of the coordinate system
CoordinateSystem.CreateFinish()

#================================================================================================================================
#  Region
#================================================================================================================================

# Start the creation of a new region
Region = CMISS.Region()
# Initialise and Create Regions
Region.CreateStart(RegionUserNumber,CMISS.WorldRegion)
# Set the region label
Region.label = "OpenCMISS"
# Set the region coordinate system as defined above
Region.coordinateSystem = CoordinateSystem
# Finish the creation of the region
Region.CreateFinish()

#================================================================================================================================
#  Bases
#================================================================================================================================

# Start the creation of SPACE bases
MeshNumberOfComponents = 2
BasisSpace = CMISS.Basis()
BasisSpace.CreateStart(BasisUserNumberSpace)
# Set the basis type (Lagrange/Simplex)
BasisSpace.type = CMISS.BasisTypes.LAGRANGE_HERMITE_TP
# Set the basis xi number
BasisSpace.numberOfXi = NumberOfDimensions
# Set the basis xi interpolation and number of Gauss points
BasisSpace.interpolationXi = [CMISS.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]
BasisSpace.quadratureNumberOfGaussXi = [BasisXiGaussSpace]
# Finish the creation of the basis
BasisSpace.CreateFinish()

# Start the creation of FLOW basis
if (BasisXiInterpolationFlow == BasisXiInterpolationSpace):
    BasisFlow = BasisSpace
    
# Start the creation of AREA basis
if (BasisXiInterpolationArea == BasisXiInterpolationSpace):
    BasisArea = BasisSpace

# Start the creation of MATERIAL bases
BasisMaterial = CMISS.Basis()
BasisMaterial.CreateStart(BasisUserNumberMaterial)
# Set the basis type (Lagrange/Simplex)
BasisMaterial.type = CMISS.BasisTypes.LAGRANGE_HERMITE_TP
# Set the basis xi number
BasisMaterial.numberOfXi = NumberOfDimensions
# Set the basis xi interpolation and number of Gauss points
BasisMaterial.interpolationXi = [CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]
BasisMaterial.quadratureNumberOfGaussXi = [BasisXiGaussSpace]
# Finish the creation of the basis
BasisMaterial.CreateFinish()

#================================================================================================================================
#  Mesh
#================================================================================================================================
  
# Start the creation of mesh nodes
Nodes = CMISS.Nodes()
Mesh = CMISS.Mesh()
Nodes.CreateStart(Region,TotalNumberOfNodes)
Nodes.CreateFinish()
# Start the creation of the mesh
Mesh.CreateStart(MeshUserNumber,Region,NumberOfDimensions)
# Set number of mesh elements
Mesh.NumberOfElementsSet(TotalNumberOfElements)
# Set number of mesh components
Mesh.NumberOfComponentsSet(MeshNumberOfComponents)
# Specify spatial mesh component
MeshElementsSpace = CMISS.MeshElements()
MeshElementsFlow = CMISS.MeshElements()
MeshElementsArea = CMISS.MeshElements()
MeshElementsMaterial = CMISS.MeshElements()
MeshComponentNumberSpace = 1
MeshComponentNumberFlow = 1
MeshComponentNumberArea = 1
MeshComponentNumberMaterial = 2

# Specify space mesh component
MeshElementsSpace.CreateStart(Mesh,MeshComponentNumberSpace,BasisSpace)
for i in range(1,TotalNumberOfElements+1):
    MeshElementsSpace.NodesSet(i,ElementNodes[i-1])
for i in range(1,NumberOfBifurcations+1):
    MeshElementsSpace.LocalElementNodeVersionSet(bifurcationElementNumber[i]+0,1,1,3) 
    MeshElementsSpace.LocalElementNodeVersionSet(bifurcationElementNumber[i]+1,2,1,1) 
    MeshElementsSpace.LocalElementNodeVersionSet(bifurcationElementNumber[i]+2,3,1,1) 
for i in range(1,NumberOfTrifurcations+1):
    MeshElementsSpace.LocalElementNodeVersionSet(trifurcationElementNumber[i]+0,1,1,3) 
    MeshElementsSpace.LocalElementNodeVersionSet(trifurcationElementNumber[i]+1,2,1,1) 
    MeshElementsSpace.LocalElementNodeVersionSet(trifurcationElementNumber[i]+2,3,1,1)
    MeshElementsSpace.LocalElementNodeVersionSet(trifurcationElementNumber[i]+3,4,1,1) 
    #(globalElementNumber,versionNumber,derivativeNumber,localElementNodeNumber)
MeshElementsSpace.CreateFinish()                        

    #   Arteries
    #                                         l femoral---l popl
    #                                        /
    #                                       /
    #                            l ext iliac
    #                           /           \
    #                          /             \
    #               l com iliac               l pro fem
    #              /           \
    #             /             \
    #            /               l int iliac
    #   abd Aorta             
    #            \               r int iliac
    #             \             /
    #              \           /
    #               r com iliac               r pro fem
    #                          \             /
    #                           \           /
    #                            r ext iliac
    #                                       \
    #                                        \
    #                                         r femoral---r popl

# Specify flow mesh component
if (BasisXiInterpolationFlow == BasisXiInterpolationSpace):
    MeshElementsFlow = MeshElementsSpace
# Specify area mesh component
if (BasisXiInterpolationArea == BasisXiInterpolationSpace):
    MeshElementsArea = MeshElementsSpace

# Specify material mesh component
MeshElementsMaterial.CreateStart(Mesh,MeshComponentNumberMaterial,BasisMaterial)
for i in range(1,TotalNumberOfElements+1):
    MeshElementsMaterial.NodesSet(i,[ElementNodes[i-1][0],ElementNodes[i-1][2]])
for i in range(1,NumberOfBifurcations+1):
    MeshElementsMaterial.LocalElementNodeVersionSet(bifurcationElementNumber[i]+0,1,1,2) 
    MeshElementsMaterial.LocalElementNodeVersionSet(bifurcationElementNumber[i]+1,2,1,1) 
    MeshElementsMaterial.LocalElementNodeVersionSet(bifurcationElementNumber[i]+2,3,1,1) 
for i in range(1,NumberOfTrifurcations+1):
    MeshElementsMaterial.LocalElementNodeVersionSet(trifurcationElementNumber[i]+0,1,1,2) 
    MeshElementsMaterial.LocalElementNodeVersionSet(trifurcationElementNumber[i]+1,2,1,1) 
    MeshElementsMaterial.LocalElementNodeVersionSet(trifurcationElementNumber[i]+2,3,1,1)
    MeshElementsMaterial.LocalElementNodeVersionSet(trifurcationElementNumber[i]+3,4,1,1) 
    #(globalElementNumber,versionNumber,derivativeNumber,localElementNodeNumber)
MeshElementsMaterial.CreateFinish()

# Finish the creation of the mesh
Mesh.CreateFinish()

#================================================================================================================================
#  Decomposition
#================================================================================================================================

# Create a decomposition
Decomposition = CMISS.Decomposition()
Decomposition.CreateStart(DecompositionUserNumber,Mesh)
# Set the decomposition to be a general decomposition with the specified number of domains
Decomposition.TypeSet = CMISS.DecompositionTypes.CALCULATED
Decomposition.NumberOfDomainsSet = DomainUserNumber
# Finish the decomposition
Decomposition.CreateFinish()

#================================================================================================================================
#  Geometric Field
#================================================================================================================================

# Start to create a default (geometric) field on the region
GeometricField = CMISS.Field()
GeometricField.CreateStart(GeometricFieldUserNumber,Region)
# Set the geometric field number of variables (U)
GeometricField.NumberOfVariablesSet(1)
# Set the field label
GeometricField.VariableLabelSet(CMISS.FieldVariableTypes.U,'Coordinates')
# Set the field type
GeometricField.TypeSet = CMISS.FieldTypes.GEOMETRIC
# Set the decomposition to use
GeometricField.meshDecomposition = Decomposition
GeometricField.ScalingTypeSet = CMISS.FieldScalingTypes.NONE
# Set the mesh component to be used by the field components.
for ComponentNumber in range(1,CoordinateSystem.dimension+1):
    GeometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,ComponentNumber,MeshComponentNumberSpace)
# Finish creating the field
GeometricField.CreateFinish()

# Initialise the equations set geometric field variables
for UserNodeNumber in range(1,NumberOfNodesSpace+1):
    GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        1,1,UserNodeNumber,1,xValues[UserNodeNumber-1][0])
    GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        1,1,UserNodeNumber,2,yValues[UserNodeNumber-1][0])
    GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        1,1,UserNodeNumber,3,zValues[UserNodeNumber-1][0])
for i in range (1,NumberOfBifurcations+1):
    for VersionIdx in range(2,BifurcationNumberOfVersions+1):
        GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            VersionIdx,1,bifurcationNodeNumber[i],1,xValues[bifurcationNodeNumber[i]-1][VersionIdx-1])
        GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            VersionIdx,1,bifurcationNodeNumber[i],2,yValues[bifurcationNodeNumber[i]-1][VersionIdx-1])
        GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            VersionIdx,1,bifurcationNodeNumber[i],3,zValues[bifurcationNodeNumber[i]-1][VersionIdx-1])
for i in range (1,NumberOfTrifurcations+1):
    for VersionIdx in range(2,TrifurcationNumberOfVersions+1):
        GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            VersionIdx,1,trifurcationNodeNumber[i],1,xValues[trifurcationNodeNumber[i]-1][VersionIdx-1])
        GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            VersionIdx,1,trifurcationNodeNumber[i],2,yValues[trifurcationNodeNumber[i]-1][VersionIdx-1])
        GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            VersionIdx,1,trifurcationNodeNumber[i],3,zValues[trifurcationNodeNumber[i]-1][VersionIdx-1])
            #(versionNumber,derivativeNumber,userNodeNumber,componentNumber,value)
                        
#================================================================================================================================
#  Equations Sets
#================================================================================================================================

# Create the equations set for Navier-Stokes
EquationsSetNavierStokes = CMISS.EquationsSet()
EquationsSetFieldNavierStokes = CMISS.Field()
# Set the equations set to be a dynamic Navier-Stokes problem
EquationsSetNavierStokes.CreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField,
    CMISS.EquationsSetClasses.FLUID_MECHANICS,CMISS.EquationsSetTypes.NAVIER_STOKES_EQUATION,
    EquationsSetSubtype,EquationsSetFieldUserNumberNavierStokes,EquationsSetFieldNavierStokes)
# Finish creating the equations set
EquationsSetNavierStokes.CreateFinish()

# Create the equations set for Characteristic
EquationsSetCharacteristic = CMISS.EquationsSet()
EquationsSetFieldCharacteristic = CMISS.Field()
# Set the equations set to be a static Nonlinear problem
EquationsSetCharacteristic.CreateStart(EquationsSetUserNumberCharacteristic,Region,GeometricField,
    CMISS.EquationsSetClasses.FLUID_MECHANICS,CMISS.EquationsSetTypes.CHARACTERISTIC_EQUATION,
    EquationsSetCharacteristicSubtype,EquationsSetFieldUserNumberCharacteristic,EquationsSetFieldCharacteristic)
# Finish creating the equations set
EquationsSetCharacteristic.CreateFinish()

#================================================================================================================================
#  Dependent Field
#================================================================================================================================

# Create the equations set dependent field variables for Navier-Stokes
DependentFieldNavierStokes = CMISS.Field()
EquationsSetCharacteristic.DependentCreateStart(DependentFieldUserNumber,DependentFieldNavierStokes)
# Set the field label
DependentFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.U,'General')
DependentFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.DELUDELN,'Derivatives')
DependentFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.V,'Characteristics')
DependentFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.U1,'calculated pressure')

# Set the mesh component to be used by the field components.
COMPONENT_NUMBER = 1 # Flow
DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,COMPONENT_NUMBER,  
     MeshComponentNumberFlow)
DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN,COMPONENT_NUMBER,  
     MeshComponentNumberFlow)
     
COMPONENT_NUMBER = 2 # Area
DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,COMPONENT_NUMBER,  
     MeshComponentNumberArea)
DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN,COMPONENT_NUMBER,  
     MeshComponentNumberArea)
       
for COMPONENT_NUMBER in range(1,3): # W(Characteristics)
    DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.V,COMPONENT_NUMBER,  
        MeshComponentNumberSpace)

if (cellmlFlag):
    COMPONENT_NUMBER = 1  # pCellML
    DependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U1,COMPONENT_NUMBER,  
        MeshComponentNumberSpace)
        
# Finish the equations set dependent field variables
EquationsSetCharacteristic.DependentCreateFinish()

# Create the equations set dependent field variables for static nonlinear Characteristic solver
EquationsSetNavierStokes.DependentCreateStart(DependentFieldUserNumber,DependentFieldNavierStokes)
# Finish the equations set dependent field variables
EquationsSetNavierStokes.DependentCreateFinish()

# Initialise the equations set dependent field variables
# (versionNumber,derivativeNumber,userNodeNumber,componentNumber,value)
for nodeIdx in range (1,NumberOfNodesSpace+1):
    DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
        CMISS.FieldParameterSetTypes.VALUES,1,1,nodeIdx,1,Q[nodeIdx])
    DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
        CMISS.FieldParameterSetTypes.VALUES,1,1,nodeIdx,2,A[nodeIdx])
for i in range (1,NumberOfBifurcations+1):
    DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
        CMISS.FieldParameterSetTypes.VALUES,2,1,bifurcationNodeNumber[i],2,A[bifurcationNodeNumber[i]+1])
    DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
        CMISS.FieldParameterSetTypes.VALUES,3,1,bifurcationNodeNumber[i],2,A[bifurcationNodeNumber[i]+3])
for i in range (1,NumberOfTrifurcations+1):
    DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
        CMISS.FieldParameterSetTypes.VALUES,2,1,trifurcationNodeNumber[i],2,A[trifurcationNodeNumber[i]+1])
    DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
        CMISS.FieldParameterSetTypes.VALUES,3,1,trifurcationNodeNumber[i],2,A[trifurcationNodeNumber[i]+3])
    DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
        CMISS.FieldParameterSetTypes.VALUES,4,1,trifurcationNodeNumber[i],2,A[trifurcationNodeNumber[i]+5])
for i in range (1,NumberOfBifurcations+1):
    for componentIdx in range(1,NumberOfComponents+1):
        for versionIdx in range(1,BifurcationNumberOfVersions+1):
            DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
                CMISS.FieldParameterSetTypes.VALUES,versionIdx,1,bifurcationNodeNumber[i],componentIdx,0.0)
for i in range (1,NumberOfTrifurcations+1):
    for componentIdx in range(1,NumberOfComponents+1):
        for versionIdx in range(1,TrifurcationNumberOfVersions+1):
            DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
                CMISS.FieldParameterSetTypes.VALUES,versionIdx,1,trifurcationNodeNumber[i],componentIdx,0.0)
                
#================================================================================================================================
#  Materials Field
#================================================================================================================================

# Create the equations set materials field variables for Navier-Stokes
MaterialsFieldNavierStokes = CMISS.Field()
EquationsSetNavierStokes.MaterialsCreateStart(MaterialsFieldUserNumber,MaterialsFieldNavierStokes)
# Set the field label
MaterialsFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.U,'Materials')
# Set the mesh component to be used by the field components.
for ComponentNumber in range(8,11):
    MaterialsFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,ComponentNumber,MeshComponentNumberMaterial)
# Finish the equations set materials field variables
EquationsSetNavierStokes.MaterialsCreateFinish()

# Create the equations set materials field variables for Characteristic
EquationsSetCharacteristic.MaterialsCreateStart(MaterialsFieldUserNumber,MaterialsFieldNavierStokes)
# Finish the equations set materials field variables
EquationsSetCharacteristic.MaterialsCreateFinish()

# Initialise the equations set materials field variables
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberMu,MU_PARAM)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberRho,RHO_PARAM)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberK,K)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberAs,As)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberRe,Re)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberFr,Fr)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
     MaterialsFieldUserNumberSt,St)

# Initialise A0,E,H0 Parameters
for i in range(1,NumberOfNodesSpace+1,2):
    MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        1,1,i,MaterialsFieldUserNumberA0,A0_PARAM[i])
    MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        1,1,i,MaterialsFieldUserNumberE,E_PARAM[i])
    MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        1,1,i,MaterialsFieldUserNumberH0,H0_PARAM[i])
for i in range(1,NumberOfBifurcations+1):
    for VersionIdx in range(2,BifurcationNumberOfVersions+1):
        MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            VersionIdx,1,bifurcationNodeNumber[i],MaterialsFieldUserNumberA0,A0_PARAM[bifurcationNodeNumber[i]+
                2*VersionIdx-2])
        MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            VersionIdx,1,bifurcationNodeNumber[i],MaterialsFieldUserNumberE,E_PARAM[bifurcationNodeNumber[i]+
                2*VersionIdx-2])
        MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            VersionIdx,1,bifurcationNodeNumber[i],MaterialsFieldUserNumberH0,H0_PARAM[bifurcationNodeNumber[i]+
                2*VersionIdx-2])
for i in range(1,NumberOfTrifurcations+1):
    for VersionIdx in range(2,TrifurcationNumberOfVersions+1):
        MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            VersionIdx,1,trifurcationNodeNumber[i],MaterialsFieldUserNumberA0,A0_PARAM[trifurcationNodeNumber[i]+
                2*VersionIdx-2])
        MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            VersionIdx,1,trifurcationNodeNumber[i],MaterialsFieldUserNumberE,E_PARAM[trifurcationNodeNumber[i]+
                2*VersionIdx-2])
        MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            VersionIdx,1,trifurcationNodeNumber[i],MaterialsFieldUserNumberH0,H0_PARAM[trifurcationNodeNumber[i]+
                2*VersionIdx-2])

if (cellmlFlag):

    # CellML Materials Field parameters setup
    #--------------------------------------------------------------------------------------------------------------------
    #  pVesselWall: a constant expressing physical features of the vessel wall
    #  pExternal: pressure external to the vessel
    #
    # WARNING: Do not change these component index values- they are what OpenCMISS uses to identify the variables
    pVesselWallComponent = 1
    pExternalComponent = 2
    #-------------------------------------------------------------------------------------------------------------------  

    # User defined variable values (feel free to change these- already initialised to 0 otherwise)
    pVesselWall = 0.0
    pExternal = 0.0
    # Set values at coupled node 
    for i in range (1,NumberOfTerminalNodes+1):
        # (versionNumber,derivativeNumber,userNodeNumber,componentNumber,value)
        MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
            CMISS.FieldParameterSetTypes.VALUES,1,1,coupledNodeNumber[i],pVesselWallComponent,pVesselWall)
        MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
            CMISS.FieldParameterSetTypes.VALUES,1,1,coupledNodeNumber[i],pExternalComponent,pExternal)

#================================================================================================================================
# Independent Field - Characteristic wave normal direction
#================================================================================================================================

# Create the equations set independent field variables for Characteristic Solver
IndependentFieldNavierStokes = CMISS.Field()
EquationsSetCharacteristic.IndependentCreateStart(IndependentFieldUserNumber,IndependentFieldNavierStokes)
IndependentFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.U,'Normal wave direction')
# Set the mesh component to be used by the field components.
IndependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,1,MeshComponentNumberSpace)
IndependentFieldNavierStokes.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,2,MeshComponentNumberSpace)
# Finish the equations set independent field variables
EquationsSetCharacteristic.IndependentCreateFinish()

# Normal Wave Direction for branches
for i in range (1,NumberOfBifurcations+1):
    nodeIdx = bifurcationNodeNumber[i]
    # Incoming(parent)
    IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
        1,1,nodeIdx,1,1.0)
    # Outgoing(daughters)
    IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
        2,1,nodeIdx,2,-1.0)
    IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
        3,1,nodeIdx,2,-1.0)
        
for i in range (1,NumberOfTrifurcations+1):
    nodeIdx = trifurcationNodeNumber[i]
    # Incoming(parent)
    IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
        1,1,nodeIdx,1,1.0)
    # Outgoing(daughters)
    IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
        2,1,nodeIdx,2,-1.0)
    IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
        3,1,nodeIdx,2,-1.0)
    IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
        4,1,nodeIdx,2,-1.0)

# Normal Wave Direction for terminal
if (cellmlFlag):
    for i in range (1,NumberOfTerminalNodes+1):
        nodeIdx = coupledNodeNumber[i]
        # Incoming normals
        IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
             1,1,nodeIdx,1,1.0)
        # Outgoing normals
        IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
            1,1,nodeIdx,2,-1.0)

# Create the equations set independent field variables for Navier-Stokes Solver
EquationsSetNavierStokes.IndependentCreateStart(IndependentFieldUserNumber,IndependentFieldNavierStokes)
# Finish the equations set independent field variables
EquationsSetNavierStokes.IndependentCreateFinish()

#================================================================================================================================
#  CellML Model Maps
#================================================================================================================================

if (cellmlFlag):

    #----------------------------------------------------------------------------------------------------------------------------
    # Description
    #----------------------------------------------------------------------------------------------------------------------------
    # A CellML OD model is used to provide the impedance from the downstream vascular bed beyond the termination
    # point of the 1D model. This is iteratively coupled with the the 1D solver. In the case of a simple resistance
    # model, P=RQ, which is analogous to Ohm's law: V=IR. A variable map copies the guess for the FlowRate, Q at 
    # the boundary from the OpenCMISS Dependent Field to the CellML equation, which then returns presssure, P.
    # The initial guess value for Q is taken from the previous time step or is 0 for t=0.  In OpenCMISS this P value is 
    # then used to compute a new Area value based on the P-A relationship and the Riemann variable W_2, which gives a
    # new value for Q until the values for Q and P converge within tolerance of the previous value.
    #----------------------------------------------------------------------------------------------------------------------------

    pCellMLComponent = 1

    # Windkessel Model
    if (windkesselFlag):
        # Create the CellML environment
        CellML = CMISS.CellML()
        CellML.CreateStart(CellMLUserNumber,Region)

        # Import an RCR windkessel model
        CellML.ModelImport("windkessel.xml",WindkesselModelIndex)    

        # - known (to OpenCMISS) variables 
        CellML.VariableSetAsKnown(WindkesselModelIndex,"equations/Q")
        CellML.VariableSetAsKnown(WindkesselModelIndex,"equations/R_p")
        CellML.VariableSetAsKnown(WindkesselModelIndex,"equations/R_d")
        CellML.VariableSetAsKnown(WindkesselModelIndex,"equations/C")
        CellML.VariableSetAsKnown(WindkesselModelIndex,"equations/t")
        # - to get from the CellML side 
        CellML.VariableSetAsWanted(WindkesselModelIndex,"equations/P")

        CellML.CreateFinish()

        # Start the creation of CellML <--> OpenCMISS field maps
        CellML.FieldMapsCreateStart()
        # Now we can set up the field variable component <--> CellML model variable mappings.

        # Map the OpenCMISS boundary flow rate values --> CellML
        # Q is component 1 of the DependentField
        CellML.CreateFieldToCellMLMap(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,1,
            CMISS.FieldParameterSetTypes.VALUES,WindkesselModelIndex,"equations/Q",CMISS.FieldParameterSetTypes.VALUES)
        CellML.CreateFieldToCellMLMap(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,1,
            CMISS.FieldParameterSetTypes.VALUES,WindkesselModelIndex,"equations/t",CMISS.FieldParameterSetTypes.VALUES)
        # Map the returned pressure values from CellML --> CMISS
        # pCellML is component 1 of the Dependent field U1 variable
        CellML.CreateCellMLToFieldMap(WindkesselModelIndex,"equations/Pressure",CMISS.FieldParameterSetTypes.VALUES,
            DependentFieldNavierStokes,CMISS.FieldVariableTypes.U1,pCellMLComponent,CMISS.FieldParameterSetTypes.VALUES)

        # Finish the creation of CellML <--> OpenCMISS field maps
        CellML.FieldMapsCreateFinish()

        # Create the CellML models field --- only 1 model here
        CellMLModelsField = CMISS.Field()
        CellML.ModelsFieldCreateStart(CellMLModelsFieldUserNumber,CellMLModelsField)
        CellML.ModelsFieldCreateFinish()

        # Create the CellML parameters field --- will be the Resistance and Flow rate
        CellMLParametersField = CMISS.Field()
        CellML.ParametersFieldCreateStart(CellMLParametersFieldUserNumber,CellMLParametersField)
        CellML.ParametersFieldCreateFinish()

        # Create the CellML intermediate field --- will be the pressure value returned from CellML to be used for 
        # recalculation of the incoming Riemann variable W(2)
        CellMLIntermediateField = CMISS.Field()
        CellML.IntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber,CellMLIntermediateField)
        CellML.IntermediateFieldCreateFinish()

        # Initialise pCellML (and previous pCellML coupling iteration values) values to 0 at the outlet nodes
        pCellML = 0.0
        pPrevious = 0.0
        for i in range (1,NumberOfTerminalNodes+1):
            DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U1,CMISS.FieldParameterSetTypes.VALUES,
                1,1,coupledNodeNumber[i],pCellMLComponent,pCellML)
                #(versionNumber,derivativeNumber,userNodeNumber,componentNumber,value)

        DependentFieldNavierStokes.ParameterSetCreate(CMISS.FieldVariableTypes.U1,CMISS.FieldParameterSetTypes.PREVIOUS_VALUES)

        for i in range (1,NumberOfTerminalNodes+1):
            DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U1,
                CMISS.FieldParameterSetTypes.PREVIOUS_VALUES,1,1,coupledNodeNumber[i],pCellMLComponent,pPrevious)

    # Resistance Model
    else: 
        # Create the CellML environment
        CellML = CMISS.CellML()
        CellML.CreateStart(CellMLUserNumber,Region)

        # Import a simple resistance model P = RQ, analogous to V = IR
        CellML.ModelImport("resistance.xml")
        ResistanceModelIndex=1
        # - known (to OpenCMISS) variables 
        CellML.VariableSetAsKnown(ResistanceModelIndex,"equations/FlowRate")
        CellML.VariableSetAsKnown(ResistanceModelIndex,"equations/Resistance")
        # - to get from the CellML side 
        CellML.VariableSetAsWanted(ResistanceModelIndex,"equations/Pressure")

        CellML.CreateFinish()

        # Start the creation of CellML <--> OpenCMISS field maps
        CellML.FieldMapsCreateStart()
        # Now we can set up the field variable component <--> CellML model variable mappings.

        # Map the OpenCMISS boundary flow rate values --> CellML
        # Q is component 1 of the DependentField
        CellML.CreateFieldToCellMLMap(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,1, 
            CMISS.FieldParameterSetTypes.VALUES,ResistanceModelIndex,"equations/FlowRate",CMISS.FieldParameterSetTypes.VALUES)
        # Map the returned pressure values from CellML --> CMISS
        # pCellML is component 2 of the Dependent field V variable
        CellML.CreateCellMLToFieldMap(ResistanceModelIndex,"equations/Pressure",CMISS.FieldParameterSetTypes.VALUES, 
            DependentFieldNavierStokes,CMISS.FieldVariableTypes.U1,pCellMLComponent,CMISS.FieldParameterSetTypes.VALUES)

        # Finish the creation of CellML <--> OpenCMISS field maps
        CellML.FieldMapsCreateFinish()

        # Create the CellML models field --- only 1 model here
        CellMLModelsField = CMISS.Field()
        CellML.ModelsFieldCreateStart(CellMLModelsFieldUserNumber,CellMLModelsField)
        CellML.ModelsFieldCreateFinish()

        # Create the CellML parameters field --- will be the Resistance and Flow rate
        CellMLParametersField = CMISS.Field()
        CellML.ParametersFieldCreateStart(CellMLParametersFieldUserNumber,CellMLParametersField)
        CellML.ParametersFieldCreateFinish()

        # Create the CellML intermediate field --- will be the pressure value returned from CellML to be used for 
        # recalculation of the incoming Riemann variable W(2)
        CellMLIntermediateField = CMISS.Field()
        CellML.IntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber,CellMLIntermediateField)
        CellML.IntermediateFieldCreateFinish()

        # Initialise pCellML (and previous pCellML coupling iteration values) values to 0 at the outlet nodes
        pCellML = 0.0
        pPrevious = 0.0
        for i in range (1,NumberOfTerminalNodes+1):
            # (field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value)
            DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U1,CMISS.FieldParameterSetTypes.VALUES,
                1,1,coupledNodeNumber[i],pCellMLComponent,pCellML)

        DependentFieldNavierStokes.ParameterSetCreate(CMISS.FieldVariableTypes.U1,CMISS.FieldParameterSetTypes.PREVIOUS_VALUES)

        for i in range (1,NumberOfTerminalNodes+1):
            DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U1, 
                CMISS.FieldParameterSetTypes.PREVIOUS_VALUES,1,1,coupledNodeNumber[i],pCellMLComponent,pPrevious)

#================================================================================================================================
#  Equations
#================================================================================================================================

# 1st Equations Set - Create the equations set Navier-Stokes
EquationsNavierStokes = CMISS.Equations()
EquationsSetNavierStokes.EquationsCreateStart(EquationsNavierStokes)
# Set the equations matrices sparsity type
EquationsNavierStokes.sparsityType = CMISS.EquationsSparsityTypes.FULL
# Set the equations lumping type
EquationsNavierStokes.lumpingType = CMISS.EquationsLumpingTypes.UNLUMPED
# Set the equations set output
# (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
EquationsNavierStokes.outputType = CMISS.EquationsOutputTypes.NONE
# Finish the equations set equations
EquationsSetNavierStokes.EquationsCreateFinish()

# 2nd Equations Set - Create the equations set Characteristic
EquationsCharacteristic = CMISS.Equations()
EquationsSetCharacteristic.EquationsCreateStart(EquationsCharacteristic)
# Set the equations matrices sparsity type
EquationsCharacteristic.sparsityType = CMISS.EquationsSparsityTypes.SPARSE
# Set the equations set output
# (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
EquationsCharacteristic.outputType = CMISS.EquationsOutputTypes.NONE
# Finish the equations set equations
EquationsSetCharacteristic.EquationsCreateFinish()

#================================================================================================================================
#  Problems
#================================================================================================================================

# Start the creation of a problem.
Problem = CMISS.Problem()
Problem.CreateStart(ProblemUserNumber)
# Set the problem to be a dynamic Navier-Stokes problem
Problem.SpecificationSet(CMISS.ProblemClasses.FLUID_MECHANICS,CMISS.ProblemTypes.NAVIER_STOKES_EQUATION,ProblemSubtype)    
# Finish the creation of a problem.
Problem.CreateFinish()

#================================================================================================================================
#  Control Loop
#================================================================================================================================

# Start the creation of the problem control loop
ControlLoop = CMISS.ControlLoop()
Problem.ControlLoopCreateStart()
# Get the control loop
Problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE],ControlLoop)
# Set the times
ControlLoop.TimesSet(DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME,DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME,  
    DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT)
# Set the output timing
ControlLoop.TimeOutputSet(DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY)
# Finish creating the problem control loop
Problem.ControlLoopCreateFinish()

#================================================================================================================================
#  Solvers
#================================================================================================================================

# Start the creation of the problem solvers
DynamicSolverNavierStokes = CMISS.Solver()
NonlinearSolverNavierStokes = CMISS.Solver()
LinearSolverNavierStokes = CMISS.Solver()
NonlinearSolverCharacteristic = CMISS.Solver()
LinearSolverCharacteristic = CMISS.Solver()

Problem.SolversCreateStart()

# 1st Solver - Get the dynamic dynamic solver
Problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],SolverNavierStokesUserNumber,DynamicSolverNavierStokes)
# Set the output type
DynamicSolverNavierStokes.OutputTypeSet(DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE)
# Set theta
DynamicSolverNavierStokes.DynamicThetaSet(DYNAMIC_SOLVER_NAVIER_STOKES_THETA)
# Get the dynamic nonlinear solver
DynamicSolverNavierStokes.DynamicNonlinearSolverGet(NonlinearSolverNavierStokes)
# Set the nonlinear Jacobian type
NonlinearSolverNavierStokes.NewtonJacobianCalculationTypeSet(CMISS.JacobianCalculationTypes.EQUATIONS) #(.FD)
# Set the output type
NonlinearSolverNavierStokes.OutputTypeSet(NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE)
# Set the solver settings
NonlinearSolverNavierStokes.NewtonAbsoluteToleranceSet(ABSOLUTE_TOLERANCE)
NonlinearSolverNavierStokes.NewtonRelativeToleranceSet(RELATIVE_TOLERANCE)
# Get the dynamic nonlinear linear solver
NonlinearSolverNavierStokes.NewtonLinearSolverGet(LinearSolverNavierStokes)
# Set the output type
LinearSolverNavierStokes.OutputTypeSet(LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE)
# Set the solver settings
if (LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG):
    LinearSolverNavierStokes.LinearTypeSet(CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
    LinearSolverNavierStokes.LibraryTypeSet(CMISS_SOLVER_MUMPS_LIBRARY)
else:
    LinearSolverNavierStokes.LinearTypeSet(CMISS.LinearSolverTypes.ITERATIVE)
    LinearSolverNavierStokes.LinearIterativeMaximumIterationsSet(MAXIMUM_ITERATIONS)
    LinearSolverNavierStokes.LinearIterativeDivergenceToleranceSet(DIVERGENCE_TOLERANCE)
    LinearSolverNavierStokes.LinearIterativeRelativeToleranceSet(RELATIVE_TOLERANCE)
    LinearSolverNavierStokes.LinearIterativeAbsoluteToleranceSet(ABSOLUTE_TOLERANCE)
    LinearSolverNavierStokes.LinearIterativeGMRESRestartSet(RESTART_VALUE)

# 2nd Solver - Get the static nonlinear solver
Problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],SolverCharacteristicUserNumber,NonlinearSolverCharacteristic)
# Set the nonlinear Jacobian type
NonlinearSolverCharacteristic.NewtonJacobianCalculationTypeSet(CMISS.JacobianCalculationTypes.EQUATIONS) #(.FD)
# Set the output type
NonlinearSolverCharacteristic.OutputTypeSet(NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE)
# Set the solver settings
NonlinearSolverCharacteristic.NewtonAbsoluteToleranceSet(ABSOLUTE_TOLERANCE)
NonlinearSolverCharacteristic.NewtonRelativeToleranceSet(RELATIVE_TOLERANCE)
# Get the nonlinear linear solver
NonlinearSolverCharacteristic.NewtonLinearSolverGet(LinearSolverCharacteristic)
# Set the output type
LinearSolverCharacteristic.OutputTypeSet(LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE)
# Set the solver settings
if (LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG):
    LinearSolverCharacteristic.LinearTypeSet(CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
    LinearSolverCharacteristic.LibraryTypeSet(CMISS_SOLVER_MUMPS_LIBRARY)
else:
    LinearSolverCharacteristic.LinearTypeSet(CMISS.LinearSolverTypes.ITERATIVE)
    LinearSolverCharacteristic.LinearIterativeMaximumIterationsSet(MAXIMUM_ITERATIONS)
    LinearSolverCharacteristic.LinearIterativeDivergenceToleranceSet(DIVERGENCE_TOLERANCE)
    LinearSolverCharacteristic.LinearIterativeRelativeToleranceSet(RELATIVE_TOLERANCE)
    LinearSolverCharacteristic.LinearIterativeAbsoluteToleranceSet(ABSOLUTE_TOLERANCE)
    LinearSolverCharacteristic.LinearIterativeGMRESRestartSet(RESTART_VALUE)

# Finish the creation of the problem solver
Problem.SolversCreateFinish()

# 3rd Solver - CellML Solver
if (cellmlFlag):

    # Create the problem solver CellML equations
    CellMLSolver = CMISS.Solver()
    CellMLEquations = CMISS.CellMLEquations()

    Problem.CellMLEquationsCreateStart()

    DynamicSolverNavierStokes.NewtonCellMLSolverGet(CellMLSolver) 
    CellMLSolver.CellMLEquationsGet(CellMLEquations)
    CellMLEquations.CellMLAdd(CellML)

    Problem.CellMLEquationsCreateFinish()

#================================================================================================================================
#  Solver Equations
#================================================================================================================================

# Start the creation of the problem solver equations
DynamicSolverNavierStokes = CMISS.Solver()
NonlinearSolverCharacteristic = CMISS.Solver()
SolverEquationsNavierStokes = CMISS.SolverEquations()
SolverEquationsCharacteristic = CMISS.SolverEquations()

Problem.SolverEquationsCreateStart()

# 1st Solver - Get the dynamic Navier-Stokes solver equations
Problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],SolverNavierStokesUserNumber,DynamicSolverNavierStokes)
DynamicSolverNavierStokes.SolverEquationsGet(SolverEquationsNavierStokes)
# Set the solver equations sparsity
SolverEquationsNavierStokes.sparsityType = CMISS.SolverEquationsSparsityTypes.SPARSE
# Add in the equations set
EquationsSetNavierStokes = SolverEquationsNavierStokes.EquationsSetAdd(EquationsSetNavierStokes)

# 2nd Solver - Get the static nonlinear solver equations
Problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],SolverCharacteristicUserNumber,NonlinearSolverCharacteristic)
NonlinearSolverCharacteristic.SolverEquationsGet(SolverEquationsCharacteristic)
# Set the solver equations sparsity
SolverEquationsCharacteristic.sparsityType = CMISS.SolverEquationsSparsityTypes.SPARSE
# Add in the equations set
EquationsSetCharacteristic = SolverEquationsCharacteristic.EquationsSetAdd(EquationsSetCharacteristic)

# Finish the creation of the problem solver equations
Problem.SolverEquationsCreateFinish()

#================================================================================================================================
#  Boundary Conditions
#================================================================================================================================

# Navier-Stokes BCs
BoundaryConditionsNavierStokes = CMISS.BoundaryConditions()
SolverEquationsNavierStokes.BoundaryConditionsCreateStart(BoundaryConditionsNavierStokes)

# Flow(inlet)
BoundaryConditionsNavierStokes.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
    1,1,1,1,CMISS.BoundaryConditionsTypes.FIXED_INLET,Q[0])
    #(versionNumber,derivativeNumber,nodeUserNumber,componentNumber,condition,value)

# Area(outlet)
for i in range (1,NumberOfTerminalNodes+1):
    BoundaryConditionsNavierStokes.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
        1,1,coupledNodeNumber[i],2,CMISS.BoundaryConditionsTypes.FIXED_OUTLET,Ae[i])

# Characteristic BCs
BoundaryConditionsCharacteristic = CMISS.BoundaryConditions()
SolverEquationsCharacteristic.BoundaryConditionsCreateStart(BoundaryConditionsCharacteristic)

# Area(outlet)
if (cellmlFlag):
    for i in range (1,NumberOfTerminalNodes+1):
        BoundaryConditionsCharacteristic.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
            1,1,coupledNodeNumber[i],2,CMISS.BoundaryConditionsTypes.FIXED_OUTLET,Ae[i])

# Finish the creation of the equations set boundary conditions
SolverEquationsNavierStokes.BoundaryConditionsCreateFinish()
SolverEquationsCharacteristic.BoundaryConditionsCreateFinish()

#================================================================================================================================
#  CellML Parameters
#================================================================================================================================

if (cellmlFlag):

    # Set CellML model parameters (Resistance/Capacitance) at boundary nodes
    resistanceComponent = CellML.FieldComponentGet(ResistanceModelIndex,CMISS.CellMLFieldTypes.PARAMETERS,"equations/Resistance")
    # Branches
    for i in range (1,NumberOfTerminalNodes+1):
        BoundaryNodeDomain = Decomposition.NodeDomainGet(coupledNodeNumber[i],1)
        if (BoundaryNodeDomain == ComputationalNodeNumber):
            CellMLParametersField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,
                coupledNodeNumber[i],resistanceComponent,resistanceProximal[i])
  
#================================================================================================================================
#  Run Solvers
#================================================================================================================================

# Solve the problem
print "Solving problem..."
Problem.Solve()
print "Problem solved!"
print "#"

#================================================================================================================================
#  Finish Program
#================================================================================================================================
