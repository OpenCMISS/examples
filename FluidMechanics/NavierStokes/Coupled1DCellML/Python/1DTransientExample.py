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
#  Start Program
#================================================================================================================================

# Program Variables
EquationsSetFieldUserNumberNavierStokes = 1337
EquationsSetFieldUserNumberCharacteristic = 1338

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

CellMLUserNumber = 13
CellMLModelsFieldUserNumber = 14
CellMLStateFieldUserNumber = 15
CellMLIntermediateFieldUserNumber = 16
CellMLParametersFieldUserNumber = 17
MaterialsFieldUserNumberCellML = 18
  
DomainUserNumber = 1
SolverCharacteristicUserNumber = 1
SolverNavierStokesUserNumber = 2
MaterialsFieldUserNumberMu = 1
MaterialsFieldUserNumberRho = 2
MaterialsFieldUserNumberK = 3
MaterialsFieldUserNumberBs = 4
MaterialsFieldUserNumberAs = 5
MaterialsFieldUserNumberRe = 6
MaterialsFieldUserNumberFr = 7
MaterialsFieldUserNumberSt = 8
MaterialsFieldUserNumberA0 = 9
MaterialsFieldUserNumberBeta = 10
MaterialsFieldUserNumberE = 11
MaterialsFieldUserNumberH0 = 12

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
CMISS.ErrorHandlingModeSet(CMISS.ErrorHandlingModes.TRAP_ERROR)
CMISS.OutputSetOn("Testing")

# Get computational nodes information
NumberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
ComputationalNodeNumber = CMISS.ComputationalNodeNumberGet()

#================================================================================================================================
#  Problem Control Panel
#================================================================================================================================

# Set geometry parameters
NumberOfDimensions = 1
TotalNumberOfElements = 6
NumberOfNodesSpace = 13
NumberOfNodesFlow = NumberOfNodesSpace
NumberOfNodesArea = NumberOfNodesSpace
TotalNumberOfNodes = NumberOfNodesSpace*3

# Set basis parameters
BasisUserNumberSpace = 1
BasisUserNumberFlow = 2
BasisUserNumberArea = 3
BasisXiInterpolationSpace = 2
BasisXiInterpolationFlow = 2
BasisXiInterpolationArea = 2

# Set interpolation parameters
BasisXiGaussSpace = 3
BasisXiGaussFlow = 3
BasisXiGaussArea = 3

# Set lumped model parameters
resistanceProximal = 1.7025E+7
resistanceDistal = 0.0
capacitance = 0.0
cellmlFlag = True
windkesselFlag = False

# Set material parameters
MU_PARAM_NAVIER_STOKES = 0.0035           # Mu (Pa.s)
RHO_PARAM_NAVIER_STOKES = 1050.0          # Rho (kg/m3)
E_PARAM_NAVIER_STOKES = 0.8E+6            # Elasticity (Pa)
H0_PARAM_NAVIER_STOKES = 0.5E-3           # Wall Thickness (m)
A0_PARAM = [0]*(TotalNumberOfElements+1)  # Area (m2)
Beta = [0]*(TotalNumberOfElements+1)      # Beta (Pa/m)
A0_PARAM[1] = 19.6e-6                 
A0_PARAM[2] = 19.6e-6
A0_PARAM[3] = 12.8e-6
A0_PARAM[4] = 12.8e-6
A0_PARAM[5] = 12.8e-6
A0_PARAM[6] = 12.8e-6
for i in range(1,TotalNumberOfElements+1):
    Beta[i] = (4.0*(3.1416**0.5)*E_PARAM_NAVIER_STOKES*H0_PARAM_NAVIER_STOKES)/(3.0*A0_PARAM[i])
    
# Set reference values
Qs = 10.0e-6                              # Flow (m3/s)
As = 19.6e-6                              # Area (m2)
Xs = 0.1                                  # Distance (m)
Ts = 0.1                                  # Time (s)
K = 4.0/3.0                               # Parabolic flow section
Bs = (4.0*(3.1416**0.5)*E_PARAM_NAVIER_STOKES*H0_PARAM_NAVIER_STOKES)/(3.0*As)  # Beta
Re = 8.0*3.1416*(MU_PARAM_NAVIER_STOKES*Xs)/(Qs*RHO_PARAM_NAVIER_STOKES)        # Reynolds number
Fr = ((As**2.5)/(Qs**2))*(Bs/(2.0*RHO_PARAM_NAVIER_STOKES))                     # Froude number
St = (As*Xs)/(Ts*Qs)                                                            # Strouhal number

# Set initial conditions
Q1 = 7.0
Q2 = 3.5
Q3 = 3.5
A1 = 1.0
A2 = 0.653
A3 = 0.653

# Set terminal nodes 
# (terminal node with either versions set so that method of characteristics may
# be used for 1D-0D coupling or just a regular terminal boundary condition)
coupledNodeNumber1 = 11
coupledNodeNumber2 = 13

# Set output parameters
# (NONE/PROGRESS/TIMING/SOLVER/MATRIX)
DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE = CMISS.SolverOutputTypes.MATRIX
NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE = CMISS.SolverOutputTypes.MATRIX
LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE = CMISS.SolverOutputTypes.MATRIX
# (NONE/TIMING/MATRIX/ELEMENT)
EQUATIONS_NAVIER_STOKES_OUTPUT = CMISS.SolverOutputTypes.MATRIX

# Set time parameters
DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME = 0.0
DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME = 0.01
DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT = 0.01
DYNAMIC_SOLVER_NAVIER_STOKES_THETA = [1.0/2.0]

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
#  Reading Mesh
#================================================================================================================================

# User numbers
NumberOfXi = 3
MeshNumberOfComponents = 1
MeshTotalNumberOfElements = 1
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
REGION.ReadMesh(BasisUserNumberSpace,MeshUserNumber,GeometricFieldUserNumber,
    "CMISS","Input/Geometry","Input/Geometry")
    
# Extracting geometry data
FieldVariable = 1
ElementNodes = [0,0,0]*(TotalNumberOfElements)
xValues = [0,0]*(NumberOfNodesSpace-1)
yValues = [0,0]*(NumberOfNodesSpace-1)
zValues = [0,0]*(NumberOfNodesSpace-1)
for i in range(1,TotalNumberOfElements+1):
    ElementNodes[i-1] = REGION.MESHES.MeshElementsNodesGet(MeshUserNumber,MeshNumberOfComponents,i)
for i in range(1,NumberOfNodesSpace+1):
    xValues[i-1] = REGION.FIELDS.FieldParameterSetNodeValuesGet(GeometricFieldUserNumber,FieldVariable,i,1)
    yValues[i-1] = REGION.FIELDS.FieldParameterSetNodeValuesGet(GeometricFieldUserNumber,FieldVariable,i,2)
    zValues[i-1] = REGION.FIELDS.FieldParameterSetNodeValuesGet(GeometricFieldUserNumber,FieldVariable,i,3)

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
MeshNumberOfComponents = 1
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
MeshComponentNumberSpace = 1
MeshComponentNumberFlow = 1
MeshComponentNumberArea = 1
# Specify space mesh component
MeshElementsSpace.CreateStart(Mesh,MeshComponentNumberSpace,BasisSpace)
for i in range(1,TotalNumberOfElements+1):
    MeshElementsSpace.NodesSet(i,ElementNodes[i-1])
    
    #   NODES
    #               7-10-11
    #              /
    #             6
    #            /
    #   1-2-3-4-5
    #            \   
    #             8 
    #              \
    #               9-12-13
    #  
    #   ELEMENTS
    #               --5--
    #              /
    #             3
    #            /
    #   --1---2--
    #            \   
    #             4 
    #              \
    #               --6--
    
# Set versions at bifurcation
# (globalElementNumber,versionNumber,derivativeNumber,localElementNodeNumber)
MeshElementsSpace.LocalElementNodeVersionSet(2,1,1,3) 
MeshElementsSpace.LocalElementNodeVersionSet(3,2,1,1) 
MeshElementsSpace.LocalElementNodeVersionSet(4,3,1,1) 
    
MeshElementsSpace.CreateFinish()
# Specify flow mesh component
if (BasisXiInterpolationFlow == BasisXiInterpolationSpace):
    MeshElementsFlow = MeshElementsSpace
# Specify area mesh component
if (BasisXiInterpolationArea == BasisXiInterpolationSpace):
    MeshElementsArea = MeshElementsSpace
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
for ComponentNumber in range(1,2):
    GeometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,ComponentNumber,MeshComponentNumberSpace)
# Finish creating the field
GeometricField.CreateFinish()

# Initialise the equations set geometric field variables
for UserNodeNumber in range(1,NumberOfNodesSpace+1):
    VersionNumber = 1
    # (field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value)
    GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        VersionNumber,1,UserNodeNumber,1,xValues[UserNodeNumber-1][VersionNumber-1])
    GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        VersionNumber,1,UserNodeNumber,2,yValues[UserNodeNumber-1][VersionNumber-1])
    GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        VersionNumber,1,UserNodeNumber,3,zValues[UserNodeNumber-1][VersionNumber-1])
        
UserNodeNumber = 5
for VersionNumber in range(2,4):
    # (field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value)
    GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        VersionNumber,1,UserNodeNumber,1,xValues[UserNodeNumber-1][VersionNumber-1])
    GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        VersionNumber,1,UserNodeNumber,2,yValues[UserNodeNumber-1][VersionNumber-1])
    GeometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        VersionNumber,1,UserNodeNumber,3,zValues[UserNodeNumber-1][VersionNumber-1])

#================================================================================================================================
#  Equations Sets
#================================================================================================================================

# Create the equations set for Navier-Stokes
EquationsSetNavierStokes = CMISS.EquationsSet()
EquationsSetFieldNavierStokes = CMISS.Field()
# Set the equations set to be a dynamic Navier-Stokes problem
EquationsSetNavierStokes.CreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField,
        CMISS.EquationsSetClasses.FLUID_MECHANICS,CMISS.EquationsSetTypes.NAVIER_STOKES_EQUATION,
        EquationsSetSubtype,EquationsSetFieldUserNumberNavierStokes,
        EquationsSetFieldNavierStokes)
# Finish creating the equations set
EquationsSetNavierStokes.CreateFinish()

# Create the equations set for Characteristic
EquationsSetCharacteristic = CMISS.EquationsSet()
EquationsSetFieldCharacteristic = CMISS.Field()
# Set the equations set to be a static Nonlinear problem
EquationsSetCharacteristic.CreateStart(EquationsSetUserNumberCharacteristic,Region,GeometricField,
        CMISS.EquationsSetClasses.FLUID_MECHANICS,CMISS.EquationsSetTypes.CHARACTERISTIC_EQUATION,
        EquationsSetCharacteristicSubtype,EquationsSetFieldUserNumberCharacteristic,
        EquationsSetFieldCharacteristic)
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
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,1,1,Q1)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,2,1,Q1)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,3,1,Q1)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,4,1,Q1)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,5,1,Q1)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,2,1,5,1,Q2)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,3,1,5,1,Q3)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,6,1,Q2)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,7,1,Q2)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,8,1,Q3)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,9,1,Q3)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,10,1,Q2)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,11,1,Q2)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,12,1,Q3)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,13,1,Q3)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,1,2,A1)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,2,2,A1)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,3,2,A1)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,4,2,A1)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,5,2,A1)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,2,1,5,2,A2)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,3,1,5,2,A3)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,6,2,A2)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,7,2,A2)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,8,2,A3)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,9,2,A3)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,10,2,A2)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,11,2,A2)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,12,2,A3)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,
    CMISS.FieldParameterSetTypes.VALUES,1,1,13,2,A3)

DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
    CMISS.FieldParameterSetTypes.VALUES,1,1,5,1,0.0)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
    CMISS.FieldParameterSetTypes.VALUES,2,1,5,1,0.0)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
    CMISS.FieldParameterSetTypes.VALUES,3,1,5,1,0.0)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
    CMISS.FieldParameterSetTypes.VALUES,1,1,5,2,0.0)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
    CMISS.FieldParameterSetTypes.VALUES,2,1,5,2,0.0)
DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
    CMISS.FieldParameterSetTypes.VALUES,3,1,5,2,0.0)

#================================================================================================================================
#  Materials Field
#================================================================================================================================

# Create the equations set materials field variables for Navier-Stokes
MaterialsFieldNavierStokes = CMISS.Field()
EquationsSetNavierStokes.MaterialsCreateStart(MaterialsFieldUserNumber,MaterialsFieldNavierStokes)
# Set the field label
MaterialsFieldNavierStokes.VariableLabelSet(CMISS.FieldVariableTypes.U,'Materials')
# Finish the equations set materials field variables
EquationsSetNavierStokes.MaterialsCreateFinish()

# Create the equations set materials field variables for Characteristic
EquationsSetCharacteristic.MaterialsCreateStart(MaterialsFieldUserNumber,MaterialsFieldNavierStokes)
# Finish the equations set materials field variables
EquationsSetCharacteristic.MaterialsCreateFinish()

# Initialise the equations set materials field variables
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberMu,MU_PARAM_NAVIER_STOKES)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberRho,RHO_PARAM_NAVIER_STOKES)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
     MaterialsFieldUserNumberE,E_PARAM_NAVIER_STOKES)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
     MaterialsFieldUserNumberH0,H0_PARAM_NAVIER_STOKES)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberK,K)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberBs,Bs)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberAs,As)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberRe,Re)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
     MaterialsFieldUserNumberFr,Fr)
MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
     MaterialsFieldUserNumberSt,St)
for i in range(1,TotalNumberOfElements+1):
    MaterialsFieldNavierStokes.ParameterSetUpdateElementDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        i,MaterialsFieldUserNumberA0,A0_PARAM[i])
    MaterialsFieldNavierStokes.ParameterSetUpdateElementDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        i,MaterialsFieldUserNumberBeta,Beta[i])

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
    # Set values at coupled node 1
    # (versionNumber,derivativeNumber,userNodeNumber,componentNumber,value)
    MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
        CMISS.FieldParameterSetTypes.VALUES,1,1,coupledNodeNumber1,pVesselWallComponent,pVesselWall)
    MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
        CMISS.FieldParameterSetTypes.VALUES,1,1,coupledNodeNumber1,pExternalComponent,pExternal)
    # Set values at coupled node 2
    # (versionNumber,derivativeNumber,userNodeNumber,componentNumber,value)
    MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
        CMISS.FieldParameterSetTypes.VALUES,1,1,coupledNodeNumber2,pVesselWallComponent,pVesselWall)
    MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.V,
        CMISS.FieldParameterSetTypes.VALUES,1,1,coupledNodeNumber2,pExternalComponent,pExternal)

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

# Normal Wave Direction for bifurcation
nodeIdx = 5
# 1 inlet/parent
componentIdx = 1 # Incoming
versionIdx = 1
VALUE = 1.0
# (versionNumber,derivativeNumber,userNodeNumber,componentNumber,value)
IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
    versionIdx,1,nodeIdx,componentIdx,VALUE)

# 2 outlet/daughters
componentIdx = 2 # Outgoing
VALUE = -1.0
for versionIdx in range(2,4):
   IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,  
       versionIdx,1,nodeIdx,componentIdx,VALUE)

# Normal Wave Direction for coupling
if (cellmlFlag):

    # Incoming normals
    nodeIdx = 11
    componentIdx = 1
    VALUE = 1.0
    versionIdx = 1
    IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        versionIdx,1,nodeIdx,componentIdx,VALUE)
    nodeIdx = 13
    IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        versionIdx,1,nodeIdx,componentIdx,VALUE)

    # Outgoing normals
    nodeIdx = 11
    componentIdx = 2
    VALUE = -1.0
    versionIdx = 1
    IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        versionIdx,1,nodeIdx,componentIdx,VALUE)
    nodeIdx = 13
    IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
        versionIdx,1,nodeIdx,componentIdx,VALUE)

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
        # (field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value)
        DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U1,CMISS.FieldParameterSetTypes.VALUES,
            1,1,coupledNodeNumber1,pCellMLComponent,pCellML)
        DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U1,CMISS.FieldParameterSetTypes.VALUES,
            1,1,coupledNodeNumber2,pCellMLComponent,pCellML)
        DependentFieldNavierStokes.ParameterSetCreate(CMISS.FieldVariableTypes.U1,CMISS.FieldParameterSetTypes.PREVIOUS_VALUES)
        DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U1,
            CMISS.FieldParameterSetTypes.PREVIOUS_VALUES,1,1,coupledNodeNumber1,pCellMLComponent,pPrevious)
        DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U1,
            CMISS.FieldParameterSetTypes.PREVIOUS_VALUES,1,1,coupledNodeNumber2,pCellMLComponent,pPrevious)

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
        # (field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value)
        DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U1,CMISS.FieldParameterSetTypes.VALUES,
            1,1,coupledNodeNumber1,pCellMLComponent,pCellML)
        DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U1,CMISS.FieldParameterSetTypes.VALUES,
            1,1,coupledNodeNumber2,pCellMLComponent,pCellML)

        DependentFieldNavierStokes.ParameterSetCreate(CMISS.FieldVariableTypes.U1,CMISS.FieldParameterSetTypes.PREVIOUS_VALUES)

        DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U1, 
            CMISS.FieldParameterSetTypes.PREVIOUS_VALUES,1,1,coupledNodeNumber1,pCellMLComponent,pPrevious)
        DependentFieldNavierStokes.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U1, 
            CMISS.FieldParameterSetTypes.PREVIOUS_VALUES,1,1,coupledNodeNumber2,pCellMLComponent,pPrevious)

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

# CellML Solver
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

# Navier-Stokes equations BCs
BoundaryConditionsNavierStokes = CMISS.BoundaryConditions()
SolverEquationsNavierStokes.BoundaryConditionsCreateStart(BoundaryConditionsNavierStokes)
# FLOW boundary conditions
CONDITION = CMISS.BoundaryConditionsTypes.FIXED_INLET
VALUE = Q1
BoundaryConditionsNavierStokes.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
    1,1,1,1,CONDITION,VALUE)
# AREA boundary conditions
CONDITION = CMISS.BoundaryConditionsTypes.FIXED_OUTLET
versionIdx = 1
VALUE = A2
BoundaryConditionsNavierStokes.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
    versionIdx,1,coupledNodeNumber1,2,CONDITION,VALUE)
BoundaryConditionsNavierStokes.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
    versionIdx,1,coupledNodeNumber2,2,CONDITION,VALUE)

# Characteristic equations BCs
BoundaryConditionsCharacteristic = CMISS.BoundaryConditions()
SolverEquationsCharacteristic.BoundaryConditionsCreateStart(BoundaryConditionsCharacteristic)
# AREA boundary conditions
if (cellmlFlag):
    CONDITION = CMISS.BoundaryConditionsTypes.FIXED_OUTLET
    versionIdx = 1
    VALUE = A2
    BoundaryConditionsCharacteristic.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
        versionIdx,1,coupledNodeNumber1,2,CONDITION,VALUE)
    VALUE = A3
    BoundaryConditionsCharacteristic.SetNode(DependentFieldNavierStokes,CMISS.FieldVariableTypes.U,
        versionIdx,1,coupledNodeNumber2,2,CONDITION,VALUE)

# (versionNumber,derivativeNumber,nodeUserNumber,componentNumber,condition,value)
# Finish the creation of the equations set boundary conditions
SolverEquationsNavierStokes.BoundaryConditionsCreateFinish()
SolverEquationsCharacteristic.BoundaryConditionsCreateFinish()

#================================================================================================================================
#  CellML Parameters
#================================================================================================================================

if (cellmlFlag):
    # Set CellML model parameters (Resistance/Capacitance) at boundary nodes
    resistanceComponent = CellML.FieldComponentGet(ResistanceModelIndex,CMISS.CellMLFieldTypes.PARAMETERS,"equations/Resistance")
    
    # Branch 1
    BoundaryNodeDomain = Decomposition.NodeDomainGet(coupledNodeNumber1,1)
    if (BoundaryNodeDomain == ComputationalNodeNumber):
        CellMLParametersField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,
            coupledNodeNumber1,resistanceComponent,resistanceProximal)

    # Branch 2 
    BoundaryNodeDomain = Decomposition.NodeDomainGet(coupledNodeNumber2,1)
    if (BoundaryNodeDomain == ComputationalNodeNumber):
        resistanceProximal = resistanceProximal*0.5
        CellMLParametersField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,
            coupledNodeNumber2,resistanceComponent,resistanceProximal)
  
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
