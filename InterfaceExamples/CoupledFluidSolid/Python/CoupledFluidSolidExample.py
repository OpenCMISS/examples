#> \file
#> \author Soroush Safaei
#> \brief This is an example program which solves a weakly coupled 
#>  FiniteElasticity-ALENavierStokes equation using OpenCMISS calls.
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
#> Contributor(s): Andreas Hessenthaler
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
#> OpenCMISS/examples/InterfaceExamples/CoupledFluidSolid/Python/CoupledFluidSolidExample.py
#>

#================================================================================================================================
#  Start Program
#================================================================================================================================

Plate2D = 2
Plate3D = 3

SolidCoordinateSystemUserNumber     = 1
FluidCoordinateSystemUserNumber     = 2
InterfaceCoordinateSystemUserNumber = 4
  
SolidRegionUserNumber = 6
FluidRegionUserNumber = 7
InterfaceUserNumber   = 40
  
SolidMeshUserNumber     = 19
FluidMeshUserNumber     = 20
InterfaceMeshUserNumber = 22
MovingMeshUserNumber    = 107
  
SolidDecompositionUserNumber     = 24
FluidDecompositionUserNumber     = 25
InterfaceDecompositionUserNumber = 27
  
SolidGeometricFieldUserNumber     = 29
FluidGeometricFieldUserNumber     = 30
InterfaceGeometricFieldUserNumber = 32
  
SolidEquationsSetUserNumber  = 34
FluidEquationsSetUserNumber  = 35
InterfaceConditionUserNumber = 42
  
SolidDependentFieldUserNumber = 37
FluidDependentFieldUserNumber = 38
LagrangeFieldUserNumber  = 44

CoupledProblemUserNumber = 46

SolidEquationsSetFieldUserNumber = 49
FluidEquationsSetFieldUserNumber = 50
  
BasisSpaceSolidUserNumber   = 100
BasisDisplacementUserNumber = 101
BasisHydrostaticPressureUserNumber = 102
BasisSpaceFluidUserNumber = 103
BasisVelocityUserNumber   = 104
BasisPressureUserNumber   = 105
InterfaceBasisUserNumber  = 12
InterfaceMappingBasisUserNumber = 47
  
FibreFieldUserNumber = 106
MovingMeshFieldUserNumber = 108
MovingMeshEquationsSetUserNumber = 109
FluidMaterialFieldUserNumber = 117
SolidMaterialFieldUserNumber = 110
EquationsSetFieldMovingMeshUserNumber = 111
SourceFieldUserNumber = 112
  
DynamicSolverIndex = 1
LinearSolverMovingMeshIndex = 2
  
FluidMaterialFieldComponentMu  = 1
FluidMaterialFieldComponentRho = 2
MaterialFieldMovingMeshUserNumberK    = 1
IndependentFieldMovingMeshUserNumberK = 1
DependentFieldMovingMeshUserNumber = 118
MaterialFieldMovingMeshUserNumber  = 119
IndependentField2UserNumber = 120
IndependentFieldMovingMeshUserNumber      = 121
LinearSolverMovingMeshEquationsUserNumber = 122

SolidEquationsSetIndex  = 1
FluidEquationsSetIndex  = 2
InterfaceConditionIndex = 1
SolidMeshIndex = 1
FluidMeshIndex = 2

derivIdx   = 1
versionIdx = 1

FileReadDiagnostics            = True
ExampleFileProgressDiagnostics = True
GeometryCheck                  = False
GravityFlag                    = False
CheckWithoutInterfaceCondition = False
SetupOutput                    = True

# Quadratic interpolation for fluid space and velocities and solid displacement
InterpolationTypeSpace        = 2
InterpolationTypeVelocity     = 2
InterpolationTypeDisplacement = 2
# Linear interpolation for pressures
InterpolationTypePressure     = 1
# Interpolation type on interface matches interpolation type for fluid velocities, solid displacements
InterpolationTypeInterface = InterpolationTypeDisplacement
    
#================================================================================================================================
#  Initialise OpenCMISS
#================================================================================================================================

# Import the libraries (OpenCMISS,python,numpy,scipy)
import numpy,csv,time,sys,os,pdb
sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))
from opencmiss import iron

# Diagnostics
#iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",[""])
#iron.ErrorHandlingModeSet(iron.ErrorHandlingModes.TRAP_ERROR)
#iron.OutputSetOn("Testing")

# Get the computational nodes info
NumberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
ComputationalNodeNumber    = iron.ComputationalNodeNumberGet()
        
#================================================================================================================================
#  Initial Data & Default Values
#================================================================================================================================

# Set initial values
InitialFieldNavierStokes = []
InitialFieldMovingMesh   = []
InitialFieldNavierStokes.append(0.0)
InitialFieldNavierStokes.append(0.0)
InitialFieldNavierStokes.append(0.0)
InitialFieldMovingMesh.append(0.0)
InitialFieldMovingMesh.append(0.0)
InitialFieldMovingMesh.append(0.0)
# (NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
LinearSolverMovingMesh_OutputType=iron.SolverOutputTypes.NONE
DynamicSolver_OutputType=iron.SolverOutputTypes.NONE
LinearSolver_OutputType=iron.SolverOutputTypes.NONE
NonlinearSolver_OutputType=iron.SolverOutputTypes.NONE
# (NoOutput/TimingOutput/MatrixOutput/ElementOutput)
EquationsNavierStokesOutput=iron.SolverOutputTypes.NONE
# Set result output parameter
OutputFrequency = 1

# Choose 2D or 3D case
MaterialSpecification = Plate2D
# Set solver parameters
RelativeTolerance      = 1.0E-4    #default: 1.0E-05
AbsoluteTolerance      = 1.0E-4    #default: 1.0E-10
DivergenceTolerance    = 1.0E5     #default: 1.0E5
MaximumIterations      = 100000000 #default: 100000
MaxFunctionEvaluations = 100000
RestartValue           = 30        #default: 30
LinesearchAlpha        = 1.0
MovingMeshParameterK   = 1.0       #default
DynamicSolver_Theta    = [1.0]
StartTime    = 0.0
StopTime     = 20.0
TimeStepSize = 1.0
if (MaterialSpecification == Plate2D):
    StopTime     = 100.0
    TimeStepSize = 0.05

# Material properties
# NOTE: USE OF SI UNITS unless comment
# Low density fluid, rubber-like solid
FluidDynamicViscosity = 0.05  # kg / (m s)
FluidDensity  = 100           # kg m^-3
SolidDensity  = 300           # kg m^-3
YoungsModulus = 2.3E4         # Pa
PoissonsRatio = 0.49          # [.]
# Neo-Hookean material law
ShearModulus  = YoungsModulus/(2.0*(1.0 + PoissonsRatio))     # N / m^2
BulkModulus   = YoungsModulus/(3.0*(1.0-2.0 * PoissonsRatio))
MooneyRivlin1 = 0.5 * ShearModulus                            # N / m^2
MooneyRivlin2 = 0.0

# Set geometric dimension n gravity
if (MaterialSpecification == Plate2D):
    NumberOfDimensions = 2
    Gravity = numpy.zeros(NumberOfDimensions)
    Gravity = [0.0,9.81] #in m s^-2
elif (MaterialSpecification == Plate3D):
    NumberOfDimensions = 3
    CheckWithoutInterfaceCondition = True #if set to true we remove all Lagrange field dofs by setting them as zero dirichlet BC
    if (CheckWithoutInterfaceCondition):
        GravityFlag = True
        Gravity = numpy.zeros(NumberOfDimensions)
        Gravity = [0.0,9.81,0.0] #in m s^-2
        
#================================================================================================================================
#  Reading geometry from files
#================================================================================================================================

print "Reading geometry from files.."
# Read solid node numbers
if (MaterialSpecification == Plate2D):
    f = './Input/2D/Solid/Plate2DSolidNodes.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/Solid/Plate3DSolidNodes.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            NumberOfSolidNodes = int(row[1])
            SolidNodeNumbers = NumberOfSolidNodes*[0]
        else:
            SolidNodeNumbers[rownum-1] = int(row[0])
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Solid nodes: ", SolidNodeNumbers

#------------------------------------------------------------------

# Read fluid node numbers
if (MaterialSpecification == Plate2D):
    f = './Input/2D/Fluid/Plate2DFluidNodes.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/Fluid/Plate3DFluidNodes.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            NumberOfFluidNodes = int(row[1])
            FluidNodeNumbers = NumberOfFluidNodes*[0]
        else:
            FluidNodeNumbers[rownum-1] = int(row[0])
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Fluid nodes: ", FluidNodeNumbers

#------------------------------------------------------------------

# Read interface node numbers
if (MaterialSpecification == Plate2D):
    f = './Input/2D/Interface/Plate2DInterfaceNodes.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/Interface/Plate3DInterfaceNodes.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            NumberOfInterfaceNodes = int(row[1])
            InterfaceNodeNumbersForGeometry = NumberOfInterfaceNodes*[0]
        else:
            InterfaceNodeNumbersForGeometry[rownum-1] = int(row[0])
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Interface nodes: ", InterfaceNodeNumbersForGeometry

#------------------------------------------------------------------

# Read geometry
if (NumberOfDimensions == 3):
    f = open('./Input/3D/Solid/Plate3DSolidX.csv')
    line = f.read().splitlines()
    SolidGeometryX = len(line)*[0]
    for s in range(0,len(line)):
        SolidGeometryX[s]=float(line[s])
    f.close()
    if (FileReadDiagnostics): print "Solid x coordinates: ", SolidGeometryX

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/Solid/Plate2DSolidGeometry.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/Solid/Plate3DSolidGeometry.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            SolidGeometryY = NumberOfSolidNodes*[0]
            SolidGeometryZ = NumberOfSolidNodes*[0]
        else:
            SolidGeometryY[rownum-1] = float(row[0])
            SolidGeometryZ[rownum-1] = float(row[1])
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Solid y coordinates: ", SolidGeometryY
if (FileReadDiagnostics): print "Solid z coordinates: ", SolidGeometryZ

#------------------------------------------------------------------

if (NumberOfDimensions == 3):
    f = open('./Input/3D/Fluid/Plate3DFluidX.csv')
    line = f.read().splitlines()
    FluidGeometryX = len(line)*[0]
    for s in range(0,len(line)):
        FluidGeometryX[s]=float(line[s])
    f.close()
    if (FileReadDiagnostics): print "Fluid x coordinates: ", FluidGeometryX

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/Fluid/Plate2DFluidGeometry.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/Fluid/Plate3DFluidGeometry.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            FluidGeometryY = NumberOfFluidNodes*[0]
            FluidGeometryZ = NumberOfFluidNodes*[0]
        else:
            FluidGeometryY[rownum-1] = float(row[0])
            FluidGeometryZ[rownum-1] = float(row[1])
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Fluid y coordinates: ", FluidGeometryY
if (FileReadDiagnostics): print "Fluid z coordinates: ", FluidGeometryZ

#------------------------------------------------------------------

if (NumberOfDimensions == 3):
    f = open('./Input/3D/Interface/Plate3DInterfaceX.csv')
    line = f.read().splitlines()
    InterfaceGeometryX = len(line)*[0]
    for s in range(0,len(line)):
        InterfaceGeometryX[s]=float(line[s])
    f.close()
    if (FileReadDiagnostics): print "Interface x coordinates: ", InterfaceGeometryX

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/Interface/Plate2DInterfaceGeometry.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/Interface/Plate3DInterfaceGeometry.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            InterfaceGeometryY = NumberOfInterfaceNodes*[0]
            InterfaceGeometryZ = NumberOfInterfaceNodes*[0]
        else:
            InterfaceGeometryY[rownum-1] = float(row[0])
            InterfaceGeometryZ[rownum-1] = float(row[1])
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Interface y coordinates: ", InterfaceGeometryY
if (FileReadDiagnostics): print "Interface z coordinates: ", InterfaceGeometryZ

#------------------------------------------------------------------

# Read element information
if (MaterialSpecification == Plate2D):
    f = './Input/2D/Solid/Plate2DSolidElementNodeNumbers.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/Solid/Plate3DSolidElementNodeNumbers.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            NumberOfSolidElements = int(row[0])
            SolidElementNodes = NumberOfSolidElements*[(3**NumberOfDimensions)*[0]]
        else:
            row = [int(i) for i in row]
            SolidElementNodes[rownum-1] = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Solid element node numbers: ", SolidElementNodes

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/Fluid/Plate2DFluidElementNodeNumbers.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/Fluid/Plate3DFluidElementNodeNumbers.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            NumberOfFluidElements = int(row[0])
            FluidElementNodes = NumberOfFluidElements*[(3**NumberOfDimensions)*[0]]
        else:
            row = [int(i) for i in row]
            FluidElementNodes[rownum-1] = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Fluid element node numbers: ", FluidElementNodes

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/Interface/Plate2DInterfaceElementNodeNumbers.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/Interface/Plate3DInterfaceElementNodeNumbers.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            NumberOfInterfaceElements = int(row[0])
            InterfaceElementNodes = NumberOfInterfaceElements*[(3**(NumberOfDimensions-1))*[0]]
        else:
            row = [int(i) for i in row]
            InterfaceElementNodes[rownum-1] = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Interface element node numbers: ", InterfaceElementNodes

#------------------------------------------------------------------

# Read interface nodes for solid/fluid/interface
if (MaterialSpecification == Plate2D):
    # Read xi position's
    f = './Input/2D/Solid/Plate2DSolidXi2.csv'
    with open(f,'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        rownum = 0
        for row in reader:
            if (rownum == 0):
                # Read the header row
                header = row
                SolidXi2 = int(row[0])*[3*[0]]
            else:
                SolidXi2[rownum-1] = row
            # Next line
            rownum+=1
        if (FileReadDiagnostics): print "Solid xi2: ", SolidXi2
#------------------------------------------------------------------
    f = './Input/2D/Solid/Plate2DSolidXi3.csv'
    with open(f,'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        rownum = 0
        for row in reader:
            if (rownum == 0):
                # Read the header row
                header = row
                SolidXi3 = int(row[0])*[3*[0]]
            else:
                SolidXi3[rownum-1] = row
            # Next line
            rownum+=1
        if (FileReadDiagnostics): print "Solid xi3: ", SolidXi3
#------------------------------------------------------------------
    f = './Input/2D/Fluid/Plate2DFluidXi2.csv'
    with open(f,'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        rownum = 0
        for row in reader:
            if (rownum == 0):
                # Read the header row
                header = row
                FluidXi2 = int(row[0])*[3*[0]]
            else:
                FluidXi2[rownum-1] = row
            # Next line
            rownum+=1
        if (FileReadDiagnostics): print "Fluid xi2: ", FluidXi2
#------------------------------------------------------------------
    f = './Input/2D/Fluid/Plate2DFluidXi3.csv'
    with open(f,'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        rownum = 0
        for row in reader:
            if (rownum == 0):
                # Read the header row
                header = row
                FluidXi3 = int(row[0])*[3*[0]]
            else:
                FluidXi3[rownum-1] = row
            # Next line
            rownum+=1
        if (FileReadDiagnostics): print "Fluid xi3: ", FluidXi3
#------------------------------------------------------------------
    f = './Input/2D/Interface/Plate2DInterfaceSolidElements.csv'
    with open(f,'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        rownum = 0
        for row in reader:
            if (rownum == 0):
                # Read the header row
                header = row
                SolidInterfaceElements = int(row[0])*[0]
            else:
                SolidInterfaceElements = row
            # Next line
            rownum+=1
    if (FileReadDiagnostics): print "Solid interface element numbers: ", SolidInterfaceElements
#------------------------------------------------------------------
    f = './Input/2D/Interface/Plate2DInterfaceFluidElements.csv'
    with open(f,'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        rownum = 0
        for row in reader:
            if (rownum == 0):
                # Read the header row
                header = row
                FluidInterfaceElements = int(row[0])*[0]
            else:
                FluidInterfaceElements = row
            # Next line
            rownum+=1
    if (FileReadDiagnostics): print "Fluid interface element numbers: ", FluidInterfaceElements
#------------------------------------------------------------------
elif (MaterialSpecification == Plate3D):

    f = './Input/3D/Solid/Plate3DSolidInterfaceNodeInformationNodeElement.csv'
    line = f.read().splitlines()
    SolidInterfaceNodeInformationNE = len(line)*[2*[0]]
    for s in range(0,len(line)):
        line1 = line[s].split(" ")
        SolidInterfaceNodeInformationNE[s]=line1[:]
    f.close()
    if (FileReadDiagnostics): print "Fluid interface element numbers: ", SolidInterfaceNodeInformationNE

    f = './Input/3D/Fluid/Plate3DFluidInterfaceNodeInformationNodeElement.csv'
    line = f.read().splitlines()
    FluidInterfaceNodeInformationNE = len(line)*[2*[0]]
    for s in range(0,len(line)):
        line1 = line[s].split(" ")
        FluidInterfaceNodeInformationNE[s]=line1[:]
    f.close()
    if (FileReadDiagnostics): print "Fluid interface element numbers: ", FluidInterfaceNodeInformationNE
    
    f = './Input/3D/Interface/Plate3DIInterfaceNodeInformationNodeElement.csv'
    line = f.read().splitlines()
    InterfaceInterfaceNodeInformationNE = len(line)*[2*[0]]
    for s in range(0,len(line)):
        line1 = line[s].split(" ")
        InterfaceInterfaceNodeInformationNE[s]=line1[:]
    f.close()
    if (FileReadDiagnostics): print "Fluid interface element numbers: ", InterfaceInterfaceNodeInformationNE
    
    f = './Input/3D/Solid/Plate3DSolidInterfaceNodeInformationXi.csv'
    line = f.read().splitlines()
    SolidInterfaceNodeInformationXi = len(line)*[3*[0]]
    for s in range(0,len(line)):
        line1 = line[s].split(" ")
        SolidInterfaceNodeInformationXi[s]=line1[:]
    f.close()
    if (FileReadDiagnostics): print "Fluid interface element numbers: ", SolidInterfaceNodeInformationXi

    f = './Input/3D/Fluid/Plate3DFluidInterfaceNodeInformationXi.csv'
    line = f.read().splitlines()
    FluidInterfaceNodeInformationXi = len(line)*[3*[0]]
    for s in range(0,len(line)):
        line1 = line[s].split(" ")
        FluidInterfaceNodeInformationXi[s]=line1[:]
    f.close()
    if (FileReadDiagnostics): print "Fluid interface element numbers: ", FluidInterfaceNodeInformationXi
        
    f = './Input/3D/BC/Plate3DlagrangeNodes.csv'
    LagrangeNodes = len(line)*[0]
    for s in range(0,len(line)):
        line1 = line[s].split(" ")
        LagrangeNodes[s]=line1[:]
    f.close()
    if (FileReadDiagnostics): print "Fluid interface element numbers: ", LagrangeNodes

#------------------------------------------------------------------

# Read connected solid/fluid nodes for sorted interface nodes
if (MaterialSpecification == Plate2D):
    f = './Input/2D/Interface/Plate2DSortedInterfaceNodes.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/Interface/Plate3DSortedInterfaceNodes.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            ConnectedInterfaceNodes = 2*[int(row[0])*[0]]
        else:
            ConnectedInterfaceNodes[rownum-1] = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Connected nodes: ", ConnectedInterfaceNodes

#------------------------------------------------------------------

# Read boundary conditions
if (MaterialSpecification == Plate2D):
    f = './Input/2D/BC/Plate2DdisplacementBC.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/BC/Plate3DdisplacementBC.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            NoDisplacementNodes = int(row[0])*[0]
        else:
            row = [int(i) for i in row]
            NoDisplacementNodes = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Solid node numbers with no displacement: ", NoDisplacementNodes  

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/BC/Plate2DfixedNodesBC.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/BC/Plate3DfixedNodesBC.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            FixedNodes = int(row[0])*[0]
        else:
            row = [int(i) for i in row]
            FixedNodes = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "In both coordinates fixed mesh nodes: ", FixedNodes  
  
#------------------------------------------------------------------
  
if (MaterialSpecification == Plate2D):
    f = './Input/2D/BC/Plate2DfixedZNodesBC.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/BC/Plate3DfixedZNodesBC.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            FixedZNodes = int(row[0])*[0]
        else:
            row = [int(i) for i in row]
            FixedZNodes = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "In 2nd coordinates fixed mesh nodes: ", FixedZNodes  

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/BC/Plate2DmovedNodesBC.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/BC/Plate3DmovedNodesBC.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            MovedNodes = int(row[0])*[0]
        else:
            row = [int(i) for i in row]
            MovedNodes = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "In both coordinates moving mesh nodes: ", MovedNodes  

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/BC/Plate2DmovedYNodesBC.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/BC/Plate3DmovedYNodesBC.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            MovedYNodes = int(row[0])*[0]
        else:
            row = [int(i) for i in row]
            MovedYNodes = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "In 1st coordinate moving mesh nodes: ", MovedYNodes 

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/BC/Plate2DinletBC.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/BC/Plate3DinletBC.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            InletNodes = int(row[0])*[0]
        else:
            row = [int(i) for i in row]
            InletNodes = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Inlet nodes: ", InletNodes 

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/BC/Plate2DnoSlipBC.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/BC/Plate3DnoSlipBC.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            NoSlipNodes = int(row[0])*[0]
        else:
            row = [int(i) for i in row]
            NoSlipNodes = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "No-slip nodes: ", NoSlipNodes 

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/BC/Plate2DslipBC.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/BC/Plate3DslipTop.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            SlipNodesTop = int(row[0])*[0]
        else:
            row = [int(i) for i in row]
            SlipNodesTop = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Slip nodes top: ", SlipNodesTop 

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/BC/Plate2DslipBC.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/BC/Plate3DslipRightLeft.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            SlipNodesRightLeft = int(row[0])*[0]
        else:
            row = [int(i) for i in row]
            SlipNodesRightLeft = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Slip nodes r/l: ", SlipNodesRightLeft 

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/BC/Plate2DpressureBC.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/BC/Plate3DpressureBC.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            OutletNodes = int(row[0])*[0]
        else:
            row = [int(i) for i in row]
            OutletNodes = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Outlet nodes: ", OutletNodes

#------------------------------------------------------------------

if (MaterialSpecification == Plate2D):
    f = './Input/2D/BC/Plate2DslipBC.csv'
elif (MaterialSpecification == Plate3D):
    f = './Input/3D/BC/Plate3DslipBC.csv'
with open(f,'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum = 0
    for row in reader:
        if (rownum == 0):
            # Read the header row
            header = row
            SlipNodes = int(row[0])*[0]
        else:
            row = [int(i) for i in row]
            SlipNodes = row
        # Next line
        rownum+=1
if (FileReadDiagnostics): print "Slip BC nodes: ", SlipNodes

#------------------------------------------------------------------

print  "Finished reading geometry."

if (SetupOutput):
    print "SUMMARY"
    print " "
    print "Start time: ", StartTime
    print "Stop time: ", StopTime
    print "Time increment: ", TimeStepSize
    print " "
    print "FLUID:"
    print "Dynamic viscosity: ", FluidDynamicViscosity ,"kg/(m s)"
    print "Density: ", FluidDensity ,"kg m^-3"
    if (NumberOfDimensions == 3):
        print "Domain:", FluidGeometryX[NumberOfFluidNodes-1] ,"x", FluidGeometryY[NumberOfFluidNodes-1] ,"x", FluidGeometryZ[NumberOfFluidNodes-1] ," m^3"
    else:
        print "Domain:", FluidGeometryY[NumberOfFluidNodes-1] ,"x", FluidGeometryZ[NumberOfFluidNodes-1] ," m^2"
    print  "Number of nodes: ", NumberOfFluidNodes
    print " "
    print "SOLID:"
    print "Density: ", SolidDensity ,"kg m^-3"
    print "Young's modulus: ", YoungsModulus ,"Pa"
    print "Poisson's ratio: ", PoissonsRatio
    print "Neo-Hookean constant: ", MooneyRivlin1
    if (NumberOfDimensions == 3):
        print "Domain:", SolidGeometryX[NumberOfSolidNodes-1]-SolidGeometryX[1] ,"x", SolidGeometryY[NumberOfSolidNodes-1]-SolidGeometryY[1] ,"x", SolidGeometryZ[NumberOfSolidNodes-1]-SolidGeometryZ[1] ,"m^3"
    else:
        print "Domain:", SolidGeometryY[NumberOfSolidNodes-1]-SolidGeometryY[1] ,"x", SolidGeometryZ[NumberOfSolidNodes-1]-SolidGeometryZ[1] ,"m^2"
    print  "Number of nodes: ", NumberOfSolidNodes
    print " "
    print  "Number of interface nodes: ", NumberOfInterfaceNodes

#================================================================================================================================
#  Coordinate System
#================================================================================================================================

# Create a RC coordinate system for the solid region
if (ExampleFileProgressDiagnostics):
    print " == >> SOLID COORDINATE SYSTEM << == "
SolidCoordinateSystem = iron.CoordinateSystem()
SolidCoordinateSystem.CreateStart(SolidCoordinateSystemUserNumber)
SolidCoordinateSystem.DimensionSet(NumberOfDimensions)
SolidCoordinateSystem.CreateFinish()
# Create a RC coordinate system for the fluid region
if (ExampleFileProgressDiagnostics):
    print " == >> FLUID COORDINATE SYSTEM << == "
FluidCoordinateSystem = iron.CoordinateSystem()
FluidCoordinateSystem.CreateStart(FluidCoordinateSystemUserNumber)
FluidCoordinateSystem.DimensionSet(NumberOfDimensions)
FluidCoordinateSystem.CreateFinish()
# Create a RC coordinate system for the interface region
if (ExampleFileProgressDiagnostics):
    print " == >> INTERFACE COORDINATE SYSTEM << == "
InterfaceCoordinateSystem = iron.CoordinateSystem()
InterfaceCoordinateSystem.CreateStart(InterfaceCoordinateSystemUserNumber)
InterfaceCoordinateSystem.DimensionSet(NumberOfDimensions)
InterfaceCoordinateSystem.CreateFinish()
  
#================================================================================================================================
#  Region
#================================================================================================================================

# Create a solid region
if (ExampleFileProgressDiagnostics):
    print " == >> SOLID REGION << == "
Region1 = iron.Region()
Region1.CreateStart(SolidRegionUserNumber,iron.WorldRegion)
Region1.label = "SolidRegion"
Region1.coordinateSystem = SolidCoordinateSystem
Region1.CreateFinish()
# Create a fluid region
if (ExampleFileProgressDiagnostics):
    print " == >> FLUID REGION << == "
Region2 = iron.Region()
Region2.CreateStart(FluidRegionUserNumber,iron.WorldRegion)
Region2.label = "FluidRegion"
Region2.coordinateSystem = FluidCoordinateSystem
Region2.CreateFinish()

#================================================================================================================================
#  Bases
#================================================================================================================================

# Create the solid basis
if (ExampleFileProgressDiagnostics):
    print " == >> BASIS FOR SOLID: DISPLACEMENT, HYDROSTATIC PRESSURE << == "
# Create a basis for the dependent field variable displacement
BasisDisplacement = iron.Basis()
BasisDisplacement.CreateStart(BasisDisplacementUserNumber)
if (InterpolationTypeDisplacement == 1 or InterpolationTypeDisplacement == 2 or 
    InterpolationTypeDisplacement == 3 or InterpolationTypeDisplacement == 4):
    BasisDisplacement.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
elif (InterpolationTypeDisplacement == 7 or InterpolationTypeDisplacement == 8 or 
      InterpolationTypeDisplacement == 9):
    BasisDisplacement.type = iron.BasisTypes.SIMPLEX_TYPE
if (InterpolationTypeDisplacement == iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE):
    NumberOfGaussXi = 2
    PressureMeshComponent = 1
elif(InterpolationTypeDisplacement == iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE):
    NumberOfGaussXi = 3
    PressureMeshComponent = 2
    InterpolationTypeHydrostaticPressure = 1
elif(InterpolationTypeDisplacement == iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE or
     InterpolationTypeDisplacement == iron.BasisInterpolationSpecifications.CUBIC_HERMITE):
    NumberOfGaussXi = 4
    PressureMeshComponent = 2
    InterpolationTypeHydrostaticPressure = 2
else:
    NumberOfGaussXi = 0
    PressureMeshComponent = 1
BasisDisplacement.numberOfXi = NumberOfDimensions
if (NumberOfDimensions == 2):
    BasisDisplacement.interpolationXi = [InterpolationTypeDisplacement,
                                         InterpolationTypeDisplacement]
    BasisDisplacement.quadratureNumberOfGaussXi = [NumberOfGaussXi,
                                                   NumberOfGaussXi]
else:
    BasisDisplacement.interpolationXi = [InterpolationTypeDisplacement,
                                         InterpolationTypeDisplacement,
                                         InterpolationTypeDisplacement]
    BasisDisplacement.quadratureNumberOfGaussXi = [NumberOfGaussXi,
                                                   NumberOfGaussXi,
                                                   NumberOfGaussXi]
BasisDisplacement.CreateFinish()
# Use the displacement basis as a space basis
BasisSpaceSolid = BasisDisplacement
# Create a basis for the dependent field variable hydrostatic pressure
BasisHydrostaticPressure = iron.Basis()
BasisHydrostaticPressure.CreateStart(BasisHydrostaticPressureUserNumber)
if (InterpolationTypeHydrostaticPressure == 1 or InterpolationTypeHydrostaticPressure == 2 or 
    InterpolationTypeHydrostaticPressure == 3 or InterpolationTypeHydrostaticPressure == 4):
    BasisHydrostaticPressure.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
elif (InterpolationTypeHydrostaticPressure == 7 or InterpolationTypeHydrostaticPressure == 8 or 
      InterpolationTypeHydrostaticPressure == 9):
    BasisHydrostaticPressure.type = iron.BasisTypes.SIMPLEX
BasisHydrostaticPressure.numberOfXi = NumberOfDimensions
if (NumberOfDimensions == 2):
    BasisHydrostaticPressure.interpolationXi = [InterpolationTypeHydrostaticPressure,
                                                InterpolationTypeHydrostaticPressure]
    BasisHydrostaticPressure.quadratureNumberOfGaussXi = [NumberOfGaussXi,
                                                          NumberOfGaussXi]
else:
    BasisHydrostaticPressure.interpolationXi = [InterpolationTypeHydrostaticPressure,
                                                InterpolationTypeHydrostaticPressure,
                                                InterpolationTypeHydrostaticPressure]
    BasisHydrostaticPressure.quadratureNumberOfGaussXi = [NumberOfGaussXi,
                                                          NumberOfGaussXi,
                                                          NumberOfGaussXi]
BasisHydrostaticPressure.CreateFinish()

#---------------------------------------

# Create the fluid basis
if (ExampleFileProgressDiagnostics):
    print " == >> BASIS FOR FLUID: SPACE, VELOCITY, PRESSURE << == "
# Create a basis for the fluid domain
MeshNumberOfComponents = 1
BasisSpaceFluid = iron.Basis()
BasisSpaceFluid.CreateStart(BasisSpaceFluidUserNumber)
if (InterpolationTypeSpace == 1 or InterpolationTypeSpace == 2 or 
    InterpolationTypeSpace == 3 or InterpolationTypeSpace == 4):
    BasisSpaceFluid.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
elif (InterpolationTypeSpace == 7 or InterpolationTypeSpace == 8 or 
      InterpolationTypeSpace == 9):
    BasisSpaceFluid.type = iron.BasisTypes.SIMPLEX
if (InterpolationTypeSpace == 2):
    NumberOfGaussXiSpace = 3
BasisSpaceFluid.numberOfXi = NumberOfDimensions
if (NumberOfDimensions == 2):
    BasisSpaceFluid.interpolationXi = [InterpolationTypeSpace,
                                       InterpolationTypeSpace]
    BasisSpaceFluid.quadratureNumberOfGaussXi = [NumberOfGaussXiSpace,
                                                 NumberOfGaussXiSpace]
else:
    BasisSpaceFluid.interpolationXi = [InterpolationTypeSpace,
                                       InterpolationTypeSpace,
                                       InterpolationTypeSpace]
    BasisSpaceFluid.quadratureNumberOfGaussXi = [NumberOfGaussXiSpace,
                                                 NumberOfGaussXiSpace,
                                                 NumberOfGaussXiSpace]
BasisSpaceFluid.CreateFinish()
# Create a basis for the dependent field variable velocity
if (InterpolationTypeVelocity == InterpolationTypeSpace):
    BasisVelocity = BasisSpaceFluid
else:
    MeshNumberOfComponents = MeshNumberOfComponents+1
    BasisVelocity = iron.Basis()
    BasisVelocity.CreateStart(BasisVelocityUserNumber)
    BasisVelocity.type = iron.BasisTypes.CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE
    NumberOfGaussXiVelocity = 3
    BasisVelocity.numberOfXi = NumberOfDimensions
    if (NumberOfDimensions == 2):
        BasisVelocity.interpolationXi = [InterpolationTypeVelocity,
                                         InterpolationTypeVelocity]
        BasisVelocity.quadratureNumberOfGaussXi = [NumberOfGaussXiVelocity,
                                                   NumberOfGaussXiVelocity]
    else:
        BasisVelocity.interpolationXi = [InterpolationTypeVelocity,
                                         InterpolationTypeVelocity,
                                         InterpolationTypeVelocity]
        BasisVelocity.quadratureNumberOfGaussXi = [NumberOfGaussXiVelocity,
                                                   NumberOfGaussXiVelocity,
                                                   NumberOfGaussXiVelocity]
    BasisVelocity.CreateFinish()
# Create a basis for the dependent field variable pressure
if (InterpolationTypePressure == InterpolationTypeSpace):
    BasisPressure=BasisSpaceFluid
elif (InterpolationTypePressure == InterpolationTypeVelocity):
    BasisPressure = BasisVelocity
else:
    MeshNumberOfComponents = MeshNumberOfComponents+1
    NumberOfGaussXiPressure = 3
    BasisPressure = iron.Basis()
    BasisPressure.CreateStart(BasisPressureUserNumber)
    BasisPressure.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
    BasisPressure.numberOfXi = NumberOfDimensions
    if (NumberOfDimensions == 2):
        BasisPressure.interpolationXi = [InterpolationTypePressure,
                                         InterpolationTypePressure]
        BasisPressure.quadratureNumberOfGaussXi = [NumberOfGaussXiPressure,
                                                   NumberOfGaussXiPressure]
    else:
        BasisPressure.interpolationXi = [InterpolationTypePressure,
                                         InterpolationTypePressure,
                                         InterpolationTypePressure]
        BasisPressure.quadratureNumberOfGaussXi = [NumberOfGaussXiPressure,
                                                   NumberOfGaussXiPressure,
                                                   NumberOfGaussXiPressure] 
    BasisPressure.CreateFinish()
  
#================================================================================================================================
#  Mesh
#================================================================================================================================

# Create the solid mesh
if (ExampleFileProgressDiagnostics):
    print " == >> SOLID MESH << == "
# Start the creation of mesh nodes
SolidNodes = iron.Nodes()
SolidNodes.CreateStart(Region1,NumberOfSolidNodes)
SolidNodes.CreateFinish()

Mesh1 = iron.Mesh()
Mesh1.CreateStart(SolidMeshUserNumber,Region1,NumberOfDimensions)
Mesh1.NumberOfElementsSet(NumberOfSolidElements)
Mesh1.NumberOfComponentsSet(2)

SolidMeshElementsSpace               = iron.MeshElements()
SolidMeshElementsDisplacement        = iron.MeshElements()
SolidMeshElementsHydrostaticPressure = iron.MeshElements()

Mesh1ComponentNumberSpace               = 1
Mesh1ComponentNumberDisplacement        = 1
Mesh1ComponentNumberHydrostaticPressure = 1

SolidMeshElementsSpace.CreateStart(Mesh1,Mesh1ComponentNumberSpace,BasisSpaceSolid)
for ElementIndex in range(1,NumberOfSolidElements+1):
    SolidMeshElementsSpace.NodesSet(ElementIndex,SolidElementNodes[ElementIndex-1])
SolidMeshElementsSpace.CreateFinish()

SolidMeshElementsDisplacement = SolidMeshElementsSpace
Mesh1ComponentNumberHydrostaticPressure = Mesh1ComponentNumberDisplacement+1

SolidMeshElementsHydrostaticPressure.CreateStart(Mesh1,Mesh1ComponentNumberHydrostaticPressure,BasisHydrostaticPressure)
if (MaterialSpecification == Plate2D):
    for ElementIndex in range(1,NumberOfSolidElements+1):
        SolidMeshElementsHydrostaticPressure.NodesSet(ElementIndex,[SolidElementNodes[ElementIndex-1][0],
                                                                    SolidElementNodes[ElementIndex-1][2],
                                                                    SolidElementNodes[ElementIndex-1][6],
                                                                    SolidElementNodes[ElementIndex-1][8]])
elif (MaterialSpecification == Plate3D):
    for ElementIndex in range(1,NumberOfSolidElements+1):
        SolidMeshElementsHydrostaticPressure.NodesSet(ElementIndex,[SolidElementNodes[ElementIndex-1][0],
                                                                    SolidElementNodes[ElementIndex-1][2],
                                                                    SolidElementNodes[ElementIndex-1][6],
                                                                    SolidElementNodes[ElementIndex-1][8],
                                                                    SolidElementNodes[ElementIndex-1][18],
                                                                    SolidElementNodes[ElementIndex-1][20],
                                                                    SolidElementNodes[ElementIndex-1][24],
                                                                    SolidElementNodes[ElementIndex-1][26]])
SolidMeshElementsHydrostaticPressure.CreateFinish()
Mesh1.CreateFinish()

#-----------------------------------

# Create the fluid mesh
if (ExampleFileProgressDiagnostics):
    print " == >> FLUID MESH << == "
# Start the creation of mesh nodes
FluidNodes = iron.Nodes()
FluidNodes.CreateStart(Region2,NumberOfFluidNodes)
FluidNodes.CreateFinish()

Mesh2 = iron.Mesh()
Mesh2.CreateStart(FluidMeshUserNumber,Region2,NumberOfDimensions)
Mesh2.NumberOfElementsSet(NumberOfFluidElements)
Mesh2.NumberOfComponentsSet(2)

FluidMeshElementsSpace    = iron.MeshElements()
FluidMeshElementsVelocity = iron.MeshElements()
FluidMeshElementsPressure = iron.MeshElements()

Mesh2ComponentNumberSpace        = 1
FluidMeshComponentNumberVelocity = 1
Mesh2ComponentNumberPressure     = 1

FluidMeshElementsSpace.CreateStart(Mesh2,Mesh2ComponentNumberSpace,BasisSpaceFluid)
for ElementIndex in range(1,NumberOfFluidElements+1):
    FluidMeshElementsSpace.NodesSet(ElementIndex,FluidElementNodes[ElementIndex-1])
FluidMeshElementsSpace.CreateFinish()

FluidMeshElementsVelocity    = FluidMeshElementsSpace
Mesh2ComponentNumberPressure = FluidMeshComponentNumberVelocity+1

FluidMeshElementsPressure.CreateStart(Mesh2,Mesh2ComponentNumberPressure,BasisPressure)
if (MaterialSpecification == Plate2D):
    for ElementIndex in range(1,NumberOfFluidElements+1):
        FluidMeshElementsPressure.NodesSet(ElementIndex,[FluidElementNodes[ElementIndex-1][0],
                                                         FluidElementNodes[ElementIndex-1][2],
                                                         FluidElementNodes[ElementIndex-1][6],
                                                         FluidElementNodes[ElementIndex-1][8]])
elif (MaterialSpecification == Plate3D):
    for ElementIndex in range(1,NumberOfFluidElements+1):
        FluidMeshElementsPressure.NodesSet(ElementIndex,[FluidElementNodes[ElementIndex-1][0],
                                                         FluidElementNodes[ElementIndex-1][2],
                                                         FluidElementNodes[ElementIndex-1][6],
                                                         FluidElementNodes[ElementIndex-1][8],
                                                         FluidElementNodes[ElementIndex-1][18],
                                                         FluidElementNodes[ElementIndex-1][20],
                                                         FluidElementNodes[ElementIndex-1][24],
                                                         FluidElementNodes[ElementIndex-1][26]])
FluidMeshElementsPressure.CreateFinish()
Mesh2.CreateFinish()

#================================================================================================================================
#  Interface
#================================================================================================================================

# Create an interface between the two meshes
if (ExampleFileProgressDiagnostics):
    print " == >> INTERFACE << == "
Interface1 = iron.Interface()
Interface1.CreateStart(InterfaceUserNumber,iron.WorldRegion)
Interface1.LabelSet("Interface1")
# Add in the two meshes
Interface1.MeshAdd(Mesh1)
Interface1.MeshAdd(Mesh2)
Interface1.CoordinateSystemSet(InterfaceCoordinateSystem)
Interface1.CreateFinish()

# Create a (bi)-quadratic-Lagrange basis (3D: faces // 2D: lines)
if (ExampleFileProgressDiagnostics):
    print " == >> INTERFACE BASIS << == "
InterfaceBasis1 = iron.Basis()
InterfaceBasis1.CreateStart(InterfaceBasisUserNumber)
InterfaceBasis1.NumberOfXiSet(NumberOfDimensions-1)
if (NumberOfDimensions == 2):
    InterfaceBasis1.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]
else:
    InterfaceBasis1.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE,
                                       iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]
InterfaceBasis1.CreateFinish()

# Create a (bi)-quadratic-Lagrange basis for the interface mapping (3D: faces // 2D: lines)
if (ExampleFileProgressDiagnostics):
    print " == >> INTERFACE MAPPING BASIS << == "
InterfaceMappingBasis1 = iron.Basis()
InterfaceMappingBasis1.CreateStart(InterfaceMappingBasisUserNumber)
InterfaceMappingBasis1.NumberOfXiSet(NumberOfDimensions-1)
if (NumberOfDimensions == 2):
    InterfaceMappingBasis1.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]
else:
    InterfaceMappingBasis1.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE,
                                              iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]
InterfaceMappingBasis1.CreateFinish()
  
#================================================================================================================================
#  Interface Mesh
#================================================================================================================================

# Create an interface mesh
if (ExampleFileProgressDiagnostics):
    print " == >> INTERFACE MESH << == "
InterfaceNodes = iron.Nodes()
InterfaceNodes.CreateStartInterface(Interface1,NumberOfInterfaceNodes)
InterfaceNodes.CreateFinish()

InterfaceMesh1 = iron.Mesh()
InterfaceMesh1.CreateStartInterface(InterfaceMeshUserNumber,Interface1,NumberOfDimensions-1)
InterfaceMesh1.NumberOfElementsSet(NumberOfInterfaceElements)
InterfaceMesh1.NumberOfComponentsSet(1)

InterfaceMeshElements = iron.MeshElements()
InterfaceMeshComponentNumber = 1
InterfaceMeshElements.CreateStart(InterfaceMesh1,InterfaceMeshComponentNumber,InterfaceBasis1)
for ElementIndex in range(1,NumberOfInterfaceElements+1):
    InterfaceMeshElements.NodesSet(ElementIndex,InterfaceElementNodes[ElementIndex-1])
InterfaceMeshElements.CreateFinish()

InterfaceMesh1.CreateFinish()

#================================================================================================================================
#  Mesh Connectivity
#================================================================================================================================

# Couple the interface meshes
if (ExampleFileProgressDiagnostics):
    print " == >> INTERFACE MESH CONNECTIVITY << == "
InterfaceMeshConnectivity1 = iron.InterfaceMeshConnectivity()
InterfaceMeshConnectivity1.CreateStart(Interface1,InterfaceMesh1)
InterfaceMeshConnectivity1.BasisSet(InterfaceMappingBasis1)

if (NumberOfDimensions == 2):
    for ElementIndex in range(1,NumberOfInterfaceElements+1):
        # Map the interface element to the elements in mesh 1
        InterfaceMeshConnectivity1.ElementNumberSet(ElementIndex,SolidMeshIndex,int(SolidInterfaceElements[ElementIndex-1]))
        for LocalNodeIndex in range(1,3**(NumberOfDimensions-1)+1):
            XI2 = [SolidXi2[ElementIndex-1][LocalNodeIndex-1],SolidXi3[ElementIndex-1][LocalNodeIndex-1]]
            InterfaceMeshConnectivity1.ElementXiSet(ElementIndex,SolidMeshIndex,int(SolidInterfaceElements[ElementIndex-1]),
                                                    LocalNodeIndex,1,XI2)
        # Map the interface element to the elements in mesh 2
        InterfaceMeshConnectivity1.ElementNumberSet(ElementIndex,FluidMeshIndex,int(FluidInterfaceElements[ElementIndex-1]))
        for LocalNodeIndex in range(1,3**(NumberOfDimensions-1)+1):
            XI2 = [FluidXi2[ElementIndex-1][LocalNodeIndex-1],FluidXi3[ElementIndex-1][LocalNodeIndex-1]]
            InterfaceMeshConnectivity1.ElementXiSet(ElementIndex,FluidMeshIndex,int(FluidInterfaceElements[ElementIndex-1]),
                                                    LocalNodeIndex,1,XI2)
    InterfaceMeshConnectivity1.NodeNumberSet(InterfaceNodeNumbersForGeometry,SolidMeshIndex,ConnectedInterfaceNodes[0],
                                             FluidMeshIndex,ConnectedInterfaceNodes[1])
else:
    LocalNodeIndex=0
    for ElementIndex in range(1,9*NumberOfInterfaceElements+1):
        LocalNodeIndex = LocalNodeIndex+1
        #pdb.set_trace()
        # Map the interface element to the elements in the solid mesh
        InterfaceMeshConnectivity1.ElementNumberSet(int(InterfaceInterfaceNodeInformationNE[ElementIndex][1]),SolidMeshIndex,
                                                    int(SolidInterfaceNodeInformationNE[ElementIndex][1]))
        XI3 = [SolidInterfaceNodeInformationXi[ElementIndex][0],SolidInterfaceNodeInformationXi[ElementIndex][1],
               SolidInterfaceNodeInformationXi[ElementIndex][2]]
        InterfaceMeshConnectivity1.ElementXiSet(int(InterfaceInterfaceNodeInformationNE[ElementIndex][1]),SolidMeshIndex,
                                                int(SolidInterfaceNodeInformationNE[ElementIndex][1]),LocalNodeIndex,1,XI3)
        # Map the interface element to the elements in the fluid mesh
        InterfaceMeshConnectivity1.ElementNumberSet(int(InterfaceInterfaceNodeInformationNE[ElementIndex][1]),FluidMeshIndex,
                                                    int(FluidInterfaceNodeInformationNE[ElementIndex][1]))
        XI3 = [FluidInterfaceNodeInformationXi[ElementIndex][0],FluidInterfaceNodeInformationXi[ElementIndex][1],
               FluidInterfaceNodeInformationXi[ElementIndex][2]]
        InterfaceMeshConnectivity1.ElementXiSet(int(InterfaceInterfaceNodeInformationNE[ElementIndex][1]),FluidMeshIndex,
                                                int(FluidInterfaceNodeInformationNE[ElementIndex][1]),LocalNodeIndex,1,XI3)
        if (LocalNodeIndex == 9): 
            LocalNodeIndex = 0
    InterfaceMeshConnectivity1.NodeNumberSet(InterfaceNodeNumbersForGeometry[1],SolidMeshIndex,ConnectedInterfaceNodes[1],
                                             FluidMeshIndex,ConnectedInterfaceNodes[2])
InterfaceMeshConnectivity1.CreateFinish()

#================================================================================================================================
#  Decomposition
#================================================================================================================================

# Create a decomposition for the solid mesh
if (ExampleFileProgressDiagnostics):
    print " == >> SOLID MESH DECOMPOSITION << == "
SolidDecomposition = iron.Decomposition()
SolidDecomposition.CreateStart(SolidDecompositionUserNumber,Mesh1)
SolidDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
SolidDecomposition.NumberOfDomainsSet(NumberOfComputationalNodes)
SolidDecomposition.CalculateFacesSet(True)
SolidDecomposition.CreateFinish()

# Create a decomposition for the fluid mesh
if (ExampleFileProgressDiagnostics):
    print " == >> FLUID MESH DECOMPOSITION << == "
FluidDecomposition = iron.Decomposition()
FluidDecomposition.CreateStart(FluidDecompositionUserNumber,Mesh2)
FluidDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
FluidDecomposition.NumberOfDomainsSet(NumberOfComputationalNodes)
FluidDecomposition.CalculateFacesSet(True)
FluidDecomposition.CreateFinish()

# Create a decomposition for the interface mesh
if (ExampleFileProgressDiagnostics):
    print " == >> INTERFACE DECOMPOSITION << == "
InterfaceDecomposition = iron.Decomposition()
InterfaceDecomposition.CreateStart(SolidDecompositionUserNumber,InterfaceMesh1)
InterfaceDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
InterfaceDecomposition.NumberOfDomainsSet(NumberOfComputationalNodes)
InterfaceDecomposition.CreateFinish()

#================================================================================================================================
#  Geometric Field
#================================================================================================================================

# Start to create a default (geometric) field on the solid region
if (ExampleFileProgressDiagnostics):
    print " == >> SOLID MESH GEOMETRIC FIELD << == "
GeometricField1 = iron.Field()
GeometricField1.CreateStart(SolidGeometricFieldUserNumber,Region1)
# Set the decomposition to use
GeometricField1.meshDecomposition = SolidDecomposition
# Set the scaling to use
GeometricField1.ScalingTypeSet(iron.FieldScalingTypes.NONE)
GeometricField1.VariableLabelSet(iron.FieldVariableTypes.U,'SolidGF')
# Set the domain to be used by the field components.
GeometricField1.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,Mesh1ComponentNumberSpace)
GeometricField1.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,Mesh1ComponentNumberSpace)
if (NumberOfDimensions == 3):
    GeometricField1.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,Mesh1ComponentNumberSpace)
# Finish creating the first field
GeometricField1.CreateFinish()

# Start to create a default (geometric) field on the fluid region
if (ExampleFileProgressDiagnostics):
    print " == >> FLUID MESH GEOMETRIC FIELD << == "  
GeometricField2 = iron.Field()
GeometricField2.CreateStart(FluidGeometricFieldUserNumber,Region2)
# Set the decomposition to use
GeometricField2.MeshDecompositionSet(FluidDecomposition)
# Set the scaling to use
GeometricField2.ScalingTypeSet(iron.FieldScalingTypes.NONE)
GeometricField2.VariableLabelSet(iron.FieldVariableTypes.U,'FluidGF')
# Set the domain to be used by the field components.
GeometricField2.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,Mesh2ComponentNumberSpace)
GeometricField2.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,Mesh2ComponentNumberSpace)
if (NumberOfDimensions == 3):
    GeometricField2.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,Mesh2ComponentNumberSpace)
# Finish creating the second field
GeometricField2.CreateFinish()

# Start to create a default (geometric) field on the Interface
if (ExampleFileProgressDiagnostics):
    print " == >> INTERFACE GEOMETRIC FIELD << == "  
InterfaceGeometricField1 = iron.Field()
InterfaceGeometricField1.CreateStartInterface(InterfaceGeometricFieldUserNumber,Interface1)
# Set the decomposition to use
InterfaceGeometricField1.MeshDecompositionSet(InterfaceDecomposition)
# Set the scaling to use
InterfaceGeometricField1.ScalingTypeSet(iron.FieldScalingTypes.NONE)
InterfaceGeometricField1.VariableLabelSet(iron.FieldVariableTypes.U,'InterfaceGF')
# Set the domain to be used by the field components.
InterfaceGeometricField1.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,InterfaceMeshComponentNumber)
InterfaceGeometricField1.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,InterfaceMeshComponentNumber)
if (NumberOfDimensions == 3):
    InterfaceGeometricField1.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,InterfaceMeshComponentNumber)
# Finish creating the first field
InterfaceGeometricField1.CreateFinish()

# Update the geometric field parameters (solid)
for NodeIndex in range(1,NumberOfSolidNodes+1):
    NodeDomain = SolidDecomposition.NodeDomainGet(SolidNodeNumbers[NodeIndex-1],Mesh1ComponentNumberSpace)
    if (NodeDomain == ComputationalNodeNumber):
        if (NumberOfDimensions == 3):
            Components = [2,3,1]
        else:
            Components = [1,2,0]
        GeometricField1.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                 versionIdx,derivIdx,SolidNodeNumbers[NodeIndex-1],
                                                 Components[0],SolidGeometryY[NodeIndex-1])
        GeometricField1.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                 versionIdx,derivIdx,SolidNodeNumbers[NodeIndex-1],
                                                 Components[1],SolidGeometryZ[NodeIndex-1])
        if (NumberOfDimensions == 3):
            GeometricField1.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                     versionIdx,derivIdx,SolidNodeNumbers[NodeIndex-1],
                                                     Components[2],SolidGeometryX[NodeIndex])
# Update the geometric field parameters (fluid)
for NodeIndex in range(1,NumberOfFluidNodes+1):
    NodeDomain = FluidDecomposition.NodeDomainGet(FluidNodeNumbers[NodeIndex-1],Mesh2ComponentNumberSpace)
    if (NodeDomain == ComputationalNodeNumber):
        GeometricField2.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                 versionIdx,derivIdx,FluidNodeNumbers[NodeIndex-1],
                                                 Components[0],FluidGeometryY[NodeIndex-1])
        GeometricField2.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                 versionIdx,derivIdx,FluidNodeNumbers[NodeIndex-1],
                                                 Components[1],FluidGeometryZ[NodeIndex-1])
        if (NumberOfDimensions == 3):
            GeometricField2.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                     versionIdx,derivIdx,FluidNodeNumbers[NodeIndex],
                                                     Components[2],FluidGeometryX[NodeIndex-1])
# Update the geometric field parameters (interface)
for NodeIndex in range(1,NumberOfInterfaceNodes+1):
    NodeDomain = InterfaceDecomposition.NodeDomainGet(int(InterfaceNodeNumbersForGeometry[NodeIndex-1]),
                                                      InterfaceMeshComponentNumber)
    if (NodeDomain == ComputationalNodeNumber):
        InterfaceGeometricField1.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                          versionIdx,derivIdx,int(InterfaceNodeNumbersForGeometry[NodeIndex-1]),
                                                          Components[0],InterfaceGeometryY[NodeIndex-1])
        InterfaceGeometricField1.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                          versionIdx,derivIdx,int(InterfaceNodeNumbersForGeometry[NodeIndex-1]),
                                                          Components[1],InterfaceGeometryZ[NodeIndex-1])
        if (NumberOfDimensions == 3):
            InterfaceGeometricField1.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                              versionIdx,derivIdx,int(InterfaceNodeNumbersForGeometry[1][NodeIndex-1]),
                                                              Components[2],InterfaceGeometryX[NodeIndex-1])

GeometricField1.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
GeometricField1.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
GeometricField2.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
GeometricField2.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
InterfaceGeometricField1.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
InterfaceGeometricField1.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

#================================================================================================================================
#  Equations Set
#================================================================================================================================

# Create the equations set for the solid region - Finite Elasticity Mooney-Rivlin
if (ExampleFileProgressDiagnostics):
    print " == >> SOLID EQUATION SET << == "
EquationsSetField1 = iron.Field()
SolidEquationsSet = iron.EquationsSet()
SolidEquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                  iron.EquationsSetTypes.FINITE_ELASTICITY,
                                  iron.EquationsSetSubtypes.MOONEY_RIVLIN]
SolidEquationsSet.CreateStart(SolidEquationsSetUserNumber,Region1,GeometricField1,
                              SolidEquationsSetSpecification,SolidEquationsSetFieldUserNumber,
                              EquationsSetField1)
SolidEquationsSet.CreateFinish()

# Create the equations set for the fluid region - ALE Navier-Stokes
if (ExampleFileProgressDiagnostics):
    print " == >> FLUID EQUATION SET << == "
EquationsSetField2 = iron.Field()
FluidEquationsSet = iron.EquationsSet()
FluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                  iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                  iron.EquationsSetSubtypes.ALE_NAVIER_STOKES]
FluidEquationsSet.CreateStart(FluidEquationsSetUserNumber,Region2,GeometricField2,
                              FluidEquationsSetSpecification,FluidEquationsSetFieldUserNumber,
                              EquationsSetField2)
FluidEquationsSet.CreateFinish()

# Create the equations set for the moving mesh
if (ExampleFileProgressDiagnostics):
    print " == >> MOVING MESH EQUATION SET << == "
EquationsSetFieldMovingMesh = iron.Field()
MovingMeshEquationsSet = iron.EquationsSet()
MovingMeshEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.LAPLACE_EQUATION,
                                       iron.EquationsSetSubtypes.MOVING_MESH_LAPLACE]
MovingMeshEquationsSet.CreateStart(MovingMeshEquationsSetUserNumber,Region2,GeometricField2,
                                   MovingMeshEquationsSetSpecification,EquationsSetFieldMovingMeshUserNumber,
                                   EquationsSetFieldMovingMesh)
MovingMeshEquationsSet.CreateFinish()

#================================================================================================================================
#  Dependent Field
#================================================================================================================================

# Create the equations set dependent field variables for the first equations set
if (ExampleFileProgressDiagnostics):
    print " == >> SOLID DEPENDENT FIELD << == "
# Create the dependent field
DependentField1 = iron.Field()
SolidEquationsSet.DependentCreateStart(SolidDependentFieldUserNumber,DependentField1)
DependentField1.VariableLabelSet(iron.FieldVariableTypes.U,'SolidDF')
for component_idx in range(1,NumberOfDimensions+1):
    DependentField1.ComponentMeshComponentSet(iron.FieldVariableTypes.U,component_idx,Mesh1ComponentNumberSpace)
    DependentField1.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,component_idx,Mesh1ComponentNumberSpace)
DependentField1.ComponentMeshComponentSet(iron.FieldVariableTypes.U,NumberOfDimensions+1,PressureMeshComponent)
DependentField1.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,NumberOfDimensions+1,PressureMeshComponent)
if (PressureMeshComponent == 1):
    DependentField1.ComponentInterpolationSet(iron.FieldVariableTypes.U,NumberOfDimensions+1,
     iron.FieldInterpolationTypes.ELEMENT_BASED)
    DependentField1.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,NumberOfDimensions+1,
     iron.FieldInterpolationTypes.ELEMENT_BASED)
else:
    DependentField1.ComponentInterpolationSet(iron.FieldVariableTypes.U,NumberOfDimensions+1,
     iron.FieldInterpolationTypes.NODE_BASED)
    DependentField1.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,NumberOfDimensions+1,
     iron.FieldInterpolationTypes.NODE_BASED)
DependentField1.ScalingTypeSet(iron.FieldScalingTypes.NONE)
SolidEquationsSet.DependentCreateFinish()

# Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
GeometricField1.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,
 iron.FieldParameterSetTypes.VALUES,1,DependentField1,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
GeometricField1.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,
 iron.FieldParameterSetTypes.VALUES,2,DependentField1,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
if (NumberOfDimensions == 3):
    GeometricField1.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,
     iron.FieldParameterSetTypes.VALUES,3,DependentField1,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
     
DependentField1.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,NumberOfDimensions+1,-MooneyRivlin1)
DependentField1.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
DependentField1.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# Create the equations set dependent field variables for dynamic Navier-Stokes
if (ExampleFileProgressDiagnostics):
    print " == >> FLUID DEPENDENT FIELD << == "
DependentField2 = iron.Field()
FluidEquationsSet.DependentCreateStart(FluidDependentFieldUserNumber,DependentField2)
DependentField2.VariableLabelSet(iron.FieldVariableTypes.U,'FluidDF')
# Set the mesh component to be used by the field components.
for ComponentNumber in range(1,NumberOfDimensions+1):
    DependentField2.ComponentMeshComponentSet(iron.FieldVariableTypes.U,ComponentNumber,Mesh2ComponentNumberSpace)
    DependentField2.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,ComponentNumber,Mesh2ComponentNumberSpace)
DependentField2.ComponentMeshComponentSet(iron.FieldVariableTypes.U,NumberOfDimensions+1,Mesh2ComponentNumberPressure)
DependentField2.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,NumberOfDimensions+1,Mesh2ComponentNumberPressure)
# Finish the equations set dependent field variables
FluidEquationsSet.DependentCreateFinish()
# Initialise dependent field
for ComponentNumber in range(1,NumberOfDimensions+1):
    DependentField2.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                ComponentNumber,InitialFieldNavierStokes[ComponentNumber-1])
# Initialise pressure component
DependentField2.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                            NumberOfDimensions+1,0.1)
     
# Create the equations set dependent field variables for moving mesh
if (ExampleFileProgressDiagnostics):
    print " == >> MOVING MESH DEPENDENT FIELD << == "
DependentFieldMovingMesh = iron.Field()
MovingMeshEquationsSet.DependentCreateStart(DependentFieldMovingMeshUserNumber,DependentFieldMovingMesh)
DependentFieldMovingMesh.VariableLabelSet(iron.FieldVariableTypes.U,'MovingMeshDF')
# Set the mesh component to be used by the field components.
for ComponentNumber in range(1,NumberOfDimensions+1):
    DependentFieldMovingMesh.ComponentMeshComponentSet(iron.FieldVariableTypes.U,ComponentNumber,Mesh2ComponentNumberSpace)
    DependentFieldMovingMesh.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,ComponentNumber,Mesh2ComponentNumberSpace)
# Finish the equations set dependent field variables
MovingMeshEquationsSet.DependentCreateFinish()
# Initialise dependent field moving mesh
for ComponentNumber in range(1,NumberOfDimensions+1):
    DependentFieldMovingMesh.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                         ComponentNumber,InitialFieldMovingMesh[ComponentNumber-1])
     
#================================================================================================================================
#  Materials Field
#================================================================================================================================

# Create the material field
if (ExampleFileProgressDiagnostics):
    print " == >> SOLID MATERIAL FIELD << == "
MaterialField1 = iron.Field()
SolidEquationsSet.MaterialsCreateStart(SolidMaterialFieldUserNumber,MaterialField1)
MaterialField1.VariableLabelSet(iron.FieldVariableTypes.U,'Material1')
MaterialField1.VariableLabelSet(iron.FieldVariableTypes.V,'SolidDensity')
SolidEquationsSet.MaterialsCreateFinish()
# Set Mooney-Rivlin constants c10 and c01 (default?2.0 and 6.0) respectively
MaterialField1.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,MooneyRivlin1)
MaterialField1.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,MooneyRivlin2)
MaterialField1.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,1,SolidDensity)

# Create the equations set materials field variables for dynamic Navier-Stokes
if (ExampleFileProgressDiagnostics):
    print " == >> FLUID MATERIAL FIELD << == "
MaterialField2 = iron.Field()
FluidEquationsSet.MaterialsCreateStart(FluidMaterialFieldUserNumber,MaterialField2)
# Finish the equations set materials field variables
FluidEquationsSet.MaterialsCreateFinish()
MaterialField2.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
 FluidMaterialFieldComponentMu,FluidDynamicViscosity)
MaterialField2.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
 FluidMaterialFieldComponentRho,FluidDensity)

# Create the equations set materials field variables for moving mesh
if (ExampleFileProgressDiagnostics):
    print " == >> MATERIAL FIELD MOVING MESH << == "
MaterialFieldMovingMesh = iron.Field()
MovingMeshEquationsSet.MaterialsCreateStart(MaterialFieldMovingMeshUserNumber,MaterialFieldMovingMesh)
# Finish the equations set materials field variables
MovingMeshEquationsSet.MaterialsCreateFinish()
MaterialFieldMovingMesh.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
 MaterialFieldMovingMeshUserNumberK,MovingMeshParameterK)
    
#================================================================================================================================
#  Source Field
#================================================================================================================================

if (GravityFlag):
    if (ExampleFileProgressDiagnostics):
        print " == >> SOURCE FIELD - GRAVITY << == "
    #Create the source field with the gravity vector
    SourceField1 = iron.Field()
    SolidEquationsSet.SourceCreateStart(SourceFieldUserNumber,SourceField1)
    SourceField1.ScalingTypeSet(iron.FieldScalingTypes.NONE)
    SolidEquationsSet.SourceCreateFinish()
    for component_idx in range(1,NumberOfDimensions+1):
        SourceField1.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
         component_idx,Gravity[component_idx-1])
         
#================================================================================================================================
# Independent Field
#================================================================================================================================

# Mesh velocity for fluid equations set and mesh stiffness for moving mesh equations set
if (ExampleFileProgressDiagnostics):
    print " == >> FLUID INDEPENDENT FIELD << == "
# Create the equations set independent field variables for ALE Navier-Stokes
IndependentField2 = iron.Field()
FluidEquationsSet.IndependentCreateStart(IndependentField2UserNumber,IndependentField2)
IndependentField2.VariableLabelSet(iron.FieldVariableTypes.U,'FluidInDF')
# Set the mesh component to be used by the field components.
for ComponentNumber in range(1,NumberOfDimensions+1):
    IndependentField2.ComponentMeshComponentSet(iron.FieldVariableTypes.U,ComponentNumber,Mesh2ComponentNumberSpace)
# Finish the equations set independent field variables
FluidEquationsSet.IndependentCreateFinish()
  
if (ExampleFileProgressDiagnostics):
    print " == >> INDEPENDENT FIELD MOVING MESH << == "
# Create the equations set independent field variables for moving mesh
IndependentFieldMovingMesh = iron.Field()
MovingMeshEquationsSet.IndependentCreateStart(IndependentFieldMovingMeshUserNumber,IndependentFieldMovingMesh)
IndependentFieldMovingMesh.VariableLabelSet(iron.FieldVariableTypes.U,'MovingMeshInDF')
# Set the mesh component to be used by the field components.
for ComponentNumber in range(1,NumberOfDimensions+1):
    IndependentFieldMovingMesh.ComponentMeshComponentSet(iron.FieldVariableTypes.U,ComponentNumber,Mesh2ComponentNumberSpace)    
# Finish the equations set independent field variables
MovingMeshEquationsSet.IndependentCreateFinish()
# Initialise independent field moving mesh
IndependentFieldMovingMesh.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
 iron.FieldParameterSetTypes.VALUES,IndependentFieldMovingMeshUserNumberK,MovingMeshParameterK)
 
#================================================================================================================================
#  Equations
#================================================================================================================================

# 1st Equations Set
if (ExampleFileProgressDiagnostics):
    print " == >> SOLID EQUATIONS << == "
Equations1 = iron.Equations()
SolidEquationsSet.EquationsCreateStart(Equations1)
Equations1.sparsityType = iron.EquationsSparsityTypes.SPARSE
# (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
Equations1.outputType = iron.EquationsOutputTypes.NONE
SolidEquationsSet.EquationsCreateFinish()

# 2nd Equations Set
if (ExampleFileProgressDiagnostics):
    print " == >> FLUID EQUATIONS << == "
Equations2 = iron.Equations()
FluidEquationsSet.EquationsCreateStart(Equations2)
Equations2.sparsityType = iron.EquationsSparsityTypes.SPARSE
# (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
Equations2.outputType = iron.EquationsOutputTypes.NONE
FluidEquationsSet.EquationsCreateFinish()

# 3rd Equations Set
if (ExampleFileProgressDiagnostics):
    print " == >> MOVING MESH EQUATIONS << == "
EquationsMovingMesh = iron.Equations()
MovingMeshEquationsSet.EquationsCreateStart(EquationsMovingMesh)
EquationsMovingMesh.sparsityType = iron.EquationsSparsityTypes.SPARSE
# (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
EquationsMovingMesh.outputType = iron.EquationsOutputTypes.NONE
MovingMeshEquationsSet.EquationsCreateFinish()

#================================================================================================================================
#  Interface Condition
#================================================================================================================================

# Create an interface condition between the two meshes
if (ExampleFileProgressDiagnostics):
    print " == >> INTERFACE CONDITIONS << == "
InterfaceCondition = iron.InterfaceCondition()
InterfaceCondition.CreateStart(InterfaceConditionUserNumber,Interface1,InterfaceGeometricField1)
# Specify the method for the interface condition
InterfaceCondition.MethodSet(iron.InterfaceConditionMethods.LAGRANGE_MULTIPLIERS)
# Specify the type of interface condition operator
InterfaceCondition.OperatorSet(iron.InterfaceConditionOperators.SOLID_FLUID)
# Add in the dependent variables from the equations sets
InterfaceCondition.DependentVariableAdd(SolidMeshIndex,SolidEquationsSet,iron.FieldVariableTypes.U)
InterfaceCondition.DependentVariableAdd(FluidMeshIndex,FluidEquationsSet,iron.FieldVariableTypes.U)
# Finish creating the interface condition
InterfaceCondition.CreateFinish()

# Create the Lagrange multipliers field
if (ExampleFileProgressDiagnostics):
    print " == >> INTERFACE LAGRANGE FIELD << == "
LagrangeField1 = iron.Field()
InterfaceCondition.LagrangeFieldCreateStart(LagrangeFieldUserNumber,LagrangeField1)
LagrangeField1.VariableLabelSet(iron.FieldVariableTypes.U,"InterfaceLF")
# Finish the Lagrange multipliers field
InterfaceCondition.LagrangeFieldCreateFinish()
for ComponentNumber in range(1,NumberOfDimensions+1):
    LagrangeField1.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,ComponentNumber,0.0)
LagrangeField1.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
LagrangeField1.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# Create the interface condition equations
if (ExampleFileProgressDiagnostics):
    print " == >> INTERFACE EQUATIONS << == "
InterfaceEquations = iron.InterfaceEquations()
InterfaceCondition.EquationsCreateStart(InterfaceEquations)
# Set the interface equations sparsity
InterfaceEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
# Set the interface equations output
InterfaceEquations.outputType = iron.EquationsOutputTypes.NONE
# Finish creating the interface equations
InterfaceCondition.EquationsCreateFinish()

#================================================================================================================================
#  Problem
#================================================================================================================================

# Create a problem
if (ExampleFileProgressDiagnostics):
    print " == >> COUPLED PROBLEM << == "
CoupledProblem = iron.Problem()
ProblemSpecification = [iron.ProblemClasses.MULTI_PHYSICS,
                        iron.ProblemTypes.FINITE_ELASTICITY_NAVIER_STOKES,
                        iron.ProblemSubtypes.FINITE_ELASTICITY_NAVIER_STOKES_ALE]
CoupledProblem.CreateStart(CoupledProblemUserNumber,ProblemSpecification)
CoupledProblem.CreateFinish()

#================================================================================================================================
#  Control Loop
#================================================================================================================================

# Create the problem control loop
if (ExampleFileProgressDiagnostics):
    print " == >> PROBLEM CONTROL LOOP << == "
ControlLoop = iron.ControlLoop()
CoupledProblem.ControlLoopCreateStart()
CoupledProblem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],ControlLoop)
ControlLoop.LabelSet('TimeLoop')
ControlLoop.TimesSet(StartTime,StopTime,TimeStepSize)
ControlLoop.TimeInputSet(MaterialSpecification)
ControlLoop.TimeOutputSet(OutputFrequency)
CoupledProblem.ControlLoopCreateFinish()

#================================================================================================================================
#  Solvers
#================================================================================================================================

# Create the problem solver
if (ExampleFileProgressDiagnostics):
    print " == >> PROBLEM SOLVERS << == "
LinearSolverMovingMesh = iron.Solver()
DynamicSolver = iron.Solver()
NonlinearSolver = iron.Solver()
LinearSolver = iron.Solver()

CoupledProblem.SolversCreateStart()
# Linear solver for moving mesh
CoupledProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],LinearSolverMovingMeshIndex,LinearSolverMovingMesh)
LinearSolverMovingMesh.OutputTypeSet(LinearSolverMovingMesh_OutputType)
# Solvers for coupled FiniteElasticity NavierStokes problem
# Get the dynamic ALE solver
CoupledProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],DynamicSolverIndex,DynamicSolver)
DynamicSolver.OutputTypeSet(DynamicSolver_OutputType)
DynamicSolver.DynamicThetaSet(DynamicSolver_Theta)
# Get the dynamic nonlinear solver
DynamicSolver.DynamicNonlinearSolverGet(NonlinearSolver)
NonlinearSolver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.LINEAR)
NonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
NonlinearSolver.NewtonMaximumFunctionEvaluationsSet(MaxFunctionEvaluations)
NonlinearSolver.OutputTypeSet(NonlinearSolver_OutputType)
NonlinearSolver.NewtonAbsoluteToleranceSet(AbsoluteTolerance)
NonlinearSolver.NewtonMaximumIterationsSet(MaximumIterations)
NonlinearSolver.NewtonRelativeToleranceSet(RelativeTolerance)
NonlinearSolver.NewtonLineSearchAlphaSet(LinesearchAlpha)
# Get the dynamic nonlinear linear solver
NonlinearSolver.NewtonLinearSolverGet(LinearSolver)
# Choose type of linear solver
if (False):
    LinearSolver.LinearTypeSet(iron.LinearSolverTypes.ITERATIVE)
    LinearSolver.LinearIterativeMaximumIterationsSet(MaximumIterations)
    LinearSolver.LinearIterativeDivergenceToleranceSet(DivergenceTolerance)
    LinearSolver.LinearIterativeRelativeToleranceSet(RelativeTolerance)
    LinearSolver.LinearIterativeAbsoluteToleranceSet(AbsoluteTolerance)
else:
    LinearSolver.LinearTypeSet(iron.LinearSolverTypes.ITERATIVE)
LinearSolver.OutputTypeSet(LinearSolver_OutputType)
# Finish the creation of the problem solver
CoupledProblem.SolversCreateFinish()

#================================================================================================================================
#  Solver Equations
#================================================================================================================================

# Create the problem solver equations
if (ExampleFileProgressDiagnostics):
    print " == >> MOVING MESH SOLVER EQUATIONS << == "
LinearSolverMovingMesh = iron.Solver()
LinearSolverMovingMeshEquations = iron.SolverEquations()
CoupledProblem.SolverEquationsCreateStart()
# Get the linear solver equations
CoupledProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],LinearSolverMovingMeshIndex,LinearSolverMovingMesh)
LinearSolverMovingMesh.SolverEquationsGet(LinearSolverMovingMeshEquations)
LinearSolverMovingMeshEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
# Add in the equations set
MovingMeshEquationsSet = LinearSolverMovingMeshEquations.EquationsSetAdd(MovingMeshEquationsSet)

# Get the dynamic solver equations
if (ExampleFileProgressDiagnostics):
    print " == >> SOLVER EQUATIONS << == "
DynamicSolver = iron.Solver()
CoupledSolverEquations = iron.SolverEquations()
CoupledProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],DynamicSolverIndex,DynamicSolver)
DynamicSolver.SolverEquationsGet(CoupledSolverEquations)
CoupledSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
CoupledSolverEquationsSet = CoupledSolverEquations.EquationsSetAdd(SolidEquationsSet)
CoupledSolverEquationsSet = CoupledSolverEquations.EquationsSetAdd(FluidEquationsSet)
CoupledSolverEquationsSet = CoupledSolverEquations.InterfaceConditionAdd(InterfaceCondition)
# Set the time dependence of the interface matrix to determine the interface matrix coefficient in the solver matrix
# (basiy position in big coupled matrix system)
InterfaceMatrices.TimeDependenceTypeSet(SolidEquationsSetIndex,True,[CMISS_INTERFACE_MATRIX_STATIC,CMISS_INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC])
InterfaceMatrices.TimeDependenceTypeSet(FluidEquationsSetIndex,True,[CMISS_INTERFACE_MATRIX_STATIC,CMISS_INTERFACE_MATRIX_STATIC])
# Finish the creation of the problem solver equations
CoupledProblem.SolverEquationsCreateFinish()
    
#================================================================================================================================
#  Boundary Conditions
#================================================================================================================================

# Start the creation of the equations set boundary conditions
if (ExampleFileProgressDiagnostics):
    print " == >> BOUNDARY CONDITIONS << == "
BoundaryConditions = iron.BoundaryConditions()
CoupledSolverEquations.BoundaryConditionsCreateStart(BoundaryConditions)
# No displacement boundary for solid
for S in range(1,len(NoDisplacementNodes)+1):
    NodeNumber = NoDisplacementNodes[S-1]
    SolidDecomposition.NodeDomainGet(NodeNumber,1)
    if (NodeDomain == ComputationalNodeNumber):
        BoundaryConditions.AddNode(DependentField1,iron.FieldVariableTypes.U,1,1,NodeNumber,1,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        BoundaryConditions.AddNode(DependentField1,iron.FieldVariableTypes.U,1,1,NodeNumber,2,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        if (NumberOfDimensions == 3):
            BoundaryConditions.AddNode(DependentField1,iron.FieldVariableTypes.U,1,1,NodeNumber,3,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
# Set outlet (zero) pressure nodes
for S in range(1,len(OutletNodes)+1):
    NodeNumber = int(OutletNodes[S-1])
    FluidDecomposition.NodeDomainGet(NodeNumber,1)
    if (NodeDomain == ComputationalNodeNumber):
        BoundaryConditions.SetNode(DependentField2,iron.FieldVariableTypes.U,1,1,NodeNumber,NumberOfDimensions+1,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
# Inlet velocity nodes, small starting velocity in 1st coordinate direction
for S in range(1,len(InletNodes)+1):
    NodeNumber = int(InletNodes[S-1])
    FluidDecomposition.NodeDomainGet(NodeNumber,1)
    if (NodeDomain == ComputationalNodeNumber):
        BoundaryConditions.SetNode(DependentField2,iron.FieldVariableTypes.U,1,1,NodeNumber,1,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        BoundaryConditions.SetNode(DependentField2,iron.FieldVariableTypes.U,1,1,NodeNumber,2,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        if (NumberOfDimensions == 3):
            BoundaryConditions.SetNode(DependentField2,iron.FieldVariableTypes.U,1,1,NodeNumber,3,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
# Set no-slip BC
for S in range(1,len(NoSlipNodes)+1):
    NodeNumber = int(NoSlipNodes[S-1])
    FluidDecomposition.NodeDomainGet(NodeNumber,1)
    if (NodeDomain == ComputationalNodeNumber):
        BoundaryConditions.SetNode(DependentField2,iron.FieldVariableTypes.U,1,1,NodeNumber,1,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        BoundaryConditions.SetNode(DependentField2,iron.FieldVariableTypes.U,1,1,NodeNumber,2,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        if (NumberOfDimensions == 3):
            BoundaryConditions.SetNode(DependentField2,iron.FieldVariableTypes.U,1,1,NodeNumber,3,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
# Set slip BC
for S in range(1,len(SlipNodesTop)+1):
    NodeNumber = int(SlipNodesTop[S-1])
    FluidDecomposition.NodeDomainGet(NodeNumber,1)
    if (NodeDomain == ComputationalNodeNumber):
        BoundaryConditions.SetNode(DependentField2,iron.FieldVariableTypes.U,1,1,NodeNumber,NumberOfDimensions,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
# Set slip BC
if (NumberOfDimensions == 3):
    for S in range(1,len(SlipNodesRightLeft)+1):
        NodeNumber = int(SlipNodesRightLeft[S-1])
        FluidDecomposition.NodeDomainGet(NodeNumber,1)
        if (NodeDomain == ComputationalNodeNumber):
            BoundaryConditions.SetNode(DependentField2,iron.FieldVariableTypes.U,1,1,NodeNumber,1,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
if (NumberOfDimensions == 2):
    # Remove dof's at nodes where solid displacement and zero velocity is set (first n last interface node)
    BoundaryConditions.SetNode(LagrangeField1,iron.FieldVariableTypes.U,1,1,1,1,
                               iron.BoundaryConditionsTypes.FIXED,0.0)
    BoundaryConditions.SetNode(LagrangeField1,iron.FieldVariableTypes.U,1,1,1,2,
                               iron.BoundaryConditionsTypes.FIXED,0.0)
    BoundaryConditions.SetNode(LagrangeField1,iron.FieldVariableTypes.U,1,1,NumberOfInterfaceNodes,1,
                               iron.BoundaryConditionsTypes.FIXED,0.0)
    BoundaryConditions.SetNode(LagrangeField1,iron.FieldVariableTypes.U,1,1,NumberOfInterfaceNodes,2,
                               iron.BoundaryConditionsTypes.FIXED,0.0)
else:
    if (CheckWithoutInterfaceCondition):
        for S in range(1,len(InterfaceNodeNumbersForGeometry)+1):
            NodeNumber = int(InterfaceNodeNumbersForGeometry[S-1])
            InterfaceDecomposition.NodeDomainGet(NodeNumber,1)
            if (NodeDomain == ComputationalNodeNumber):
                BoundaryConditions.SetNode(LagrangeField1,iron.FieldVariableTypes.U,1,1,NodeNumber,1,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                BoundaryConditions.SetNode(LagrangeField1,iron.FieldVariableTypes.U,1,1,NodeNumber,2,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                BoundaryConditions.SetNode(LagrangeField1,iron.FieldVariableTypes.U,1,1,NodeNumber,3,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
    else:
        for S in range(1,len(LagrangeNodes)+1):
            NodeNumber = int(LagrangeNodes[S-1])
            InterfaceDecomposition.NodeDomainGet(NodeNumber,1)
            if (NodeDomain == ComputationalNodeNumber):
                BoundaryConditions.SetNode(LagrangeField1,iron.FieldVariableTypes.U,1,1,NodeNumber,1,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                BoundaryConditions.SetNode(LagrangeField1,iron.FieldVariableTypes.U,1,1,NodeNumber,2,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                BoundaryConditions.SetNode(LagrangeField1,iron.FieldVariableTypes.U,1,1,NodeNumber,3,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
# Finish equations set boundary conditions
CoupledSolverEquations.BoundaryConditionsCreateFinish()
  
# Start the creation of the moving mesh boundary conditions
BoundaryConditionsMovingMesh = iron.BoundaryConditions()
LinearSolverMovingMeshEquations.BoundaryConditionsCreateStart(BoundaryConditionsMovingMesh)
# Fixed boundary nodes. May be used to move nodes..
for S in range(1,len(MovedYNodes)+1):
    NodeNumber = int(MovedYNodes[S-1])
    FluidDecomposition.NodeDomainGet(NodeNumber,1)
    if (NodeDomain == ComputationalNodeNumber):
        BoundaryConditionsMovingMesh.SetNode(DependentFieldMovingMesh,iron.FieldVariableTypes.U,1,1,NodeNumber,1,
                                             iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
        BoundaryConditionsMovingMesh.SetNode(DependentFieldMovingMesh,iron.FieldVariableTypes.U,1,1,NodeNumber,2,
                                             iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
        if (NumberOfDimensions == 3):
            BoundaryConditionsMovingMesh.SetNode(DependentFieldMovingMesh,iron.FieldVariableTypes.U,1,1,NodeNumber,3,
                                                 iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
# Mesh nodes that are moving wall nodes
for S in range(1,len(MovedNodes)+1):
    NodeNumber = int(MovedNodes[S-1])
    FluidDecomposition.NodeDomainGet(NodeNumber,1)
    if (NodeDomain == ComputationalNodeNumber):
        BoundaryConditionsMovingMesh.SetNode(DependentFieldMovingMesh,iron.FieldVariableTypes.U,1,1,NodeNumber,1,
                                             iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
        BoundaryConditionsMovingMesh.SetNode(DependentFieldMovingMesh,iron.FieldVariableTypes.U,1,1,NodeNumber,2,
                                             iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
        if (NumberOfDimensions == 3):
            BoundaryConditionsMovingMesh.SetNode(DependentFieldMovingMesh,iron.FieldVariableTypes.U,1,1,NodeNumber,3,
                                                 iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
# Mesh nodes that are fixed in space
for S in range(1,len(FixedNodes)+1):
    NodeNumber = int(FixedNodes[S-1])
    FluidDecomposition.NodeDomainGet(NodeNumber,1)
    if (NodeDomain == ComputationalNodeNumber):
        BoundaryConditionsMovingMesh.SetNode(DependentFieldMovingMesh,iron.FieldVariableTypes.U,1,1,NodeNumber,1,
                                             iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
        BoundaryConditionsMovingMesh.SetNode(DependentFieldMovingMesh,iron.FieldVariableTypes.U,1,1,NodeNumber,2,
                                             iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
        if (NumberOfDimensions == 3):
            BoundaryConditionsMovingMesh.SetNode(DependentFieldMovingMesh,iron.FieldVariableTypes.U,1,1,NodeNumber,3,
                                                 iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
# Finish moving mesh boundary conditions
LinearSolverMovingMeshEquations.BoundaryConditionsCreateFinish()

#================================================================================================================================
#  Run Solvers
#================================================================================================================================

# Solve the problem
print "Solving problem..."
start = time.time()
CoupledProblem.Solve()
end = time.time()
elapsed = end - start
print "Calculation Time = %3.4f" %elapsed
print "Problem solved!"
print "#"

#================================================================================================================================
#  Finish Program
#================================================================================================================================
