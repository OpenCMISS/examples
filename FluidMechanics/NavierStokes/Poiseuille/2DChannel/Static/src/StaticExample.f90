!> \file
!> \author David Ladd
!> \brief This is an example program to numerically solve a 2D Poiseuille 
!>  flow problem with the Navier-Stokes equations and validate against the
!>  exact analytic form using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example FluidMechanics/Poiseuille/Static/src/StaticExample.f90
!! Example program to solve a static Poiseuille equation using openCMISS calls.
!<

!> Main program
PROGRAM STATICPOISEUILLEEXAMPLE

  USE OPENCMISS
  USE FLUID_MECHANICS_IO_ROUTINES
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField

  !Test program parameters


!==============================================================

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=10.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: INLET_VELOCITY=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: OUTLET_PRESSURE=0.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: VISCOSITY=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: DENSITY=1.0_CMISSDP

!==============================================================


  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=14

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_DIMENSIONS,INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS
  INTEGER(CMISSIntg) :: component_idx
  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  REAL(CMISSDP) :: POSITION,TOTAL_LENGTH
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT

  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber,FirstNodeDomain,LastNodeDomain

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_SPACE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_PRESSURE
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_SPACE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_PRESSURE
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
!   INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_NAVIER_STOKES
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_NODES_NAVIER_STOKES
  INTEGER(CMISSIntg) :: NUMBER_OF_OUTLET_NODES_NAVIER_STOKES

  INTEGER(CMISSIntg) :: EQUATIONS_NAVIER_STOKES_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_COUNTER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_NAVIER_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_NODES_NAVIER_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: OUTLET_NODES_NAVIER_STOKES

  REAL(CMISSDP) :: INITIAL_FIELD_NAVIER_STOKES(2)
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_NAVIER_STOKES(2)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE
  REAL(CMISSDP) :: MU_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: RHO_PARAM_NAVIER_STOKES

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG
  LOGICAL :: FIXED_WALL_NODES_NAVIER_STOKES_FLAG
  LOGICAL :: LID_NODES_NAVIER_STOKES_FLAG


  !CMISS variables


  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSBasisType) :: BasisSpace,BasisVelocity,BasisPressure

  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField,AnalyticField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSProblemType) :: Problem

  TYPE(CMISSSolverType) :: LinearSolver,NonlinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSFieldType) :: EquationsSetField

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables

  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,BoundaryNodeDomain
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: Err

#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

!   !Get input arguments
!   NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
!   IF(NUMBER_OF_ARGUMENTS >= 2) THEN
!     !If we have enough arguments then use the first four for setting up the problem. The subsequent arguments may be used to
!     !pass flags to, say, PETSc.
!     CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
!     IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
!     READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_X_ELEMENTS
!     IF(NUMBER_GLOBAL_X_ELEMENTS<=0) CALL HANDLE_ERROR("Invalid number of X elements.")
!     CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
!     IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 2.")
!     READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) INTERPOLATION_TYPE
!     IF(INTERPOLATION_TYPE<=0) CALL HANDLE_ERROR("Invalid Interpolation specification.")
!   ELSE
!     !If there are not enough arguments default the problem specification
!     NUMBER_GLOBAL_X_ELEMENTS=2
!     INTERPOLATION_TYPE=1
!   ENDIF

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NUMBER_GLOBAL_X_ELEMENTS=10
  NUMBER_GLOBAL_Y_ELEMENTS=10
  NUMBER_OF_DOMAINS=NumberOfComputationalNodes

  !Set all diganostic levels on for testing

  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,2,Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !Start the creation of a basis (bilinear lagrange)
  CALL CMISSBasisTypeInitialise(Basis,Err)
  CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
  CALL CMISSBasisNumberOfXiSet(Basis,2,Err)
  CALL CMISSBasisCreateFinish(Basis,Err)

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)
  !Define the mesh on the region
  CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/LENGTH,HEIGHT/),Err)
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/),Err)
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,1,Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

  !Create the equations_set
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetFluidmechanicsClass, &
    & CMISSEquationsSetPoiseuilleEquationType,CMISSEquationsSetStaticPoiseuilleSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field variables
  CALL CMISSFieldTypeInitialise(MaterialsField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)
  PIPE_LENGTH=SQRT(DOT_PRODUCT(POSITION,POSITION))
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,PIPE_LENGTH,Err)

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemFluidmechanicsClass,CMISSProblemPoiseuilleEquationType, &
    & CMISSProblemStaticPoiseuilleSubtype,Err)
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
  !Set solver parameters
  !CALL CMISSSolverLinearTypeSet(Solver,CMISSSolverLinearIterativeSolveType,Err)
  CALL CMISSSolverLinearTypeSet(Solver,CMISSSolverLinearDirectSolveType,Err)
  !CALL CMISSSolverLinearIterativeRelativeToleranceSet(Solver,1.0E-8_CMISSDP,Err)
  !CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(Solver,1.0E-8_CMISSDP,Err)
  !CALL CMISSSolverLinearIterativeMaximumIterationsSet(Solver,10000,Err)
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Set up the boundary conditions
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the fixed boundary conditions at the first node and last nodes
  FirstNodeNumber=1
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSRegionNodesGet(Region,Nodes,Err)
  CALL CMISSNodesNumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL CMISSDecompositionNodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL CMISSDecompositionNodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,FirstNodeNumber,1, &
      & CMISSBoundaryConditionFixed,100.0_CMISSDP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,LastNodeNumber,1, &
      & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  ENDIF
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)

  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"Poiseuille","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"Poiseuille","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR

END PROGRAM STATICPOISEUILLEEXAMPLE
