!> \file
!> \author David Ladd
!> \brief This is an example program to solve a transient viscous burgers equation using OpenCMISS calls.
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

!> \example FluidMechanics/Burgers/Burgers/src/BurgersExample.f90
!! Example program to solve a burgers equation using openCMISS calls.
!! \htmlinclude FluidMechanics/Burgers/Burgers/history.html
!<

!> Main program
PROGRAM BURGERSEXAMPLE


  USE OPENCMISS
  USE MPI


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField

  !Test program parameters
  
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopNode=0
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: SolverUserNumber=1
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER

  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
  
  INTEGER(CMISSIntg) :: MPI_IERROR

  INTEGER(CMISSIntg) :: NONLINEAR_SOLVER_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: EQUATIONS_OUTPUT

  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA

  LOGICAL :: LINEAR_SOLVER_DIRECT_FLAG

  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_OUTPUT_FREQUENCY
  REAL(CMISSDP) :: DYNAMIC_SOLVER_START_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_STOP_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_THETA
  REAL(CMISSDP) :: DYNAMIC_SOLVER_TIME_INCREMENT

  REAL(CMISSDP) :: NU_PARAM
  REAL(CMISSDP) :: LENGTH

  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_COUNTER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber,FirstNodeDomain,LastNodeDomain
  
  !Program types

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField,AnalyticField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: DynamicSolver,NonlinearSolver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

  LOGICAL :: EXPORT_FIELD

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

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------
  
  ! Set number of elements for FEM discretization
  NUMBER_GLOBAL_X_ELEMENTS=4
  NUMBER_OF_DOMAINS=NumberOfComputationalNodes

  ! Set viscous coefficient
  NU_PARAM = 1.0_CMISSDP
  ! Set length of domain
  LENGTH = 1.0_CMISSDP

  !Set output parameters
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  DYNAMIC_SOLVER_OUTPUT_TYPE=CMISSSolverNoOutput
  NONLINEAR_SOLVER_OUTPUT_TYPE=CMISSSolverNoOutput
  LINEAR_SOLVER_OUTPUT_TYPE=CMISSSolverNoOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_OUTPUT=CMISSEquationsNoOutput

  !Set time parameter
  DYNAMIC_SOLVER_START_TIME=0.0_CMISSDP
  DYNAMIC_SOLVER_STOP_TIME=5.0_CMISSDP 
  DYNAMIC_SOLVER_TIME_INCREMENT=1.0_CMISSDP
  DYNAMIC_SOLVER_THETA=1.0_CMISSDP

  !Set solver parameters
  LINEAR_SOLVER_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-6_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-6_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E5 !default: 1.0E5
  MAXIMUM_ITERATIONS=100000 !default: 100000
  RESTART_VALUE=3000 !default: 30
  LINESEARCH_ALPHA=1.0


  !Set all diganostic levels on for testing
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !-----------------------------------------------------------------------------------------------------------
  !COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 1D
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,1,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err) 

  !-----------------------------------------------------------------------------------------------------------
  !REGION
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionLabelSet(Region,"BurgersRegion",Err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)
  
  !-----------------------------------------------------------------------------------------------------------
  !BASIS
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of a basis
  CALL CMISSBasisTypeInitialise(Basis,Err)
  CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
  CALL CMISSBasisTypeSet(Basis,CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(Basis,1,Err)
  !Set the basis xi interpolation and number of Gauss points
  CALL CMISSBasisInterpolationXiSet(Basis,(/CMISSBasisLinearLagrangeInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis,(/2/),Err)
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(Basis,Err)

  !-----------------------------------------------------------------------------------------------------------
  !MESH
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/LENGTH/),Err)
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS/),Err)
  !Finish the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !-----------------------------------------------------------------------------------------------------------  
  !GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL CMISSFieldScalingTypeSet(GeometricField,CMISSFieldNoScaling,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,1,Err)
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)  
  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

  !-----------------------------------------------------------------------------------------------------------  
  !EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the equations_set for a dynamic nonlinear burgers equation
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetFluidMechanicsClass, &
    & CMISSEquationsSetBurgersEquationType,CMISSEquationsSetBurgersSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Set the mesh component to be used by the field components.
  COMPONENT_NUMBER = 1
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
    & COMPONENT_NUMBER,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
    & COMPONENT_NUMBER,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)
  !Initialise dependent field
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & COMPONENT_NUMBER,0.0_CMISSDP,Err)

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIALS FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set material field variables
  CALL CMISSFieldTypeInitialise(MaterialsField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set material field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)
  !Initialise materials field
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,NU_PARAM,Err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set analytic field variables
  !CALL CMISSFieldTypeInitialise(AnalyticField,Err)
  !CALL CMISSEquationsSetAnalyticCreateStart(EquationsSet,CMISSEquationsSetBurgersOneDim1,AnalyticFieldUserNumber, & 
  ! & AnalyticField,Err)
  !Finish the equations set analytic field variables
  !CALL CMISSEquationsSetAnalyticCreateFinish(EquationsSet,Err)

  !-----------------------------------------------------------------------------------------------------------  
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type (Sparse/Full)
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsFullMatrices,Err)
  !Set the equations set output (NoOutput/TimingOutput/MatrixOutput/SolverMatrix/ElementMatrixOutput)
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)

  !Create the equations set boundary conditions
  !CALL CMISSEquationsSetBoundaryConditionsAnalytic(EquationsSet,Err)

  !-----------------------------------------------------------------------------------------------------------
  !PROBLEM
  !-----------------------------------------------------------------------------------------------------------
  !Create the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a static Burgers problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemFluidMechanicsClass,CMISSProblemBurgersEquationType, &
    & CMISSProblemDynamicBurgersSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Create the problem control
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)

  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,DYNAMIC_SOLVER_START_TIME,DYNAMIC_SOLVER_STOP_TIME, & 
    & DYNAMIC_SOLVER_TIME_INCREMENT,Err)
  !Set the output timing
  CALL CMISSControlLoopTimeOutputSet(ControlLoop,DYNAMIC_SOLVER_OUTPUT_FREQUENCY,Err)

  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !-----------------------------------------------------------------------------------------------------------
  !SOLVER
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(DynamicSolver,Err)
  CALL CMISSSolverTypeInitialise(NonlinearSolver,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)

  !Get the dymamic solver
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverUserNumber,DynamicSolver,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(DynamicSolver,DYNAMIC_SOLVER_OUTPUT_TYPE,Err)
  !Set theta
  CALL CMISSSolverDynamicThetaSet(DynamicSolver,DYNAMIC_SOLVER_THETA,Err)

  !Get the dynamic nonlinear solver
  CALL CMISSSolverDynamicNonlinearSolverGet(DynamicSolver,NonlinearSolver,Err)
  !Set the nonlinear Jacobian type
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(NonlinearSolver,CMISSSolverNewtonJacobianAnalyticCalculated,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(NonlinearSolver,NONLINEAR_SOLVER_OUTPUT_TYPE,Err)
  !Set the solver settings
  CALL CMISSSolverNewtonAbsoluteToleranceSet(NonlinearSolver,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolverNewtonRelativeToleranceSet(NonlinearSolver,RELATIVE_TOLERANCE,Err)
  !Get the dynamic nonlinear linear solver
  CALL CMISSSolverNewtonLinearSolverGet(NonlinearSolver,LinearSolver,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(LinearSolver,LINEAR_SOLVER_OUTPUT_TYPE,Err)
  !Set the solver settings

  IF(LINEAR_SOLVER_DIRECT_FLAG) THEN
    CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearDirectSolveType,Err)
    CALL CMISSSolverLibraryTypeSet(LinearSolver,CMISSSolverMUMPSLibrary,Err)
  ELSE
    CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearIterativeSolveType,Err)
    CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolver,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolverLinearIterativeDivergenceToleranceSet(LinearSolver,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeRelativeToleranceSet(LinearSolver,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(LinearSolver,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeGMRESRestartSet(LinearSolver,RESTART_VALUE,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)


  !-----------------------------------------------------------------------------------------------------------
  !SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------
  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the dynamic solver equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,DynamicSolver,Err)
  CALL CMISSSolverSolverEquationsGet(DynamicSolver,SolverEquations,Err)
  !Set the solver equations sparsity (Sparse/Full)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)  
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !-----------------------------------------------------------------------------------------------------------
  !BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------
  !Set up the boundary conditions
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the fixed boundary conditions at the first node and last nodes
  FirstNodeNumber=1
  COMPONENT_NUMBER=1
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSRegionNodesGet(Region,Nodes,Err)
  CALL CMISSNodesNumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL CMISSDecompositionNodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL CMISSDecompositionNodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1, &
      & CMISSNoGlobalDerivative,FirstNodeNumber,COMPONENT_NUMBER,CMISSBoundaryConditionFixed, &
      & 1.0_CMISSDP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1, &
      & CMISSNoGlobalDerivative,LastNodeNumber,COMPONENT_NUMBER,CMISSBoundaryConditionFixed, &
      & 0.0_CMISSDP,Err)
  ENDIF
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)
  !-----------------------------------------------------------------------------------------------------------
  !SOLVE
  !-----------------------------------------------------------------------------------------------------------
  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)

  !-----------------------------------------------------------------------------------------------------------
  !OUTPUT
  !-----------------------------------------------------------------------------------------------------------
  !Output Analytic analysis
  !Call CMISSAnalyticAnalysisOutput(DependentField,"BurgersAnalytics_1D",Err)

  !export fields
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"Burgers_1D","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"Burgers_1D","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
  ENDIF
  
  !Output timing summary
  !CALL TIMING_SUMMARY_OUTPUT(ERR,ERROR,*999)

  !Calculate the stop times and write out the elapsed user and system times
!   CALL CPU_TIMER(USER_CPU,STOP_USER_TIME,ERR,ERROR,*999)
!   CALL CPU_TIMER(SYSTEM_CPU,STOP_SYSTEM_TIME,ERR,ERROR,*999)
! 
!   CALL WRITE_STRING_TWO_VALUE(GENERAL_OUTPUT_TYPE,"User time = ",STOP_USER_TIME(1)-START_USER_TIME(1),", System time = ", &
!     & STOP_SYSTEM_TIME(1)-START_SYSTEM_TIME(1),ERR,ERROR,*999)
!   
  !CALL CMISS_FINALISE(ERR,ERROR,*999)
  !CALL CMISSFinalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  

  STOP
  
END PROGRAM BURGERSEXAMPLE
