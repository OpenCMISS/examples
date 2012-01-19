!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a linear Poisson equation using openCMISS calls.
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

!> \example ClassicalField/NonlinearPoisson/src/NonlinearPoissonExample.f90
!! Example program to solve a nonlinear Poisson equation using openCMISS calls.
!! \htmlinclude ClassicalField/NonlinearPoisson/history.html
!<

!> Main program
PROGRAM LINEARPOISSONEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField


  !Test program parameters

! Setup used for quads
  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP/2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=1.0_CMISSDP/2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=0.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=13
 
  INTEGER(CMISSIntg), PARAMETER :: CONSTANT_MATERIALS_COMPONENT=3 !! For 3D this is 4
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  REAL(CMISSDP) :: C_PARAM_POISSON

  INTEGER(CMISSIntg) :: MPI_IERROR
  LOGICAL :: EXPORT_FIELD
  

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField,AnalyticField,SourceField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE


#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: Err
  
!  INTEGER(CMISSIntg) :: DIAG_LEVEL_LIST(1)
!  CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(1)

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

  !Set material parameter (constant added to all dofs)
  C_PARAM_POISSON=8.0_CMISSDP  !For 2D


  !Initialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Set diagnostics
!  DIAG_LEVEL_LIST(1)=1
!  DIAG_ROUTINE_LIST(1)="FIELD_MAPPINGS_CALCULATE"

  !CMISS_ALL_DIAG_TYPE/CMISS_IN_DIAG_TYPE/CMISS_FROM_DIAG_TYPE
!  CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,DIAG_LEVEL_LIST,"LinearPoissonDiagnostics",DIAG_ROUTINE_LIST,Err)

  NUMBER_GLOBAL_X_ELEMENTS=2
  NUMBER_GLOBAL_Y_ELEMENTS=2
  NUMBER_GLOBAL_Z_ELEMENTS=0

  RELATIVE_TOLERANCE=1.0E-13_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-13_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E5_CMISSDP !default: 1.0E5
  MAXIMUM_ITERATIONS=10000000_CMISSIntg !default: 100000
  RESTART_VALUE=30_CMISSIntg !default: 30

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMBER_OF_DOMAINS,MPI_IERROR)

  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
 
   !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bilinear Lagrange basis
    CALL CMISSBasis_NumberOfXiSet(Basis,2,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis,(/2,2/),Err) !Quad, linear is default
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,(/3,3/),Err) !3 Gauss per Xi direction per element
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis,(/2,2,2/),Err) !Quad, linear is default
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,(/3,3,3/),Err) !3 Gauss per Xi direction per element
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(Basis,Err)
   
  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/WIDTH,HEIGHT/),Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/),Err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/WIDTH,HEIGHT,LENGTH/),Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS/),Err)
  ENDIF    
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err) !!!What is wrong with this
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)
       
  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
  
  !Create the equations_set
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
    CALL CMISSField_Initialise(EquationsSetField,Err)
CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_POISSON_EQUATION_TYPE,CMISS_EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Set the equations set to be a standard Laplace problem
  
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field variables
  CALL CMISSField_Initialise(MaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,C_PARAM_POISSON,Err)


  !Create the equations set analytic field variables
  CALL CMISSField_Initialise(AnalyticField,Err)
  CALL CMISSEquationsSet_AnalyticCreateStart(EquationsSet,CMISS_TEST_CASE_MIXED_NEUMANN_DIRICHLET_2, &
    & AnalyticFieldUserNumber,AnalyticField,Err)

!Finish the equations set analytic field variables
  CALL CMISSEquationsSet_AnalyticCreateFinish(EquationsSet,Err)
  

  !Create the equations set source field variables for Poisson
  CALL CMISSField_Initialise(SourceField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(EquationsSet,SourceFieldUserNumber, & 
    & SourceField,Err)
!   !Set the mesh component to be used by the field components.
!   CALL CMISSField_ComponentMeshComponentSet(SourceFieldPoisson,CMISS_FIELD_U_VARIABLE_TYPE,MESH_COMPONENT_NUMBER_SPACE, & 
!     & MESH_COMPONENT_NUMBER_PRESSURE,Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSet_SourceCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
!   CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS,CMISS_PROBLEM_POISSON_EQUATION_TYPE, &
    & CMISS_PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)
 
  !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  !Set the solver output
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_SOLVER_OUTPUT,Err)
!   CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  !Set the Jacobian type
  !CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  !CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  !Get the associated linear solver
  !CALL CMISSSolver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  !CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolver,300,Err)
  !Finish the creation of the problem solver

!  CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
!  CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_MUMPS_LIBRARY,Err)

    CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL CMISSSolver_LinearIterativeMaximumIterationsSet(Solver,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolver_LinearIterativeDivergenceToleranceSet(Solver,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(Solver,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(Solver,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeGMRESRestartSet(Solver,RESTART_VALUE,Err)

  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)

  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)

  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Set up the boundary conditions as per the analytic solution
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsAnalytic(SolverEquations,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblem_Solve(Problem,Err)

   !Output Analytic analysis
   Call CMISSAnalyticAnalysisOutput(DependentField,"AnalyticalPoisson",Err)
     
   EXPORT_FIELD=.TRUE. 
   IF(EXPORT_FIELD) THEN
     CALL CMISSFields_Initialise(Fields,Err)
     CALL CMISSFields_Create(Region,Fields,Err)
     CALL CMISSFields_NodesExport(Fields,"NonlinearPoisson","FORTRAN",Err)
     CALL CMISSFields_ElementsExport(Fields,"NonlinearPoisson","FORTRAN",Err)
     CALL CMISSFields_Finalise(Fields,Err)
   ENDIF

  !Finalise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM LINEARPOISSONEXAMPLE
