!> \file
!> \author David Ladd
!> \brief This is an example program to solve a generalised Burgers's equation using OpenCMISS calls.
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

!> \example ClassicalField/Burgers/src/GeneralisedBurgersExample.f90
!! Example program to solve a generalised Burgers equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/Burgers/history.html
!<

!> Main program
PROGRAM GENERALISEDBURGERSEXAMPLE

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Test program parameters
  
  REAL(CMISSRP), PARAMETER :: LENGTH=3.0_CMISSRP
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_GLOBAL_X_ELEMENTS=6
  REAL(CMISSRP), PARAMETER :: START_TIME=0.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: STOP_TIME=0.1_CMISSRP
  REAL(CMISSRP), PARAMETER :: TIME_INCREMENT=0.01_CMISSRP
    
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: SolverUserNumber=1
  
  !Program variables

  !Program types

  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) ::  AnalyticField,DependentField,EquationsSetField,GeometricField,MaterialsField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh  
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: DynamicSolver,NonlinearSolver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: Err
  LOGICAL :: LINEAR_SOLVER_DIRECT_FLAG
  
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
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)


  CALL cmfe_OutputSetOn("Burgers1DAnalytic",Err)
  
  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !-----------------------------------------------------------------------------------------------------------
  !COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 1D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,1,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err) 

  !-----------------------------------------------------------------------------------------------------------
  !REGION
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"BurgersRegion",Err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)
  
  !-----------------------------------------------------------------------------------------------------------
  !BASIS
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of a basis
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(Basis,1,Err)
  !Set the basis xi interpolation and number of Gauss points
  CALL cmfe_Basis_InterpolationXiSet(Basis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[3],Err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis,Err)

  !-----------------------------------------------------------------------------------------------------------
  !MESH
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH],Err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS],Err)
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !-----------------------------------------------------------------------------------------------------------  
  !GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)
  
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL cmfe_Field_ScalingTypeSet(GeometricField,CMFE_FIELD_NO_SCALING,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)  
  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !-----------------------------------------------------------------------------------------------------------  
  !EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the equations_set for a Generalised Burgers's equation
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & CMFE_EQUATIONS_SET_BURGERS_EQUATION_TYPE,CMFE_EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIALS FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set material field variables
  CALL cmfe_Field_Initialise(MaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set material field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)
  !Initialise materials field
  !Set A
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 1,1.0_CMISSRP,Err)
  !Set B
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 2,-1.0_CMISSRP,Err)
  !Set C
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 3,1.0_CMISSRP,Err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set analytic field variables
  CALL cmfe_Field_Initialise(AnalyticField,Err)
  CALL cmfe_EquationsSet_AnalyticCreateStart(EquationsSet,CMFE_EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_2, &
    & AnalyticFieldUserNumber,AnalyticField,Err)
  !Finish the equations set analytic field variables
  CALL cmfe_EquationsSet_AnalyticCreateFinish(EquationsSet,Err)

  !-----------------------------------------------------------------------------------------------------------  
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type (Sparse/Full)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output (NoOutput/TimingOutput/MatrixOutput/SolverMatrix/ElementMatrixOutput)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !-----------------------------------------------------------------------------------------------------------
  !PROBLEM
  !-----------------------------------------------------------------------------------------------------------
  !Create the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS,CMFE_PROBLEM_BURGERS_EQUATION_TYPE, &
    & CMFE_PROBLEM_DYNAMIC_BURGERS_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,START_TIME,STOP_TIME,TIME_INCREMENT,Err)
  !Set the output timing
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,1,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !-----------------------------------------------------------------------------------------------------------
  !SOLVER
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(DynamicSolver,Err)
  CALL cmfe_Solver_Initialise(NonlinearSolver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)

  !Get the dymamic solver
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverUserNumber,DynamicSolver,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(DynamicSolver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  !Set theta
  CALL cmfe_Solver_DynamicThetaSet(DynamicSolver,0.5_CMISSRP,Err)

  !Get the dynamic nonlinear solver
  CALL cmfe_Solver_DynamicNonlinearSolverGet(DynamicSolver,NonlinearSolver,Err)
  !Set the nonlinear Jacobian type
  !CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(NonlinearSolver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(NonlinearSolver,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  !Set the line search
  CALL cmfe_Solver_NewtonLineSearchTypeSet(NonlinearSolver,CMFE_SOLVER_NEWTON_LINESEARCH_LINEAR,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(NonlinearSolver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  !Get the dynamic nonlinear linear solver
  CALL cmfe_Solver_NewtonLinearSolverGet(NonlinearSolver,LinearSolver,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(LinearSolver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  !Set the solver settings

  LINEAR_SOLVER_DIRECT_FLAG=.FALSE.
  IF(LINEAR_SOLVER_DIRECT_FLAG) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LibraryTypeSet(LinearSolver,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolver,10000,Err)
    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(LinearSolver,50,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)


  !-----------------------------------------------------------------------------------------------------------
  !SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------
  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the dynamic solver equations
  CALL cmfe_Solver_Initialise(DynamicSolver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,DynamicSolver,Err)
  CALL cmfe_Solver_SolverEquationsGet(DynamicSolver,SolverEquations,Err)
  !Set the solver equations sparsity (Sparse/Full)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !-----------------------------------------------------------------------------------------------------------
  !BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------
  !Set up the boundary conditions

  !Create the equations set boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(SolverEquations,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)
  !-----------------------------------------------------------------------------------------------------------
  !SOLVE
  !-----------------------------------------------------------------------------------------------------------
  !Solve the problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !-----------------------------------------------------------------------------------------------------------
  !OUTPUT
  !-----------------------------------------------------------------------------------------------------------
  !Output Analytic analysis
  Call cmfe_AnalyticAnalysis_Output(DependentField,"BurgersAnalytic_1D",Err)

  !export fields
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"Burgers_1D","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"Burgers_1D","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP
  
END PROGRAM GENERALISEDBURGERSEXAMPLE
