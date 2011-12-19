!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a Laplace equation using OpenCMISS calls with objects accessed by user number.
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

!> \example ClassicalField/Laplace/NumberLaplace/src/NumberLaplaceExample.f90
!! Example program to solve a Laplace equation using OpenCMISS calls with objects accessed by user number.
!! \htmlinclude ClassicalField/Laplace/NumberLaplace/history.html
!!
!<

!> Main program
PROGRAM NUMBERLAPLACEEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=3.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=11
  
 
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  
  INTEGER(CMISSIntg) :: MPI_IERROR
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: WorldCoordinateSystemUserNumber
  INTEGER(CMISSIntg) :: WorldRegionUserNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
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

  !Intialise cmiss
  CALL CMISSInitialise(WorldCoordinateSystemUserNumber,WorldRegionUserNumber,Err)
 
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  NUMBER_GLOBAL_X_ELEMENTS=2
  NUMBER_GLOBAL_Y_ELEMENTS=2
  NUMBER_GLOBAL_Z_ELEMENTS=0
  NUMBER_OF_DOMAINS=1
    
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystemUserNumber,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystemUserNumber,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystemUserNumber,Err)

  !Start the creation of the region
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegionUserNumber,Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(RegionUserNumber,CoordinateSystemUserNumber,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(RegionUserNumber,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bilinear Lagrange basis
    CALL CMISSBasis_NumberOfXiSet(BasisUserNumber,2,Err)
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL CMISSBasis_NumberOfXiSet(BasisUserNumber,3,Err)
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(BasisUserNumber,Err)
    
  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,RegionUserNumber,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(RegionUserNumber,GeneratedMeshUserNumber,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(RegionUserNumber,GeneratedMeshUserNumber,BasisUserNumber,Err)
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(RegionUserNumber,GeneratedMeshUserNumber,(/WIDTH,HEIGHT/),Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(RegionUserNumber,GeneratedMeshUserNumber, &
      & (/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/),Err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(RegionUserNumber,GeneratedMeshUserNumber,(/WIDTH,HEIGHT,LENGTH/),Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(RegionUserNumber,GeneratedMeshUserNumber, &
      & (/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS/),Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_CreateFinish(RegionUserNumber,GeneratedMeshUserNumber,MeshUserNumber,Err)
  
  !Create a decomposition
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,RegionUserNumber,MeshUserNumber,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,RegionUserNumber,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(RegionUserNumber,GeometricFieldUserNumber,MeshUserNumber,DecompositionUserNumber,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(RegionUserNumber,GeometricFieldUserNumber,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(RegionUserNumber,GeometricFieldUserNumber,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSField_ComponentMeshComponentSet(RegionUserNumber,GeometricFieldUserNumber,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSField_CreateFinish(RegionUserNumber,GeometricFieldUserNumber,Err)
       
  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(RegionUserNumber,GeneratedMeshUserNumber,GeometricFieldUserNumber,Err)
  
  !Create the equations_set
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,RegionUserNumber,GeometricFieldUserNumber,&
     & CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
     & CMISS_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMISS_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE,EquationsSetFieldUserNumber,&
     & Err)
  !Set the equations set to be a standard Laplace problem
!   CALL CMISSEquationsSet_SpecificationSet(RegionUserNumber,EquationsSetUserNumber,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
!     & CMISS_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMISS_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(RegionUserNumber,EquationsSetUserNumber,Err)

  !Create the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateStart(RegionUserNumber,EquationsSetUserNumber,DependentFieldUserNumber,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(RegionUserNumber,EquationsSetUserNumber,Err)

  !Create the equations set equations
  CALL CMISSEquationsSet_EquationsCreateStart(RegionUserNumber,EquationsSetUserNumber,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(RegionUserNumber,EquationsSetUserNumber,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  !CALL CMISSEquations_OutputTypeSet(RegionUserNumber,EquationsSetUserNumber,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquations_OutputTypeSet(RegionUserNumber,EquationsSetUserNumber,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(RegionUserNumber,EquationsSetUserNumber,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(RegionUserNumber,EquationsSetUserNumber,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(RegionUserNumber,EquationsSetUserNumber,Err)

  !Start the creation of a problem.
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblem_SpecificationSet(ProblemUserNumber,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS,CMISS_PROBLEM_LAPLACE_EQUATION_TYPE, &
    & CMISS_PROBLEM_STANDARD_LAPLACE_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(ProblemUserNumber,Err)

  !Start the creation of the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(ProblemUserNumber,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(ProblemUserNumber,Err)
 
  !Start the creation of the problem solvers
  CALL CMISSProblem_SolversCreateStart(ProblemUserNumber,Err)
  !CALL CMISSSolver_OutputTypeSet(ProblemUserNumber,CMISS_CONTROL_LOOP_NODE,1,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(ProblemUserNumber,CMISS_CONTROL_LOOP_NODE,1,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(ProblemUserNumber,CMISS_CONTROL_LOOP_NODE,1,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(ProblemUserNumber,CMISS_CONTROL_LOOP_NODE,1,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  CALL CMISSSolver_OutputTypeSet(ProblemUserNumber,CMISS_CONTROL_LOOP_NODE,1,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(ProblemUserNumber,Err)

  !Start the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateStart(ProblemUserNumber,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquations_SparsityTypeSet(ProblemUserNumber,CMISS_CONTROL_LOOP_NODE,1,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !CALL CMISSSolverEquations_SparsityTypeSet(ProblemUserNumber,CMISS_CONTROL_LOOP_NODE,1,CMISS_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(ProblemUserNumber,CMISS_CONTROL_LOOP_NODE,1,RegionUserNumber,EquationsSetUserNumber, &
    & EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(ProblemUserNumber,Err)

  !Start the creation of the equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(ProblemUserNumber,CMISS_CONTROL_LOOP_NODE,1,Err)
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
  ELSE
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  ENDIF
  CALL CMISSBoundaryConditions_SetNode(RegionUserNumber,ProblemUserNumber,CMISS_CONTROL_LOOP_NODE,1,DependentFieldUserNumber, &
    & CMISS_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(RegionUserNumber,ProblemUserNumber,CMISS_CONTROL_LOOP_NODE,1,DependentFieldUserNumber, &
    & CMISS_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,1.0_CMISSDP,Err)
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(ProblemUserNumber,CMISS_CONTROL_LOOP_NODE,1,Err)

  !Solve the problem
  CALL CMISSProblem_Solve(ProblemUserNumber,Err)

  !Finialise CMISS
  CALL CMISSFinalise(Err)
 
  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM NUMBERLAPLACEEXAMPLE
