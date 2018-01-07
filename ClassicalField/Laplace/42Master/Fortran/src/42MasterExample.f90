!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a Laplace equation using OpenCMISS calls.
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

!> \example ClassicalField/Laplace/42Master/src/42MasterExample.f90
!! Example program to solve a Laplace equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/Laplace/42Master/history.html
!<

!> Main program

PROGRAM LAPLACEEXAMPLE

  USE OPENCMISS
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


  !Test program parameters

  REAL(CMFEDP), PARAMETER :: HEIGHT=1.0_CMFEDP
  REAL(CMFEDP), PARAMETER :: WIDTH=2.0_CMFEDP
  REAL(CMFEDP), PARAMETER :: LENGTH=3.0_CMFEDP

  LOGICAL, PARAMETER :: SOLVER_DIRECT_TYPE=.TRUE.

  INTEGER(CMFEIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMFEIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMFEIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMFEIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMFEIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMFEIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMFEIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMFEIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMFEIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMFEIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMFEIntg), PARAMETER :: ProblemUserNumber=11

  !Program types

  !Program variables

  INTEGER(CMFEIntg) :: NUMBER_OF_ARGUMENTS
  INTEGER(CMFEIntg) :: NUMBER_GLOBAL_X_ELEMENTS = -1
  INTEGER(CMFEIntg) :: NUMBER_GLOBAL_Y_ELEMENTS = -1
  INTEGER(CMFEIntg) :: NUMBER_GLOBAL_Z_ELEMENTS = -1
  INTEGER(CMFEIntg) :: INTERPOLATION_TYPE
  INTEGER(CMFEIntg) :: NUMBER_OF_GAUSS_XI
  INTEGER(CMFEIntg) :: I
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT,Filename
  CHARACTER(LEN=255) :: BUFFER

  LOGICAL :: OPTION_1D   = .FALSE.
  LOGICAL :: OPTION_2D   = .FALSE.
  LOGICAL :: OPTION_3D   = .FALSE.
  LOGICAL :: OPTION_TRI  = .FALSE.
  LOGICAL :: OPTION_TET  = .FALSE.
  LOGICAL :: OPTION_QUAD = .FALSE.
  LOGICAL :: OPTION_HEX  = .FALSE.
  LOGICAL :: OPTION_HERMITE        = .FALSE.
  LOGICAL :: OPTION_LINEARBASIS    = .FALSE.
  LOGICAL :: OPTION_QUADRATICBASIS = .FALSE.
  LOGICAL :: OPTION_CUBICBASIS     = .FALSE.


  !CMISS variables

  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables

  INTEGER(CMFEIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMFEIntg) :: EquationsSetIndex
  INTEGER(CMFEIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMFEIntg) :: FirstNodeDomain,LastNodeDomain
  INTEGER(CMFEIntg) :: Err

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

  !Get command line arguments
  I = 1
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  DO WHILE (I <= NUMBER_OF_ARGUMENTS)
    CALL GET_COMMAND_ARGUMENT(I,COMMAND_ARGUMENT)
    SELECT CASE(COMMAND_ARGUMENT)
      CASE('-1D')
        OPTION_1D = .TRUE.
      CASE('-2D')
        OPTION_2D = .TRUE.
      CASE('-3D')
        OPTION_3D = .TRUE.
      CASE('-tri')
        OPTION_TRI = .TRUE.
      CASE('-tet')
        OPTION_TET = .TRUE.
      CASE('-quad')
        OPTION_QUAD = .TRUE.
      CASE('-hex')
        OPTION_HEX = .TRUE.
      CASE('-linearbasis')
        OPTION_LINEARBASIS = .TRUE.
      CASE('-quadraticbasis')
        OPTION_QUADRATICBASIS = .TRUE.
      CASE('-cubicbasis')
        OPTION_CUBICBASIS = .TRUE.
      CASE('-hermite')
        OPTION_HERMITE = .TRUE.
      CASE('-nx')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ (BUFFER,*) NUMBER_GLOBAL_X_ELEMENTS
        I = I + 1
      CASE('-ny')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ (BUFFER,*) NUMBER_GLOBAL_Y_ELEMENTS
        I = I + 1
      CASE('-nz')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ (BUFFER,*) NUMBER_GLOBAL_Z_ELEMENTS
        I = I + 1
      CASE DEFAULT
        WRITE(*,*) 'Unknown argument: ', COMMAND_ARGUMENT
      END SELECT
      I = I + 1
  ENDDO

  !Check for nonsensical arguemnts
  !TODO
  IF (.NOT.(OPTION_3D .OR. OPTION_2D .OR. OPTION_1D)) THEN
    CALL PRINT_USAGE_AND_DIE('Must specify the dimensionality')
  ENDIF

  !option_hermite indicates quad or hex elements
  IF (.NOT.(OPTION_TRI .OR. OPTION_TET .OR. OPTION_QUAD .OR. OPTION_HEX .OR. OPTION_HERMITE .OR. OPTION_1D)) THEN
    CALL PRINT_USAGE_AND_DIE('Must specify the element type')
  ENDIF

  IF (.NOT.(OPTION_LINEARBASIS .OR. OPTION_QUADRATICBASIS .OR. OPTION_CUBICBASIS .OR. OPTION_HERMITE)) THEN
    CALL PRINT_USAGE_AND_DIE('Must specify the basis type')
  ENDIF

  IF (NUMBER_GLOBAL_X_ELEMENTS < 1 .OR. ((OPTION_2D .OR. OPTION_3D) .AND. NUMBER_GLOBAL_Y_ELEMENTS < 1) .OR.  &
      & (OPTION_3D .AND. NUMBER_GLOBAL_Z_ELEMENTS < 0)) THEN
    CALL PRINT_USAGE_AND_DIE('Must specify number of elements')
  ENDIF

  !Interpret arguments

  IF(OPTION_1D) THEN
    NUMBER_GLOBAL_Y_ELEMENTS = 0
    NUMBER_GLOBAL_Z_ELEMENTS = 0
  ELSEIF(OPTION_2D) THEN
    NUMBER_GLOBAL_Z_ELEMENTS = 0
  ENDIF

  NUMBER_OF_GAUSS_XI=0
  IF(OPTION_LINEARBASIS) THEN
    IF(OPTION_TRI .OR. OPTION_TET) THEN
      INTERPOLATION_TYPE = 7
    ELSE
      INTERPOLATION_TYPE = 1
      NUMBER_OF_GAUSS_XI=2
    ENDIF
  ELSEIF(OPTION_QUADRATICBASIS) THEN
    IF(OPTION_TRI .OR. OPTION_TET) THEN
      INTERPOLATION_TYPE = 8
    ELSE
      INTERPOLATION_TYPE = 2
      NUMBER_OF_GAUSS_XI=3
    ENDIF
  ELSEIF(OPTION_HERMITE) THEN
    INTERPOLATION_TYPE = 4
    NUMBER_OF_GAUSS_XI=4
  ELSEIF(OPTION_CUBICBASIS) THEN
    IF(OPTION_TRI .OR. OPTION_TET) THEN
      INTERPOLATION_TYPE = 9
    ELSE
      INTERPOLATION_TYPE = 3
      NUMBER_OF_GAUSS_XI=4
    ENDIF
  ELSE
    CALL HANDLE_ERROR("Could not set interploation type")
  ENDIF


  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  WRITE(Filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "Laplace",NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS,INTERPOLATION_TYPE

  CALL cmfe_OutputSetOn(Filename,Err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  IF(OPTION_1D) THEN
    !Set the coordinate system to be 1D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,1,Err)
  ELSEIF(OPTION_2D) THEN
    !Set the coordinate system to be 2D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"LaplaceRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  IF(OPTION_TRI .OR. OPTION_TET) THEN !Default is Lagrange Hermite type
    CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  ENDIF
  IF(OPTION_1D) THEN
    !Set the basis to be a linear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis,1,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE],Err)
    IF(.NOT. OPTION_TRI .AND. .NOT. OPTION_TET) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ELSEIF(OPTION_2D) THEN
    !Set the basis to be a bilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis,2,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(.NOT. OPTION_TRI .AND. .NOT. OPTION_TET) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis,3,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(.NOT. OPTION_TRI .AND. .NOT. OPTION_TET) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis,Err)

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)
  !Define the mesh on the region
  IF(OPTION_1D) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS],Err)
  ELSEIF(OPTION_2D) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

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
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  IF(NUMBER_GLOBAL_Y_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  ENDIF
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Create the Standard Laplace Equations set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMFE_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Set the DOFs to be contiguous across components
  CALL cmfe_Field_DOFOrderTypeSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,Err)
  CALL cmfe_Field_DOFOrderTypeSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Initialise the field with an initial guess
  IF(.NOT. SOLVER_DIRECT_TYPE) THEN
    CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      & 0.5_CMFEDP, &
      & Err)
  ENDIF

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS,CMFE_PROBLEM_LAPLACE_EQUATION_TYPE, &
    & CMFE_PROBLEM_STANDARD_LAPLACE_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  IF(SOLVER_DIRECT_TYPE) THEN
    CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_MUMPS_LIBRARY,Err)
    !CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_SUPERLU_LIBRARY,Err)
    !CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_PASTIX_LIBRARY,Err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-10_CMFEDP,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(Solver,1.0E-10_CMFEDP,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(Solver,100000,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Start the creation of the equations set boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Region_NodesGet(Region,Nodes,Err)
  CALL cmfe_Nodes_NumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMFEDP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMFEDP,Err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !Export results
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"Laplace","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"Laplace","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR

  SUBROUTINE PRINT_USAGE_AND_DIE(ERROR_STRING)
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE (*,*) ERROR_STRING
    WRITE (*,*) ''
    WRITE (*,*) 'Usage:'
    WRITE (*,*) 'The options fall into the following groups:'
    WRITE (*,*) '([*] indicates the default option)'
    WRITE (*,*) 'Dimension:    -1D/-2D/-3D                              [2D]'
    WRITE (*,*) 'Element Type: -tri/-tet/-quad/-hex/-hermite            [quad/hex]'
    WRITE (*,*) 'Basis Type:   -linearbasis/-quadraticbasis/-cubicbasis [linearbasis]'

    WRITE (*,*) ''
    WRITE (*,*) 'Furthermore, the user must specify the number of elements in each '
    WRITE (*,*) 'dimension, for the appropriate number of dimensions. The format is:'
    !WRITE (*,*) ''
    WRITE (*,*) '-nx # -ny # -nz #'

    CALL EXIT(1)

  END SUBROUTINE PRINT_USAGE_AND_DIE

END PROGRAM LAPLACEEXAMPLE
