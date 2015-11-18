!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve an Analytic Laplace equation using OpenCMISS calls.
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

!> \example ClassicalField/Laplace/AnalyticLaplace/src/AnalyticLaplaceExample.f90
!! Example illustrating the use of OpenCMISS to solve the Laplace problem and check with its Analytic Solution.
!! 
!! \htmlinclude ClassicalField/Laplace/AnalyticLaplace/history.html
!< 

!> Main program
PROGRAM ANALYTICLAPLACEEXAMPLE
#ifndef NOMPIMOD
  USE MPI
#endif

  USE OpenCMISS
  USE OpenCMISS_Iron

  USE TEST_FRAMEWORK_ROUTINES

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  REAL(CMISSRP), PARAMETER :: ORIGIN(2)=[-3.141592653579_CMISSRP/2.0_CMISSRP, -3.141592653579_CMISSRP/2.0_CMISSRP]
  REAL(CMISSRP), PARAMETER :: HEIGHT=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=2.0_CMISSRP

  !Program types

  !Program variables

  TYPE(cmfe_RegionType) :: WORLD_REGION
  TYPE(cmfe_CoordinateSystemType) :: WorldCoordinateSystem

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS,INTERPOLATION
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables
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
  CALL cmfe_Initialise(WorldCoordinateSystem,WORLD_REGION,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  CALL cmfe_RandomSeedsSet(9999,Err)
  
  CALL cmfe_DiagnosticsSetOn(CMFE_ALL_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",[""],Err)

  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) INTERPOLATION

    SELECT CASE(INTERPOLATION)
    CASE(1)
      CALL ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(2,6,2)
      CALL ANALYTICLAPLACE_TESTCASE_LINEAR_LAGRANGE_EXPORT(2,2,0)
    CASE(4)
      CALL ANALYTICLAPLACE_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(2,10,2)
      CALL ANALYTICLAPLACE_TESTCASE_CUBIC_HERMITE_EXPORT(2,2,0)
    CASE(7)
      CALL ANALYTICLAPLACE_TESTCASE_BILINEAR_SIMPLEX_CONVERGENCE(2,6,2)
      CALL ANALYTICLAPLACE_TESTCASE_LINEAR_SIMPLEX_EXPORT(2,2,0)
    CASE DEFAULT
      CALL HANDLE_ERROR("Invalid interpolation specified.")
    END SELECT
  ELSE
    !Run all tests
    CALL ANALYTICLAPLACE_TESTCASE_BILINEAR_SIMPLEX_CONVERGENCE(2,6,2)
    CALL ANALYTICLAPLACE_TESTCASE_LINEAR_SIMPLEX_EXPORT(2,2,0)
    CALL ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(2,6,2)
    CALL ANALYTICLAPLACE_TESTCASE_LINEAR_LAGRANGE_EXPORT(2,2,0)
    CALL ANALYTICLAPLACE_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(2,10,2)
    CALL ANALYTICLAPLACE_TESTCASE_CUBIC_HERMITE_EXPORT(2,2,0)
  ENDIF

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

CONTAINS

  !
  !================================================================================================================================
  !  

  !>Export analytic analyis
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_CUBIC_HERMITE_EXPORT(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<increment interval number of elements per axis
    !Local Variables
    TYPE(cmfe_FieldType) :: FIELD

    CALL cmfe_Field_Initialise(FIELD,Err)
    CALL ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,4, &
      & FIELD)

    CALL cmfe_AnalyticAnalysis_Output(FIELD,"AnalyticLaplaceCubicHermite",Err)
    
    CALL ANALYTICLAPLACE_GENERIC_CLEAN(1,1,1,1,1)

  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_CUBIC_HERMITE_EXPORT
  
  !
  !================================================================================================================================
  !
  
  !>Export analytic analyis
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_LINEAR_LAGRANGE_EXPORT(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<increment interval number of elements per axis
    !Local Variables
    TYPE(cmfe_FieldType) :: FIELD

    CALL cmfe_Field_Initialise(FIELD,Err)
    CALL ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,1, &
      & FIELD)

    CALL cmfe_AnalyticAnalysis_Output(FIELD,"AnalyticLaplaceLinearLagrange",Err)
    
    CALL ANALYTICLAPLACE_GENERIC_CLEAN(1,1,1,1,1)

  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_LINEAR_LAGRANGE_EXPORT
  
  !
  !================================================================================================================================
  !
  
  !>Export analytic analyis
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_LINEAR_SIMPLEX_EXPORT(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<increment interval number of elements per axis
    !Local Variables
    TYPE(cmfe_FieldType) :: FIELD

    CALL cmfe_Field_Initialise(FIELD,Err)
    CALL ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,7, &
      & FIELD)

    CALL cmfe_AnalyticAnalysis_Output(FIELD,"AnalyticLaplaceLinearSimplex",Err)
    
    CALL ANALYTICLAPLACE_GENERIC_CLEAN(1,1,1,1,1)

  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_LINEAR_SIMPLEX_EXPORT
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear Simplex interpolation is expected.
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_BILINEAR_SIMPLEX_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(CMISSRP) :: VALUE
    REAL(CMISSRP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    
    CALL ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,7,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)

    CALL TEST_FRAMEWORK_ASSERT_EQUALS(2.0_CMISSRP,VALUE,0.5_CMISSRP,ERR)
    
    WRITE(*,'(A)') "Analytic Laplace Example Testcase 1 - bilinear Simplex has successfully completed."
    
  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_BILINEAR_SIMPLEX_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(CMISSRP) :: VALUE
    REAL(CMISSRP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    
    CALL ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,1,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)

    CALL TEST_FRAMEWORK_ASSERT_EQUALS(2.0_CMISSRP,VALUE,0.5_CMISSRP,ERR)
    
    WRITE(*,'(A)') "Analytic Laplace Example Testcase 2 - bilinear Lagrange has successfully completed."
    
  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(CMISSRP) :: VALUE
    REAL(CMISSRP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)

    !Note INTERPOLATION_SPECIFICATIONS of 4 is Cubic Hermite
    CALL ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,4,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)
    !This test is superconvergent so look for a slope of 5 rather than 4. Should really test >= 4
    CALL TEST_FRAMEWORK_ASSERT_EQUALS(5.0_CMISSRP,VALUE,1.0_CMISSRP,Err)
    IF (Err/=0) THEN
      WRITE(*,'(A,F6.3)') "Analytic Laplace Example Testcase 3 - bicubic Hermite failure: Convergence should be around 4.0" &
        & //", but it was ", VALUE
    ENDIF
    WRITE(*,'(A)') "Analytic Laplace Example Testcase 3 - bicubic Hermite has successfully completed."
    
  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_BICUBIC_HERMITE_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of the specified interpolation is expected.
  SUBROUTINE ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
    & NUMBER_OF_ELEMENTS_XI_INTERVAL,INTERPOLATION_SPECIFICATIONS,X_VALUES,Y_VALUES)
  
    !Argument variables 
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: INTERPOLATION_SPECIFICATIONS !<interpolation specifications
    REAL(CMISSRP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    !Local Variables
    REAL(CMISSRP) :: VALUE
    
    INTEGER(CMISSIntg) :: i
    TYPE(cmfe_FieldType) :: FIELD
    
    ALLOCATE(X_VALUES((NUMBER_OF_ELEMENTS_XI_END-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1),STAT=ERR)
    ALLOCATE(Y_VALUES((NUMBER_OF_ELEMENTS_XI_END-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1),STAT=ERR)

    DO i = NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL

      CALL cmfe_Field_Initialise(FIELD,Err)
      CALL ANALYTICLAPLACE_GENERIC(i,i,0,INTERPOLATION_SPECIFICATIONS,FIELD)
      CALL cmfe_AnalyticAnalysis_AbsoluteErrorGetNode(FIELD,1,1,1,(i+1)**2/2+1,1,VALUE,Err)

      Y_VALUES((i-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1)=log10(VALUE)
      X_VALUES((i-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1)=log10(HEIGHT/i)
      CALL ANALYTICLAPLACE_GENERIC_CLEAN(1,1,1,1,1)
   
    ENDDO
  END SUBROUTINE ANALYTICLAPLACE_GENERIC_CONVERGENCE
  
  
  !
  !================================================================================================================================
  !   
    
  SUBROUTINE ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
    & INTERPOLATION_SPECIFICATIONS,DEPENDENT_FIELD)
    !Argument variables 
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<number of elements on x axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<number of elements on y axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<number of elements on z axis
    INTEGER(CMISSIntg), INTENT(IN) :: INTERPOLATION_SPECIFICATIONS !<the interpolation specifications
    TYPE(cmfe_FieldType) :: DEPENDENT_FIELD
    !Local Variables
    INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
    INTEGER(CMISSIntg) :: MPI_IERROR
    INTEGER(CMISSIntg) :: ANALYTIC_FUNCTION
    INTEGER(CMISSIntg) :: EquationsSetIndex

    TYPE(cmfe_BasisType) :: Basis
    TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
    TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
    TYPE(cmfe_GeneratedMeshType) :: GENERATED_MESH
    TYPE(cmfe_MeshType) :: MESH
    TYPE(cmfe_DecompositionType) :: DECOMPOSITION
    TYPE(cmfe_EquationsType) :: EQUATIONS
    TYPE(cmfe_EquationsSetType) :: EQUATIONS_SET
    TYPE(cmfe_FieldType) :: ANALYTIC_FIELD,GEOMETRIC_FIELD,EquationsSetField
    TYPE(cmfe_ProblemType) :: PROBLEM
    TYPE(cmfe_RegionType) :: REGION
    TYPE(cmfe_SolverType) :: SOLVER
    TYPE(cmfe_SolverEquationsType) :: SOLVER_EQUATIONS
    
    NUMBER_OF_DOMAINS=1

    !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
    CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(INTERPOLATION_SPECIFICATIONS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

    !Start the creation of a new RC coordinate system
    CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
    CALL cmfe_CoordinateSystem_CreateStart(1,CoordinateSystem,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the coordinate system to be 2D
      CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
    ELSE
      !Set the coordinate system to be 3D
      CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
    ENDIF
    !Finish the creation of the coordinate system
    CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

    !Start the creation of the region
    CALL cmfe_Region_Initialise(REGION,Err)
    CALL cmfe_Region_CreateStart(1,WORLD_REGION,REGION,Err)
    !Set the regions coordinate system to the 2D RC coordinate system that we have created
    CALL cmfe_Region_CoordinateSystemSet(REGION,CoordinateSystem,Err)
    !Finish the creation of the region
    CALL cmfe_Region_CreateFinish(REGION,Err)

  
    !Start the creation of a basis (default is trilinear lagrange)
    CALL cmfe_Basis_Initialise(Basis,Err)
    CALL cmfe_Basis_CreateStart(1,Basis,Err)
    SELECT CASE(INTERPOLATION_SPECIFICATIONS)
    CASE(CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
      & CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
      & CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION)
      !Do nothing
    CASE(CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION,CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
      & CMFE_BASIS_CUBIC_SIMPLEX_INTERPOLATION)
      CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_SIMPLEX_TYPE,Err)
    CASE DEFAULT
      WRITE(*,'(A)') "Invalid interpolation specification."
      STOP
    END SELECT
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the basis to be a bilinear basis
      CALL cmfe_Basis_NumberOfXiSet(Basis,2,Err)
      CALL cmfe_Basis_InterpolationXiSet(Basis,[INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS],Err)
    ELSE
      !Set the basis to be a trilinear basis
      CALL cmfe_Basis_NumberOfXiSet(Basis,3,Err)
      CALL cmfe_Basis_InterpolationXiSet(Basis,[INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS, &
          & INTERPOLATION_SPECIFICATIONS],Err)
    ENDIF
    !Set the number of Gauss points
    IF(INTERPOLATION_SPECIFICATIONS==CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[3,3],Err)
    ENDIF
    !Finish the creation of the basis
    CALL cmfe_Basis_CreateFinish(Basis,Err)

    !Start the creation of a generated mesh in the region
    CALL cmfe_GeneratedMesh_Initialise(GENERATED_MESH,Err)
    CALL cmfe_GeneratedMesh_CreateStart(1,REGION,GENERATED_MESH,Err)
    !Set up a regular 100x100 mesh
    CALL cmfe_GeneratedMesh_TypeSet(GENERATED_MESH,1,Err)
    CALL cmfe_GeneratedMesh_BasisSet(GENERATED_MESH,Basis,Err)
    !Define the mesh on the region
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      CALL cmfe_GeneratedMesh_ExtentSet(GENERATED_MESH,[WIDTH,HEIGHT],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GENERATED_MESH,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS], &
        & Err)
      CALL cmfe_GeneratedMesh_OriginSet(GENERATED_MESH,ORIGIN,Err)
    ELSE
      CALL cmfe_GeneratedMesh_ExtentSet(GENERATED_MESH,[WIDTH,HEIGHT,LENGTH],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GENERATED_MESH,[NUMBER_GLOBAL_X_ELEMENTS, &
        & NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS],Err)
    ENDIF
    !Finish the creation of a generated mesh in the region
    CALL cmfe_Mesh_Initialise(MESH,Err)
    CALL cmfe_GeneratedMesh_CreateFinish(GENERATED_MESH,1,MESH,Err)
    
    !Create a decomposition
    CALL cmfe_Decomposition_Initialise(DECOMPOSITION,Err)
    CALL cmfe_Decomposition_CreateStart(1,MESH,DECOMPOSITION,Err)
    !Set the decomposition to be a general decomposition with the specified number of domains
    CALL cmfe_Decomposition_TypeSet(DECOMPOSITION,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL cmfe_Decomposition_NumberOfDomainsSet(DECOMPOSITION,NUMBER_OF_DOMAINS,Err)
    CALL cmfe_Decomposition_CreateFinish(DECOMPOSITION,Err)

    !Start to create a default (geometric) field on the region
    CALL cmfe_Field_Initialise(GEOMETRIC_FIELD,Err)
    CALL cmfe_Field_CreateStart(1,REGION,GEOMETRIC_FIELD,Err)
    !Set the decomposition to use
    CALL cmfe_Field_MeshDecompositionSet(GEOMETRIC_FIELD,DECOMPOSITION,Err)
    !Set the domain to be used by the field components
    !NB these are needed now as the default mesh component number is 1
    CALL cmfe_Field_ComponentMeshComponentSet(GEOMETRIC_FIELD,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(GEOMETRIC_FIELD,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      CALL cmfe_Field_ComponentMeshComponentSet(GEOMETRIC_FIELD,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
    ENDIF
    IF(INTERPOLATION_SPECIFICATIONS==CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL cmfe_Field_ScalingTypeSet(GEOMETRIC_FIELD,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
    ENDIF
    !Finish creating the field
    CALL cmfe_Field_CreateFinish(GEOMETRIC_FIELD,Err)

    !Update the geometric field parameters
    CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GENERATED_MESH,GEOMETRIC_FIELD,Err)

    !Create the equations_set
    CALL cmfe_EquationsSet_Initialise(EQUATIONS_SET,Err)
    CALL cmfe_Field_Initialise(EquationsSetField,Err)
    !Set the equations set to be a standard Laplace problem
    CALL cmfe_EquationsSet_CreateStart(1,REGION,GEOMETRIC_FIELD,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
      & CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMFE_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE],8, &
      & EquationsSetField,EQUATIONS_SET,Err)
    
    !Finish creating the equations set
    CALL cmfe_EquationsSet_CreateFinish(EQUATIONS_SET,Err)
  
    !Create the equations set dependent field variables
    CALL cmfe_Field_Initialise(DEPENDENT_FIELD,Err)
    CALL cmfe_EquationsSet_DependentCreateStart(EQUATIONS_SET,2,DEPENDENT_FIELD,Err)
    !Finish the equations set dependent field variables
    CALL cmfe_EquationsSet_DependentCreateFinish(EQUATIONS_SET,Err)

    !Create the equations set analytic field variables
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      ANALYTIC_FUNCTION=CMFE_EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2
    ELSE
      ANALYTIC_FUNCTION=CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2
    ENDIF
    CALL cmfe_Field_Initialise(ANALYTIC_FIELD,Err)
    CALL cmfe_EquationsSet_AnalyticCreateStart(EQUATIONS_SET,ANALYTIC_FUNCTION,3,ANALYTIC_FIELD,Err)
    !Finish the equations set analtyic field variables
    CALL cmfe_EquationsSet_AnalyticCreateFinish(EQUATIONS_SET,Err)

    !Create the equations set equations
    CALL cmfe_Equations_Initialise(EQUATIONS,Err)
    CALL cmfe_EquationsSet_EquationsCreateStart(EQUATIONS_SET,EQUATIONS,Err)
    !Set the equations matrices sparsity type
    CALL cmfe_Equations_SparsityTypeSet(EQUATIONS,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
    CALL cmfe_EquationsSet_EquationsCreateFinish(EQUATIONS_SET,Err)

    !Create the problem
    CALL cmfe_Problem_Initialise(PROBLEM,Err)
    CALL cmfe_Problem_CreateStart(1,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS,CMFE_PROBLEM_LAPLACE_EQUATION_TYPE, &
      & CMFE_PROBLEM_STANDARD_LAPLACE_SUBTYPE],PROBLEM,Err)
    !Finish creating the problem
    CALL cmfe_Problem_CreateFinish(PROBLEM,Err)

    !Create the problem control loop
    CALL cmfe_Problem_ControlLoopCreateStart(PROBLEM,Err)
    !Finish creating the problem control
    CALL cmfe_Problem_ControlLoopCreateFinish(PROBLEM,Err)

    !Start the creation of the problem solvers
    CALL cmfe_Solver_Initialise(Solver,Err)
    CALL cmfe_Problem_SolversCreateStart(Problem,Err)
    CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
    !Set solver to direct type
    CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-12_CMISSRP,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(Solver,1.0E-12_CMISSRP,Err)
    !CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    !CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_MUMPS_LIBRARY,Err)
    !Finish the creation of the problem solver
    CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

    !Start the creation of the problem solver equations
    CALL cmfe_Solver_Initialise(Solver,Err)
    CALL cmfe_SolverEquations_Initialise(Solver_Equations,Err)
    CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
    !Get the solve equations
    CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
    CALL cmfe_Solver_SolverEquationsGet(Solver,Solver_Equations,Err)
    !Set the solver equations sparsity
    CALL cmfe_SolverEquations_SparsityTypeSet(Solver_Equations,CMFE_SOLVER_SPARSE_MATRICES,Err)
    !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)
    !Add in the equations set
    CALL cmfe_SolverEquations_EquationsSetAdd(Solver_Equations,Equations_Set,EquationsSetIndex,Err)
    !Finish the creation of the problem solver equations
    CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

    !Set up the boundary conditions as per the analytic solution
    CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
    CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(Solver_Equations,BoundaryConditions,Err)
    CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(Solver_Equations,Err)
    CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(Solver_Equations,Err)

    !Solve the problem
    CALL cmfe_Problem_Solve(PROBLEM,Err)

  END SUBROUTINE ANALYTICLAPLACE_GENERIC

  SUBROUTINE ANALYTICLAPLACE_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber,GeneratedMeshUserNumber, &
    & ProblemUserNumber)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: CoordinateSystemUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: RegionUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: BasisUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: GeneratedMeshUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: ProblemUserNumber

    CALL cmfe_Problem_Destroy(ProblemUserNumber,Err)
    CALL cmfe_GeneratedMesh_Destroy(RegionUserNumber,GeneratedMeshUserNumber,Err)
    CALL cmfe_Basis_Destroy(BasisUserNumber,Err)
    CALL cmfe_Region_Destroy(RegionUserNumber,Err)
    CALL cmfe_CoordinateSystem_Destroy(CoordinateSystemUserNumber,Err)

  END SUBROUTINE ANALYTICLAPLACE_GENERIC_CLEAN

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR

END PROGRAM ANALYTICLAPLACEEXAMPLE 
