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

  USE MPI
  USE OPENCMISS
  USE TEST_FRAMEWORK_ROUTINES

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField


  !Test program parameters

  REAL(CMISSDP), PARAMETER :: ORIGIN(2)=[-3.141592653579_CMISSDP/2, -3.141592653579_CMISSDP/2]
  REAL(CMISSDP), PARAMETER :: HEIGHT=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=2.0_CMISSDP

  !Program types

  !Program variables

  TYPE(CMISSRegionType) :: WORLD_REGION
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem

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
  CALL CMISSInitialise(WorldCoordinateSystem,WORLD_REGION,Err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  CALL CMISSRandomSeedsSet(9999,Err)
  
  CALL CMISSDiagnosticsSetOn(CMISS_ALL_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",[""],Err)

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

  CALL CMISSFinalise(Err)

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
    TYPE(CMISSFieldType) :: FIELD

    CALL CMISSField_Initialise(FIELD,Err)
    CALL ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,4, &
      & FIELD)

    CALL CMISSAnalyticAnalysisOutput(FIELD,"AnalyticLaplaceCubicHermite",Err)
    
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
    TYPE(CMISSFieldType) :: FIELD

    CALL CMISSField_Initialise(FIELD,Err)
    CALL ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,1, &
      & FIELD)

    CALL CMISSAnalyticAnalysisOutput(FIELD,"AnalyticLaplaceLinearLagrange",Err)
    
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
    TYPE(CMISSFieldType) :: FIELD

    CALL CMISSField_Initialise(FIELD,Err)
    CALL ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,7, &
      & FIELD)

    CALL CMISSAnalyticAnalysisOutput(FIELD,"AnalyticLaplaceLinearSimplex",Err)
    
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
    REAL(CMISSDP) :: VALUE
    REAL(CMISSDP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    
    CALL ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,7,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)

    CALL TEST_FRAMEWORK_ASSERT_EQUALS(2.0_CMISSDP,VALUE,0.5_CMISSDP,ERR)
    
    WRITE(*,'(A)') "Analytic Laplace Example Testcase1 - bilinear Simplex is successfully completed."
    
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
    REAL(CMISSDP) :: VALUE
    REAL(CMISSDP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    
    CALL ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,1,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)

    CALL TEST_FRAMEWORK_ASSERT_EQUALS(2.0_CMISSDP,VALUE,0.5_CMISSDP,ERR)
    
    WRITE(*,'(A)') "Analytic Laplace Example Testcase2 - bilinear lagrange is successfully completed."
    
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
    REAL(CMISSDP) :: VALUE
    REAL(CMISSDP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)

    !Note INTERPOLATION_SPECIFICATIONS of 4 is Cubic Hermite
    CALL ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,4,X_VALUES,Y_VALUES)
    
   CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)
   CALL TEST_FRAMEWORK_ASSERT_EQUALS(4.0_CMISSDP,VALUE,1.0_CMISSDP,Err)
   IF (Err/=0) THEN
     WRITE(*,'(A,F6.3)') "Analytic Laplace Example Testcase3 - bicubic Hermite failure: Convergence should be around 4.0" &
       & //", but it was ", VALUE
   ENDIF
   WRITE(*,'(A)') "Analytic Laplace Example Testcase3 - bicubic Hermite is successfully completed."

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
    REAL(CMISSDP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    !Local Variables
    REAL(CMISSDP) :: VALUE
    
    INTEGER(CMISSIntg) :: i
    TYPE(CMISSFieldType) :: FIELD
    
    ALLOCATE(X_VALUES((NUMBER_OF_ELEMENTS_XI_END-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1),STAT=ERR)
    ALLOCATE(Y_VALUES((NUMBER_OF_ELEMENTS_XI_END-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1),STAT=ERR)

    DO i = NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL

      CALL CMISSField_Initialise(FIELD,Err)
      CALL ANALYTICLAPLACE_GENERIC(i,i,0,INTERPOLATION_SPECIFICATIONS,FIELD)
      CALL CMISSAnalyticAnalysisAbsoluteErrorGetNode(FIELD,1,1,1,(i+1)**2/2+1,1,VALUE,Err)

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
    TYPE(CMISSFieldType) :: DEPENDENT_FIELD
    !Local Variables
    INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
    INTEGER(CMISSIntg) :: MPI_IERROR
    INTEGER(CMISSIntg) :: ANALYTIC_FUNCTION
    INTEGER(CMISSIntg) :: EquationsSetIndex

    TYPE(CMISSBasisType) :: Basis
    TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
    TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
    TYPE(CMISSGeneratedMeshType) :: GENERATED_MESH
    TYPE(CMISSMeshType) :: MESH
    TYPE(CMISSDecompositionType) :: DECOMPOSITION
    TYPE(CMISSEquationsType) :: EQUATIONS
    TYPE(CMISSEquationsSetType) :: EQUATIONS_SET
    TYPE(CMISSFieldType) :: ANALYTIC_FIELD,GEOMETRIC_FIELD
    TYPE(CMISSProblemType) :: PROBLEM
    TYPE(CMISSRegionType) :: REGION
    TYPE(CMISSSolverType) :: SOLVER
    TYPE(CMISSSolverEquationsType) :: SOLVER_EQUATIONS
    
    NUMBER_OF_DOMAINS=1

    !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
    CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(INTERPOLATION_SPECIFICATIONS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

    !Start the creation of a new RC coordinate system
    CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
    CALL CMISSCoordinateSystem_CreateStart(1,CoordinateSystem,Err)
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
    CALL CMISSRegion_Initialise(REGION,Err)
    CALL CMISSRegion_CreateStart(1,WORLD_REGION,REGION,Err)
    !Set the regions coordinate system to the 2D RC coordinate system that we have created
    CALL CMISSRegion_CoordinateSystemSet(REGION,CoordinateSystem,Err)
    !Finish the creation of the region
    CALL CMISSRegion_CreateFinish(REGION,Err)

  
    !Start the creation of a basis (default is trilinear lagrange)
    CALL CMISSBasis_Initialise(Basis,Err)
    CALL CMISSBasis_CreateStart(1,Basis,Err)
    SELECT CASE(INTERPOLATION_SPECIFICATIONS)
    CASE(CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
      & CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
      & CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)
      !Do nothing
    CASE(CMISS_BASIS_LINEAR_SIMPLEX_INTERPOLATION,CMISS_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
      & CMISS_BASIS_CUBIC_SIMPLEX_INTERPOLATION)
      CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_SIMPLEX_TYPE,Err)
    CASE DEFAULT
      WRITE(*,'(A)') "Invalid interpolation specification."
      STOP
    END SELECT
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the basis to be a bilinear basis
      CALL CMISSBasis_NumberOfXiSet(Basis,2,Err)
      CALL CMISSBasis_InterpolationXiSet(Basis,[INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS],Err)
    ELSE
      !Set the basis to be a trilinear basis
      CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
      CALL CMISSBasis_InterpolationXiSet(Basis,[INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS, &
          & INTERPOLATION_SPECIFICATIONS],Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasis_CreateFinish(Basis,Err)

    !Start the creation of a generated mesh in the region
    CALL CMISSGeneratedMesh_Initialise(GENERATED_MESH,Err)
    CALL CMISSGeneratedMesh_CreateStart(1,REGION,GENERATED_MESH,Err)
    !Set up a regular 100x100 mesh
    CALL CMISSGeneratedMesh_TypeSet(GENERATED_MESH,1,Err)
    CALL CMISSGeneratedMesh_BasisSet(GENERATED_MESH,Basis,Err)
    !Define the mesh on the region
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      CALL CMISSGeneratedMesh_ExtentSet(GENERATED_MESH,[WIDTH,HEIGHT],Err)
      CALL CMISSGeneratedMesh_NumberOfElementsSet(GENERATED_MESH,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS], &
        & Err)
      CALL CMISSGeneratedMesh_OriginSet(GENERATED_MESH,ORIGIN,Err)
    ELSE
      CALL CMISSGeneratedMesh_ExtentSet(GENERATED_MESH,[WIDTH,HEIGHT,LENGTH],Err)
      CALL CMISSGeneratedMesh_NumberOfElementsSet(GENERATED_MESH,[NUMBER_GLOBAL_X_ELEMENTS, &
        & NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS],Err)
    ENDIF
    !Finish the creation of a generated mesh in the region
    CALL CMISSMesh_Initialise(MESH,Err)
    CALL CMISSGeneratedMesh_CreateFinish(GENERATED_MESH,1,MESH,Err)
    
    !Create a decomposition
    CALL CMISSDecomposition_Initialise(DECOMPOSITION,Err)
    CALL CMISSDecomposition_CreateStart(1,MESH,DECOMPOSITION,Err)
    !Set the decomposition to be a general decomposition with the specified number of domains
    CALL CMISSDecomposition_TypeSet(DECOMPOSITION,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL CMISSDecomposition_NumberOfDomainsSet(DECOMPOSITION,NUMBER_OF_DOMAINS,Err)
    CALL CMISSDecomposition_CreateFinish(DECOMPOSITION,Err)

    !Start to create a default (geometric) field on the region
    CALL CMISSField_Initialise(GEOMETRIC_FIELD,Err)
    CALL CMISSField_CreateStart(1,REGION,GEOMETRIC_FIELD,Err)
    !Set the decomposition to use
    CALL CMISSField_MeshDecompositionSet(GEOMETRIC_FIELD,DECOMPOSITION,Err)
    !Set the domain to be used by the field components
    !NB these are needed now as the default mesh component number is 1
    CALL CMISSField_ComponentMeshComponentSet(GEOMETRIC_FIELD,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
    CALL CMISSField_ComponentMeshComponentSet(GEOMETRIC_FIELD,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      CALL CMISSField_ComponentMeshComponentSet(GEOMETRIC_FIELD,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
    ENDIF
    IF(INTERPOLATION_SPECIFICATIONS==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL CMISSField_ScalingTypeSet(GEOMETRIC_FIELD,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
    ENDIF
    !Finish creating the field
    CALL CMISSField_CreateFinish(GEOMETRIC_FIELD,Err)

    !Update the geometric field parameters
    CALL CMISSGeneratedMesh_GeometricParametersCalculate(GEOMETRIC_FIELD,GENERATED_MESH,Err)

    !Create the equations_set
    CALL CMISSEquationsSet_Initialise(EQUATIONS_SET,Err)
    CALL CMISSField_Initialise(EquationsSetField,Err)
    CALL CMISSEquationsSet_CreateStart(1,REGION,GEOMETRIC_FIELD,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMISS_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE,EquationsSetFieldUserNumber, &
      & EquationsSetField, &
    & EQUATIONS_SET,Err)
    !Set the equations set to be a standard Laplace problem
    
    !Finish creating the equations set
    CALL CMISSEquationsSet_CreateFinish(EQUATIONS_SET,Err)
  
    !Create the equations set dependent field variables
    CALL CMISSField_Initialise(DEPENDENT_FIELD,Err)
    CALL CMISSEquationsSet_DependentCreateStart(EQUATIONS_SET,2,DEPENDENT_FIELD,Err)
    !Finish the equations set dependent field variables
    CALL CMISSEquationsSet_DependentCreateFinish(EQUATIONS_SET,Err)

    !Create the equations set analytic field variables
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      ANALYTIC_FUNCTION=CMISS_EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2
    ELSE
      ANALYTIC_FUNCTION=CMISS_EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2
    ENDIF
    CALL CMISSField_Initialise(ANALYTIC_FIELD,Err)
    CALL CMISSEquationsSet_AnalyticCreateStart(EQUATIONS_SET,ANALYTIC_FUNCTION,3,ANALYTIC_FIELD,Err)
    !Finish the equations set analtyic field variables
    CALL CMISSEquationsSet_AnalyticCreateFinish(EQUATIONS_SET,Err)

    !Create the equations set equations
    CALL CMISSEquations_Initialise(EQUATIONS,Err)
    CALL CMISSEquationsSet_EquationsCreateStart(EQUATIONS_SET,EQUATIONS,Err)
    !Set the equations matrices sparsity type
    CALL CMISSEquations_SparsityTypeSet(EQUATIONS,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
    CALL CMISSEquationsSet_EquationsCreateFinish(EQUATIONS_SET,Err)

    !Create the problem
    CALL CMISSProblem_Initialise(PROBLEM,Err)
    CALL CMISSProblem_CreateStart(1,PROBLEM,Err)
    !Set the problem to be a standard Laplace problem
    CALL CMISSProblem_SpecificationSet(PROBLEM,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS,CMISS_PROBLEM_LAPLACE_EQUATION_TYPE, &
      & CMISS_PROBLEM_STANDARD_LAPLACE_SUBTYPE,Err)
    !Finish creating the problem
    CALL CMISSProblem_CreateFinish(PROBLEM,Err)

    !Create the problem control loop
    CALL CMISSProblem_ControlLoopCreateStart(PROBLEM,Err)
    !Finish creating the problem control
    CALL CMISSProblem_ControlLoopCreateFinish(PROBLEM,Err)

    !Start the creation of the problem solvers
    CALL CMISSSolver_Initialise(Solver,Err)
    CALL CMISSProblem_SolversCreateStart(Problem,Err)
    CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
    !Set solver to direct type
    CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-12_CMISSDP,Err)
    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(Solver,1.0E-12_CMISSDP,Err)
    !CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    !CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_MUMPS_LIBRARY,Err)
    !Finish the creation of the problem solver
    CALL CMISSProblem_SolversCreateFinish(Problem,Err)

    !Start the creation of the problem solver equations
    CALL CMISSSolver_Initialise(Solver,Err)
    CALL CMISSSolverEquations_Initialise(Solver_Equations,Err)
    CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
    !Get the solve equations
    CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
    CALL CMISSSolver_SolverEquationsGet(Solver,Solver_Equations,Err)
    !Set the solver equations sparsity
    CALL CMISSSolverEquations_SparsityTypeSet(Solver_Equations,CMISS_SOLVER_SPARSE_MATRICES,Err)
    !CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_FULL_MATRICES,Err)
    !Add in the equations set
    CALL CMISSSolverEquations_EquationsSetAdd(Solver_Equations,Equations_Set,EquationsSetIndex,Err)
    !Finish the creation of the problem solver equations
    CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

    !Set up the boundary conditions as per the analytic solution
    CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
    CALL CMISSSolverEquations_BoundaryConditionsCreateStart(Solver_Equations,BoundaryConditions,Err)
    CALL CMISSSolverEquations_BoundaryConditionsAnalytic(Solver_Equations,Err)
    CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(Solver_Equations,Err)

    !Solve the problem
    CALL CMISSProblem_Solve(PROBLEM,Err)

  END SUBROUTINE ANALYTICLAPLACE_GENERIC

  SUBROUTINE ANALYTICLAPLACE_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber,GeneratedMeshUserNumber, &
    & ProblemUserNumber)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: CoordinateSystemUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: RegionUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: BasisUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: GeneratedMeshUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: ProblemUserNumber

    CALL CMISSProblem_Destroy(ProblemUserNumber,Err)
    CALL CMISSGeneratedMesh_Destroy(RegionUserNumber,GeneratedMeshUserNumber,Err)
    CALL CMISSBasis_Destroy(BasisUserNumber,Err)
    CALL CMISSRegion_Destroy(RegionUserNumber,Err)
    CALL CMISSCoordinateSystem_Destroy(CoordinateSystemUserNumber,Err)

  END SUBROUTINE ANALYTICLAPLACE_GENERIC_CLEAN

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR

END PROGRAM ANALYTICLAPLACEEXAMPLE 
