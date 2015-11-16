!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve an Analytic Helmholtz equation using OpenCMISS calls.
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

!> \example ClassicalField/Helmholtz/AnalyticHelmholtz/src/AnalyticHelmholtzExample.f90
!! Example illustrating the use of OpenCMISS to solve the Helmholtz problem and check with its Analytic Solution.
!! 
!! \htmlinclude ClassicalField/Helmholtz/AnalyticHelmholtz/history.html
!< 

!> Main program
PROGRAM ANALYTICHELMHOLTZEXAMPLE
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

  REAL(CMISSRP), PARAMETER :: ORIGIN(2)=[-3.141592653579_CMISSRP/2, -3.141592653579_CMISSRP/2]
  REAL(CMISSRP), PARAMETER :: HEIGHT=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: k=1.0_CMISSRP

  !Program types

  !Program variables

  TYPE(cmfe_RegionType) :: WORLD_REGION
  TYPE(cmfe_CoordinateSystemType) :: WorldCoordinateSystem
  
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

  CALL ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(2,10,2)
  CALL ANALYTICHELMHOLTZ_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(2,10,2)
  CALL ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_EXPORT(2,6,0)

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

CONTAINS

  !
  !================================================================================================================================
  !  
    !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_EXPORT(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<number of elements in x direction
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<number of elements in y direction
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<number of elements in z direction
    !Local Variables
    TYPE(cmfe_FieldType) :: FIELD

    CALL ANALYTICHELMHOLTZ_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,1, &
      & FIELD)

    CALL cmfe_AnalyticAnalysis_Output(FIELD,"AnalyticHelmholtzBilinear",Err)
    
    CALL ANALYTICHELMHOLTZ_GENERIC_CLEAN(1,1,1,1,1)

  END SUBROUTINE ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_EXPORT
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(CMISSRP) :: VALUE
    REAL(CMISSRP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    
    CALL ANALYTICHELMHOLTZ_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,1,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)

    CALL TEST_FRAMEWORK_ASSERT_EQUALS(2.0_CMISSRP,VALUE,0.5_CMISSRP,ERR)
    
    WRITE(*,'(A)') "Analytic Helmholtz Example Testcase1 - bilinear lagrange is successfully completed."
    
  END SUBROUTINE ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICHELMHOLTZ_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(CMISSRP) :: VALUE
    REAL(CMISSRP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)

    CALL ANALYTICHELMHOLTZ_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,3,X_VALUES,Y_VALUES)
    
   CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)
   CALL TEST_FRAMEWORK_ASSERT_EQUALS(4.0_CMISSRP,VALUE,1.0_CMISSRP,Err)
   IF (Err/=0) THEN
     WRITE(*,'(A,F3.5)') "Analytic Helmholtz Example Testcase2 - bicubic Hermite failure: Convergence should be around 4.0" &
       & //", but it was ", VALUE
   ENDIF
   WRITE(*,'(A)') "Analytic Helmholtz Example Testcase2 - bicubic Hermite is successfully completed."

  END SUBROUTINE ANALYTICHELMHOLTZ_TESTCASE_BICUBIC_HERMITE_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICHELMHOLTZ_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
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
      
      CALL ANALYTICHELMHOLTZ_GENERIC(i,i,0,INTERPOLATION_SPECIFICATIONS,FIELD)
      CALL cmfe_AnalyticAnalysis_AbsoluteErrorGetNode(FIELD,1,1,1,(i+1)**2/2+1,1,VALUE,Err)

      Y_VALUES((i-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1)=log10(VALUE)
      X_VALUES((i-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1)=log10(HEIGHT/i)
      CALL ANALYTICHELMHOLTZ_GENERIC_CLEAN(1,1,1,1,1)
   
    ENDDO
  END SUBROUTINE ANALYTICHELMHOLTZ_GENERIC_CONVERGENCE
  
  
  !
  !================================================================================================================================
  !   
    
  SUBROUTINE ANALYTICHELMHOLTZ_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
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

    INTEGER(CMISSIntg) :: AnalyticFunction
    INTEGER(CMISSIntg) :: EquationsSetIndex

    INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
    INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
    INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=1
    INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
    INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
    INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
    INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=1
    INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=2
    INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=3
    INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=4
    INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=5
    INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=1
    INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

    TYPE(cmfe_BasisType) :: Basis
    TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
    TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
    TYPE(cmfe_MeshType) :: Mesh
    TYPE(cmfe_DecompositionType) :: Decomposition
    TYPE(cmfe_EquationsType) :: Equations
    TYPE(cmfe_EquationsSetType) :: EquationsSet
    TYPE(cmfe_FieldType) :: AnalyticField,GeometricField,EquationsSetField,MaterialsField
    TYPE(cmfe_ProblemType) :: Problem
    TYPE(cmfe_RegionType) :: Region
    TYPE(cmfe_SolverType) :: Solver
    TYPE(cmfe_SolverEquationsType) :: SolverEquations
    TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
    
    NUMBER_OF_DOMAINS=1

    !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
    CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(INTERPOLATION_SPECIFICATIONS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

    !Start the creation of a new RC coordinate system
    CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
    CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
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
    CALL cmfe_Region_Initialise(Region,Err)
    CALL cmfe_Region_CreateStart(RegionUserNumber,WORLD_REGION,Region,Err)
    !Set the regions coordinate system to the 2D RC coordinate system that we have created
    CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
    !Finish the creation of the region
    CALL cmfe_Region_CreateFinish(Region,Err)
  
    !Start the creation of a basis (default is trilinear lagrange)
    CALL cmfe_Basis_Initialise(Basis,Err)
    CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
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
    !Finish the creation of the basis
    CALL cmfe_Basis_CreateFinish(Basis,Err)

    !Start the creation of a generated mesh in the region
    CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
    CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
    !Set up a regular mesh
    CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
    CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)
    !Define the mesh on the region
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS], &
        & Err)
      CALL cmfe_GeneratedMesh_OriginSet(GeneratedMesh,ORIGIN,Err)
    ELSE
      CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS, &
        & NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS],Err)
    ENDIF
    !Finish the creation of a generated mesh in the region
    CALL cmfe_Mesh_Initialise(Mesh,Err)
    CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
    
    !Create a decomposition
    CALL cmfe_Decomposition_Initialise(Decomposition,Err)
    CALL cmfe_Decomposition_CreateStart(1,Mesh,Decomposition,Err)
    !Set the decomposition to be a general decomposition with the specified number of domains
    CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
    CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

    !Start to create a default (geometric) field on the region
    CALL cmfe_Field_Initialise(GeometricField,Err)
    CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
    !Set the decomposition to use
    CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
    !Set the domain to be used by the field components
    !NB these are needed now as the default mesh component number is 1
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
    ENDIF
    !Finish creating the field
    CALL cmfe_Field_CreateFinish(GeometricField,Err)

    !Update the geometric field parameters
    CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

    !Create the equations set
    CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
    CALL cmfe_Field_Initialise(EquationsSetField,Err)
    !Set the equations set to be a standard Helmholtz problem
    CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
      & CMFE_EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE,CMFE_EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE],EquationsSetFieldUserNumber, &
      & EquationsSetField,EquationsSet,Err)
    !Finish creating the equations set
    CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)
  
    !Create the equations set dependent field variables
    CALL cmfe_Field_Initialise(DEPENDENT_FIELD,Err)
    CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DEPENDENT_FIELD,Err)
    !Finish the equations set dependent field variables
    CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

    !Create the equations set material field variables
    CALL cmfe_Field_Initialise(MaterialsField,Err)
    CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
    CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)
    !Set wave number, k
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,k,Err)

    !Create the equations set analytic field variables
    AnalyticFunction=CMFE_EQUATIONS_SET_HELMHOLTZ_EQUATION_TWO_DIM_1
    CALL cmfe_Field_Initialise(AnalyticField,Err)
    CALL cmfe_EquationsSet_AnalyticCreateStart(EquationsSet,AnalyticFunction,AnalyticFieldUserNumber,AnalyticField,Err)
    !Finish the equations set analtyic field variables
    CALL cmfe_EquationsSet_AnalyticCreateFinish(EquationsSet,Err)

    !Create the equations set equations
    CALL cmfe_Equations_Initialise(Equations,Err)
    CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
    !Set the equations matrices sparsity type
    CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
    CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)
  
    !Create the problem
    CALL cmfe_Problem_Initialise(Problem,Err)
    CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS,CMFE_PROBLEM_HELMHOLTZ_EQUATION_TYPE, &
      & CMFE_PROBLEM_STANDARD_HELMHOLTZ_SUBTYPE],Problem,Err)
    !Finish creating the problem
    CALL cmfe_Problem_CreateFinish(Problem,Err)

    !Create the problem control loop
    CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
    !Finish creating the problem control
    CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

    !Start the creation of the problem solvers
    CALL cmfe_Solver_Initialise(Solver,Err)
    CALL cmfe_Problem_SolversCreateStart(Problem,Err)
    CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
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

    !Set up the boundary conditions as per the analytic solution
    CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
    CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
    CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(SolverEquations,Err)
    CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

    !Solve the problem
    CALL cmfe_Problem_Solve(Problem,Err)

  END SUBROUTINE ANALYTICHELMHOLTZ_GENERIC

  SUBROUTINE ANALYTICHELMHOLTZ_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber,GeneratedMeshUserNumber, &
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

  END SUBROUTINE ANALYTICHELMHOLTZ_GENERIC_CLEAN



END PROGRAM ANALYTICHELMHOLTZEXAMPLE
