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

  REAL(CMISSDP), PARAMETER :: ORIGIN(2)=(/-3.141592653579_CMISSDP/2, -3.141592653579_CMISSDP/2/)
  REAL(CMISSDP), PARAMETER :: HEIGHT=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: k=1.0_CMISSDP

  !Program types

  !Program variables

  TYPE(CMISSRegionType) :: WORLD_REGION
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  
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

  CALL ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(2,10,2)
  CALL ANALYTICHELMHOLTZ_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(2,10,2)
  CALL ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_EXPORT(2,6,0)

  CALL CMISSFinalise(Err)

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
    TYPE(CMISSFieldType) :: FIELD

    CALL ANALYTICHELMHOLTZ_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,1, &
      & FIELD)

    CALL CMISSAnalyticAnalysisOutput(FIELD,"AnalyticHelmholtzBilinear",Err)
    
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
    REAL(CMISSDP) :: VALUE
    REAL(CMISSDP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    
    CALL ANALYTICHELMHOLTZ_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,1,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)

    CALL TEST_FRAMEWORK_ASSERT_EQUALS(2.0_CMISSDP,VALUE,0.5_CMISSDP,ERR)
    
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
    REAL(CMISSDP) :: VALUE
    REAL(CMISSDP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)

    CALL ANALYTICHELMHOLTZ_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,3,X_VALUES,Y_VALUES)
    
   CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)
   CALL TEST_FRAMEWORK_ASSERT_EQUALS(4.0_CMISSDP,VALUE,1.0_CMISSDP,Err)
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
    REAL(CMISSDP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    !Local Variables
    REAL(CMISSDP) :: VALUE
    
    INTEGER(CMISSIntg) :: i
    TYPE(CMISSFieldType) :: FIELD
    
    ALLOCATE(X_VALUES((NUMBER_OF_ELEMENTS_XI_END-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1),STAT=ERR)
    ALLOCATE(Y_VALUES((NUMBER_OF_ELEMENTS_XI_END-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1),STAT=ERR)

    DO i = NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL
      
      CALL ANALYTICHELMHOLTZ_GENERIC(i,i,0,INTERPOLATION_SPECIFICATIONS,FIELD)
      CALL CMISSAnalyticAnalysisAbsoluteErrorGetNode(FIELD,1,1,1,(i+1)**2/2+1,1,VALUE,Err)

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
    TYPE(CMISSFieldType) :: DEPENDENT_FIELD
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
    INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=2
    INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=3
    INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=4
    INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=1
    INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

    TYPE(CMISSBasisType) :: Basis
    TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
    TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
    TYPE(CMISSMeshType) :: Mesh
    TYPE(CMISSDecompositionType) :: Decomposition
    TYPE(CMISSEquationsType) :: Equations
    TYPE(CMISSEquationsSetType) :: EquationsSet
    TYPE(CMISSFieldType) :: AnalyticField,GeometricField,MaterialsField
    TYPE(CMISSProblemType) :: Problem
    TYPE(CMISSRegionType) :: Region
    TYPE(CMISSSolverType) :: Solver
    TYPE(CMISSSolverEquationsType) :: SolverEquations
    TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
    
    NUMBER_OF_DOMAINS=1

    !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
    CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(INTERPOLATION_SPECIFICATIONS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

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
    CALL CMISSRegion_CreateStart(RegionUserNumber,WORLD_REGION,Region,Err)
    !Set the regions coordinate system to the 2D RC coordinate system that we have created
    CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
    !Finish the creation of the region
    CALL CMISSRegion_CreateFinish(Region,Err)

  
    !Start the creation of a basis (default is trilinear lagrange)
    CALL CMISSBasis_Initialise(Basis,Err)
    CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the basis to be a bilinear basis
      CALL CMISSBasis_NumberOfXiSet(Basis,2,Err)
      CALL CMISSBasis_InterpolationXiSet(Basis,(/INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS/),Err)
    ELSE
      !Set the basis to be a trilinear basis
      CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
      CALL CMISSBasis_InterpolationXiSet(Basis,(/INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS, &
          & INTERPOLATION_SPECIFICATIONS/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasis_CreateFinish(Basis,Err)

    !Start the creation of a generated mesh in the region
    CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
    CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
    !Set up a regular mesh
    CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
    CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)
    !Define the mesh on the region
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/WIDTH,HEIGHT/),Err)
      CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/), &
        & Err)
      CALL CMISSGeneratedMesh_OriginSet(GeneratedMesh,ORIGIN,Err)
    ELSE
      CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/WIDTH,HEIGHT,LENGTH/),Err)
      CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS, &
        & NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS/),Err)
    ENDIF
    !Finish the creation of a generated mesh in the region
    CALL CMISSMesh_Initialise(Mesh,Err)
    CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
    
    !Create a decomposition
    CALL CMISSDecomposition_Initialise(Decomposition,Err)
    CALL CMISSDecomposition_CreateStart(1,Mesh,Decomposition,Err)
    !Set the decomposition to be a general decomposition with the specified number of domains
    CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
    CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

    !Start to create a default (geometric) field on the region
    CALL CMISSField_Initialise(GeometricField,Err)
    CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
    !Set the decomposition to use
    CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
    !Set the domain to be used by the field components
    !NB these are needed now as the default mesh component number is 1
    CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
    CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
    ENDIF
    !Finish creating the field
    CALL CMISSField_CreateFinish(GeometricField,Err)

    !Update the geometric field parameters
    CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

    !Create the equations set
    CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
      CALL CMISSField_Initialise(EquationsSetField,Err)
CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE,CMISS_EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
    !Set the equations set to be a standard Helmholtz problem
    
    !Finish creating the equations set
    CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)
  
    !Create the equations set dependent field variables
    CALL CMISSField_Initialise(DEPENDENT_FIELD,Err)
    CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DEPENDENT_FIELD,Err)
    !Finish the equations set dependent field variables
    CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

    !Create the equations set material field variables
    CALL CMISSField_Initialise(MaterialsField,Err)
    CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
    CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)
    !Set wave number, k
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,k,Err)

    !Create the equations set analytic field variables
    AnalyticFunction=CMISS_EQUATIONS_SET_HELMHOLTZ_EQUATION_TWO_DIM_1
    CALL CMISSField_Initialise(AnalyticField,Err)
    CALL CMISSEquationsSet_AnalyticCreateStart(EquationsSet,AnalyticFunction,AnalyticFieldUserNumber,AnalyticField,Err)
    !Finish the equations set analtyic field variables
    CALL CMISSEquationsSet_AnalyticCreateFinish(EquationsSet,Err)

    !Create the equations set equations
    CALL CMISSEquations_Initialise(Equations,Err)
    CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
    !Set the equations matrices sparsity type
    CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
    CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)
  
    !Create the problem
    CALL CMISSProblem_Initialise(Problem,Err)
    CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
    !Set the problem to be a standard Helmholtzproblem
    CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS,CMISS_PROBLEM_HELMHOLTZ_EQUATION_TYPE, &
      & CMISS_PROBLEM_STANDARD_HELMHOLTZ_SUBTYPE,Err)
    !Finish creating the problem
    CALL CMISSProblem_CreateFinish(Problem,Err)

    !Create the problem control loop
    CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
    !Finish creating the problem control
    CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

    !Start the creation of the problem solvers
    CALL CMISSSolver_Initialise(Solver,Err)
    CALL CMISSProblem_SolversCreateStart(Problem,Err)
    CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
    !Finish the creation of the problem solver
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

  END SUBROUTINE ANALYTICHELMHOLTZ_GENERIC

  SUBROUTINE ANALYTICHELMHOLTZ_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber,GeneratedMeshUserNumber, &
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

  END SUBROUTINE ANALYTICHELMHOLTZ_GENERIC_CLEAN



END PROGRAM ANALYTICHELMHOLTZEXAMPLE
