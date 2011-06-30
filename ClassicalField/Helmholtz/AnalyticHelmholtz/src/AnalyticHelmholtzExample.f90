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
    CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
    CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the coordinate system to be 2D
      CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,2,Err)
    ELSE
      !Set the coordinate system to be 3D
      CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,3,Err)
    ENDIF
    !Finish the creation of the coordinate system
    CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

    !Start the creation of the region
    CALL CMISSRegionTypeInitialise(Region,Err)
    CALL CMISSRegionCreateStart(RegionUserNumber,WORLD_REGION,Region,Err)
    !Set the regions coordinate system to the 2D RC coordinate system that we have created
    CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
    !Finish the creation of the region
    CALL CMISSRegionCreateFinish(Region,Err)

  
    !Start the creation of a basis (default is trilinear lagrange)
    CALL CMISSBasisTypeInitialise(Basis,Err)
    CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the basis to be a bilinear basis
      CALL CMISSBasisNumberOfXiSet(Basis,2,Err)
      CALL CMISSBasisInterpolationXiSet(Basis,(/INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS/),Err)
    ELSE
      !Set the basis to be a trilinear basis
      CALL CMISSBasisNumberOfXiSet(Basis,3,Err)
      CALL CMISSBasisInterpolationXiSet(Basis,(/INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS, &
          & INTERPOLATION_SPECIFICATIONS/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(Basis,Err)

    !Start the creation of a generated mesh in the region
    CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
    CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
    !Set up a regular mesh
    CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
    CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)
    !Define the mesh on the region
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/WIDTH,HEIGHT/),Err)
      CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/), &
        & Err)
      CALL CMISSGeneratedMeshOriginSet(GeneratedMesh,ORIGIN,Err)
    ELSE
      CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/WIDTH,HEIGHT,LENGTH/),Err)
      CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS, &
        & NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS/),Err)
    ENDIF
    !Finish the creation of a generated mesh in the region
    CALL CMISSMeshTypeInitialise(Mesh,Err)
    CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
    
    !Create a decomposition
    CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
    CALL CMISSDecompositionCreateStart(1,Mesh,Decomposition,Err)
    !Set the decomposition to be a general decomposition with the specified number of domains
    CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
    CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
    CALL CMISSDecompositionCreateFinish(Decomposition,Err)

    !Start to create a default (geometric) field on the region
    CALL CMISSFieldTypeInitialise(GeometricField,Err)
    CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
    !Set the decomposition to use
    CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
    !Set the domain to be used by the field components
    !NB these are needed now as the default mesh component number is 1
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,1,Err)
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,1,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,1,Err)
    ENDIF
    !Finish creating the field
    CALL CMISSFieldCreateFinish(GeometricField,Err)

    !Update the geometric field parameters
    CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

    !Create the equations set
    CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
      CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetClassicalFieldClass, &
    & CMISSEquationsSetHelmholtzEquationType,CMISSEquationsSetStandardHelmholtzSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
    !Set the equations set to be a standard Helmholtz problem
    
    !Finish creating the equations set
    CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)
  
    !Create the equations set dependent field variables
    CALL CMISSFieldTypeInitialise(DEPENDENT_FIELD,Err)
    CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DEPENDENT_FIELD,Err)
    !Finish the equations set dependent field variables
    CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

    !Create the equations set material field variables
    CALL CMISSFieldTypeInitialise(MaterialsField,Err)
    CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
    CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)
    !Set wave number, k
    CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,k,Err)

    !Create the equations set analytic field variables
    AnalyticFunction=CMISSEquationsSetHelmholtzEquationTwoDim1
    CALL CMISSFieldTypeInitialise(AnalyticField,Err)
    CALL CMISSEquationsSetAnalyticCreateStart(EquationsSet,AnalyticFunction,AnalyticFieldUserNumber,AnalyticField,Err)
    !Finish the equations set analtyic field variables
    CALL CMISSEquationsSetAnalyticCreateFinish(EquationsSet,Err)

    !Create the equations set equations
    CALL CMISSEquationsTypeInitialise(Equations,Err)
    CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
    !Set the equations matrices sparsity type
    CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
    CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)
  
    !Create the problem
    CALL CMISSProblemTypeInitialise(Problem,Err)
    CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
    !Set the problem to be a standard Helmholtzproblem
    CALL CMISSProblemSpecificationSet(Problem,CMISSProblemClassicalFieldClass,CMISSProblemHelmholtzEquationType, &
      & CMISSProblemStandardHelmholtzSubtype,Err)
    !Finish creating the problem
    CALL CMISSProblemCreateFinish(Problem,Err)

    !Create the problem control loop
    CALL CMISSProblemControlLoopCreateStart(Problem,Err)
    !Finish creating the problem control
    CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

    !Start the creation of the problem solvers
    CALL CMISSSolverTypeInitialise(Solver,Err)
    CALL CMISSProblemSolversCreateStart(Problem,Err)
    CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
    !Finish the creation of the problem solver
    CALL CMISSProblemSolversCreateFinish(Problem,Err)

    !Start the creation of the problem solver equations
    CALL CMISSSolverTypeInitialise(Solver,Err)
    CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
    CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
    !Get the solve equations
    CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
    CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
    !Set the solver equations sparsity
    CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
    !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)
    !Add in the equations set
    CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
    !Finish the creation of the problem solver equations
    CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

    !Set up the boundary conditions as per the analytic solution
    CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
    CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
    CALL CMISSProblemSolverEquationsBoundaryConditionsAnalytic(SolverEquations,Err)
    CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

    !Solve the problem
    CALL CMISSProblemSolve(Problem,Err)

  END SUBROUTINE ANALYTICHELMHOLTZ_GENERIC

  SUBROUTINE ANALYTICHELMHOLTZ_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber,GeneratedMeshUserNumber, &
    & ProblemUserNumber)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: CoordinateSystemUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: RegionUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: BasisUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: GeneratedMeshUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: ProblemUserNumber

    CALL CMISSProblemDestroy(ProblemUserNumber,Err)
    CALL CMISSGeneratedMeshDestroy(RegionUserNumber,GeneratedMeshUserNumber,Err)
    CALL CMISSBasisDestroy(BasisUserNumber,Err)
    CALL CMISSRegionDestroy(RegionUserNumber,Err)
    CALL CMISSCoordinateSystemDestroy(CoordinateSystemUserNumber,Err)

  END SUBROUTINE ANALYTICHELMHOLTZ_GENERIC_CLEAN



END PROGRAM ANALYTICHELMHOLTZEXAMPLE
