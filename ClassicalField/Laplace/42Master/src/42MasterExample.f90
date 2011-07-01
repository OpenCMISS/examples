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
  USE MPI


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=3.0_CMISSDP

  LOGICAL, PARAMETER :: SOLVER_DIRECT_TYPE=.TRUE.

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11

  !Program types

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS = -1
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_Y_ELEMENTS = -1
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_Z_ELEMENTS = -1
  INTEGER(CMISSIntg) :: INTERPOLATION_TYPE
  INTEGER(CMISSIntg) :: NUMBER_OF_GAUSS_XI
  INTEGER(CMISSIntg) :: I
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

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,EquationsSetField,DependentField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables

  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain
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
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  WRITE(Filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "Laplace",NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS,INTERPOLATION_TYPE

  CALL CMISSOutputSetOn(Filename,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  IF(OPTION_1D) THEN
    !Set the coordinate system to be 1D
    CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,1,Err)
  ELSEIF(OPTION_2D) THEN
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
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionLabelSet(Region,"LaplaceRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasisTypeInitialise(Basis,Err)
  CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
  IF(OPTION_TRI .OR. OPTION_TET) THEN !Default is Lagrange Hermite type
    CALL CMISSBasisTypeSet(Basis,CMISSBasisSimplexType,Err)
  ENDIF
  IF(OPTION_1D) THEN
    !Set the basis to be a linear Lagrange basis
    CALL CMISSBasisNumberOfXiSet(Basis,1,Err)
    CALL CMISSBasisInterpolationXiSet(Basis,[INTERPOLATION_TYPE],Err)
    IF(.NOT. OPTION_TRI .AND. .NOT. OPTION_TET) THEN
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ELSEIF(OPTION_2D) THEN
    !Set the basis to be a bilinear Lagrange basis
    CALL CMISSBasisNumberOfXiSet(Basis,2,Err)
    CALL CMISSBasisInterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(.NOT. OPTION_TRI .AND. .NOT. OPTION_TET) THEN
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL CMISSBasisNumberOfXiSet(Basis,3,Err)
    CALL CMISSBasisInterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(.NOT. OPTION_TRI .AND. .NOT. OPTION_TET) THEN
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(Basis,Err)

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)
  !Define the mesh on the region
  IF(OPTION_1D) THEN
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[WIDTH],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS],Err)
  ELSEIF(OPTION_2D) THEN
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,1,Err)
  IF(NUMBER_GLOBAL_Y_ELEMENTS/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,1,Err)
  ENDIF
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

  !Create the Standard Laplace Equations set
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetClassicalFieldClass, &
    & CMISSEquationsSetLaplaceEquationType,CMISSEquationsSetStandardLaplaceSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Set the DOFs to be contiguous across components
  CALL CMISSFieldDOFOrderTypeSet(DependentField,CMISSFieldUVariableType,CMISSFieldContiguousComponentDOFOrder,Err)
  CALL CMISSFieldDOFOrderTypeSet(DependentField,CMISSFieldDelUDelNVariableType,CMISSFieldContiguousComponentDOFOrder,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Initialise the field with an initial guess
  IF(.NOT. SOLVER_DIRECT_TYPE) THEN
    CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,0.5_CMISSDP,Err)
  ENDIF

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  !CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsFullMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemClassicalFieldClass,CMISSProblemLaplaceEquationType, &
    & CMISSProblemStandardLaplaceSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
  IF(SOLVER_DIRECT_TYPE) THEN
    CALL CMISSSolverLinearTypeSet(Solver,CMISSSolverLinearDirectSolveType,Err)
    CALL CMISSSolverLibraryTypeSet(Solver,CMISSSolverMUMPSLibrary,Err)
    !CALL CMISSSolverLibraryTypeSet(Solver,CMISSSolverSuperLULibrary,Err)
    !CALL CMISSSolverLibraryTypeSet(Solver,CMISSSolverPaStiXLibrary,Err)
  ELSE
    CALL CMISSSolverLinearTypeSet(Solver,CMISSSolverLinearIterativeSolveType,Err)
    CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(Solver,1.0E-10_CMISSDP,Err)
    CALL CMISSSolverLinearIterativeRelativeToleranceSet(Solver,1.0E-10_CMISSDP,Err)
    CALL CMISSSolverLinearIterativeMaximumIterationsSet(Solver,100000,Err)
  ENDIF
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

  !Start the creation of the equations set boundary conditions
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSRegionNodesGet(Region,Nodes,Err)
  CALL CMISSNodesNumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL CMISSDecompositionNodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL CMISSDecompositionNodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,FirstNodeNumber,1, &
      & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,LastNodeNumber,1, &
      & CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)

  !Export results
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"Laplace","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"Laplace","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  !Finialise CMISS
  CALL CMISSFinalise(Err)

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
