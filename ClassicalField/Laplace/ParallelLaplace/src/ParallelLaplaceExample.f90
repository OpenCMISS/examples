!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a Laplace equation in parallel using OpenCMISS calls.
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

!> \example ClassicalField/Laplace/ParallelLaplace/src/ParallelLaplaceExample.f90
!! Example program to solve a Laplace equation in parallel using OpenCMISS calls.
!! \htmlinclude ClassicalField/Laplace/ParallelLaplace/history.html
!!
!<

!> Main program
PROGRAM PARALLELLAPLACEEXAMPLE

  USE OPENCMISS

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
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=12
 
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS,INTERPOLATION_TYPE
  
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT,Filename

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,EquationsSetField,AnalyticField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
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

  !Get the command arguments
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS == 4) THEN
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_X_ELEMENTS
    IF(NUMBER_GLOBAL_X_ELEMENTS<=0) CALL HANDLE_ERROR("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Y_ELEMENTS
    IF(NUMBER_GLOBAL_Y_ELEMENTS<=0) CALL HANDLE_ERROR("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Z_ELEMENTS
    IF(NUMBER_GLOBAL_Y_ELEMENTS<0) CALL HANDLE_ERROR("Invalid number of Z elements.")
    CALL GET_COMMAND_ARGUMENT(4,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 4.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) INTERPOLATION_TYPE
    IF(INTERPOLATION_TYPE<=0.OR.INTERPOLATION_TYPE>9) CALL HANDLE_ERROR("Invalid Interpolation specification.")
  ELSE
    NUMBER_GLOBAL_X_ELEMENTS=25
    NUMBER_GLOBAL_Y_ELEMENTS=25
    NUMBER_GLOBAL_Z_ELEMENTS=25
    INTERPOLATION_TYPE=1
  ENDIF
  
  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  CALL CMISSDiagnosticsSetOn(CMISSInDiagType,[1,2,3,4,5],"Diagnostics",["DISTRIBUTED_VECTOR_UPDATE_START"],Err)

  WRITE(Filename,'(A,"_",I0,"_",I0,"_",I0,"_",I0)') "Laplace",NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS,INTERPOLATION_TYPE
  
  !CALL CMISSOutputSetOn(Filename,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  
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
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasisTypeInitialise(Basis,Err)
  CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
  IF(INTERPOLATION_TYPE>6) CALL CMISSBasisTypeSet(Basis,CMISSBasisSimplexType,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSBasisNumberOfXiSet(Basis,2,Err)
    CALL CMISSBasisInterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
  ELSE
    CALL CMISSBasisNumberOfXiSet(Basis,3,Err)
    CALL CMISSBasisInterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
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
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
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
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)
  
  !Create the equations_set
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
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Create the equations set analytic field variables
  CALL CMISSFieldTypeInitialise(AnalyticField,Err)
  !Create the equations set analytic 
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSEquationsSetAnalyticCreateStart(EquationsSet,CMISSEquationsSetLaplaceEquationThreeDim2,AnalyticFieldUserNumber, &
      AnalyticField,Err)
  ELSE
    CALL CMISSEquationsSetAnalyticCreateStart(EquationsSet,CMISSEquationsSetLaplaceEquationTwoDim2,AnalyticFieldUserNumber, &
      AnalyticField,Err)
  ENDIF
  !Finish the equations set analytic
  CALL CMISSEquationsSetAnalyticCreateFinish(EquationsSet,Err)
  
  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
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
  CALL CMISSSolverLinearTypeSet(Solver,CMISSSolverLinearDirectSolveType,Err)
  CALL CMISSSolverLibraryTypeSet(Solver,CMISSSolverMUMPSLibrary,Err)
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

  !Compute Analytic analysis
  !CALL CMISSAnalyticAnalysisOutput(DependentField,"ParallelLaplace",Err)
  
  !Export result
  !CALL CMISSFieldsTypeInitialise(Fields,Err)
  !CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  !CALL CMISSFieldIONodesExport(Fields,"ParallelLaplace","FORTRAN",Err)
  !CALL CMISSFieldIOElementsExport(Fields,"ParallelLaplace","FORTRAN",Err)
  !CALL CMISSFieldsTypeFinalise(Fields,Err)
  
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
  
END PROGRAM PARALLELLAPLACEEXAMPLE

