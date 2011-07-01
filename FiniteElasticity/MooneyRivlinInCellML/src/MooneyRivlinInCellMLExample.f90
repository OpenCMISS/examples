!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a finite elasticity equation using openCMISS calls.
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
!> The Original Code is openCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): Jack Lee, Thomas Heidlauf
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

!> \example FiniteElasticity/MooneyRivlinInCellML/src/MooneyRivlinInCellMLExample.f90
!! Example program to solve a finite elasticity equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/MooneyRivlinInCellML/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/MooneyRivlinInCellML/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM MOONEYRIVLININCELLMLEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT!,Filename
!  INTEGER(CMISSIntg), PARAMETER :: NumberOfElementsInEachDirection=10
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldAutoCreate=1 ! 1=yes   0=no
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=1.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponents=2

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariables=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponentsDisplPress=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponentsStressStrain=6

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussPoints=3

  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=12
! State field is not required!
!  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=15

  !Program types


  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx
  INTEGER(CMISSIntg),ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg) :: BottomNormalXi,LeftNormalXi,RightNormalXi,FrontNormalXi

  INTEGER(CMISSIntg) :: DependentVariableTypes(4)
  
  INTEGER(CMISSIntg) :: MooneyRivlinModelIndex
  INTEGER(CMISSIntg) :: CellMLIndex
  

  !CMISS variables

  TYPE(CMISSBasisType) :: QuadraticBasis, LinearBasis, Basis(2)
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,FibreField,DependentField,EquationsSetField
  TYPE(CMISSFieldType) :: MaterialField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSControlLoopType) :: ControlLoop

  TYPE(CMISSCellMLType) :: CellML
  TYPE(CMISSCellMLEquationsType) :: CellMLEquations
  TYPE(CMISSFieldType) :: CellMLModelsField
  TYPE(CMISSFieldType) :: CellMLIntermediateField,CellMLParametersField !,CellMLStateField !SourceField 
  TYPE(CMISSSolverType) :: CellMLSolver

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Get the command arguments
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS == 3) THEN
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NumberGlobalXElements
    IF(NumberGlobalXElements<=0) CALL HANDLE_ERROR("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NumberGlobalYElements
    IF(NumberGlobalYElements<=0) CALL HANDLE_ERROR("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NumberGlobalZElements
    IF(NumberGlobalZElements<0) CALL HANDLE_ERROR("Invalid number of Z elements.")
  ELSE
    NumberGlobalXElements=1
    NumberGlobalYElements=1
    NumberGlobalZElements=1
  ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Intialise cmiss
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
!  CALL CMISSDiagnosticsSetOn(CMISSFromDiagType,[1,2,3,4,5],"Diagnostics",["PROBLEM_FINITE_ELEMENT_CALCULATE"],Err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NumberOfDomains=NumberOfComputationalNodes

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemTypeSet(CoordinateSystem,CMISSCoordinateRectangularCartesianType,Err)
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystemOriginSet(CoordinateSystem,[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP],Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Create a region and assign the CS to the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !Define basis functions - tri-Quadratic Lagrange and tri-Linear Lagrange
  CALL CMISSBasisTypeInitialise(QuadraticBasis,Err)
  CALL CMISSBasisCreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasisTypeSet(QuadraticBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(QuadraticBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(QuadraticBasis,[CMISSBasisQuadraticLagrangeInterpolation, &
    & CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation],Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL CMISSBasisCreateFinish(QuadraticBasis,Err)

  CALL CMISSBasisTypeInitialise(LinearBasis,Err)
  CALL CMISSBasisCreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasisTypeSet(LinearBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(LinearBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(LinearBasis,[CMISSBasisLinearLagrangeInterpolation, &
    & CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation],Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(LinearBasis, &
    & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL CMISSBasisCreateFinish(LinearBasis,Err)

  Basis(1)=QuadraticBasis
  Basis(2)=LinearBasis

  !Start the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the basis - Basis is now an array: [QuadraticBasis,LinearBasis]
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)
  !Define the mesh on the region
  CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[LENGTH,WIDTH,HEIGHT],Err)
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
      & NumberGlobalZElements],Err)
  !Finish the creation of the generated mesh in the region
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecompositionCalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry) - quadratic interpolation
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)  
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

  !Create a fibre field and attach it to the geometric field - linear interpolation
  CALL CMISSFieldTypeInitialise(FibreField,Err)
  CALL CMISSFieldCreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSFieldTypeSet(FibreField,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreField,CMISSFieldUVariableType,FieldFibreNumberOfComponents,Err)  
!LinearMeshComponentNumber QuadraticMeshComponentNumber ???
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,1,LinearMeshComponentNumber,Err) 
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,2,LinearMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,3,LinearMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(FibreField,Err)

  !Create a material field and attach it to the geometric field - quadratic interpolation
  CALL CMISSFieldTypeInitialise(MaterialField,Err)
  CALL CMISSFieldCreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL CMISSFieldTypeSet(MaterialField,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialField,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(MaterialField,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialField,CMISSFieldUVariableType,FieldMaterialNumberOfComponents,Err)
  ! default is CMISSFieldNodeBasedInterpolation
  CALL CMISSFieldComponentInterpolationSet(MaterialField,CMISSFieldUVariableType,1,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(MaterialField,CMISSFieldUVariableType,2,CMISSFieldConstantInterpolation,Err)
! not implemented yet!!!
!  CALL CMISSFieldComponentInterpolationSet(MaterialField,CMISSFieldUVariableType,1,CMISSFieldGaussPointBasedInterpolation,Err)
!  CALL CMISSFieldComponentInterpolationSet(MaterialField,CMISSFieldUVariableType,2,CMISSFieldGaussPointBasedInterpolation,Err)
  CALL CMISSFieldCreateFinish(MaterialField,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,6.0_CMISSDP,Err)

  !Create the dependent field with 4 variables and the respective number of components
  !   1   U_Var_Type            4 components: 3 displacement (quad interpol) + 1 pressure (lin interpol))
  !   2   DELUDELN_Var_Type     4 components: 3 displacement (quad interpol) + 1 pressure (lin interpol))
  !   3   U1_Var_Type           6 components: 6 independent components of the strain tensor (quad interpol) [independent]
  !   4   U2_Var_Type           6 components: 6 independent components of the stress tensor (quad interpol) [dependent]
  DependentVariableTypes = [CMISSFieldUVariableType,CMISSFieldDelUDelNVariableType, &
    & CMISSFieldU1VariableType,CMISSFieldU2VariableType]
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  IF(DependentFieldAutoCreate/=1) THEN
    CALL CMISSFieldCreateStart(FieldDependentUserNumber,Region,DependentField,Err)
    CALL CMISSFieldTypeSet(DependentField,CMISSFieldGeneralType,Err)
    CALL CMISSFieldMeshDecompositionSet(DependentField,Decomposition,Err)
    CALL CMISSFieldGeometricFieldSet(DependentField,GeometricField,Err)
    CALL CMISSFieldDependentTypeSet(DependentField,CMISSFieldDependentType,Err)
    CALL CMISSFieldNumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
    CALL CMISSFieldVariableTypesSet(DependentField,DependentVariableTypes,Err)
    CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldUVariableType,FieldDependentNumberOfComponentsDisplPress,Err)
    CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldDelUDelNVariableType, &
      & FieldDependentNumberOfComponentsDisplPress,Err)
    CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldU1VariableType,FieldDependentNumberOfComponentsStressStrain,Err)
    CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldU2VariableType,FieldDependentNumberOfComponentsStressStrain,Err)

    ! set interpolation for the components of the field variables. 
    !   Default is Node Based Interpolation - so there is nothing to be done for U_Variable_Type and DelUDelN_Variable_Type
    CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldU1VariableType,1,CMISSFieldGaussPointBasedInterpolation,Err)
    CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldU1VariableType,2,CMISSFieldGaussPointBasedInterpolation,Err)
    CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldU1VariableType,3,CMISSFieldGaussPointBasedInterpolation,Err)
    CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldU1VariableType,4,CMISSFieldGaussPointBasedInterpolation,Err)
    CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldU1VariableType,5,CMISSFieldGaussPointBasedInterpolation,Err)
    CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldU1VariableType,6,CMISSFieldGaussPointBasedInterpolation,Err)

    CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldU2VariableType,1,CMISSFieldGaussPointBasedInterpolation,Err)
    CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldU2VariableType,2,CMISSFieldGaussPointBasedInterpolation,Err)
    CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldU2VariableType,3,CMISSFieldGaussPointBasedInterpolation,Err)
    CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldU2VariableType,4,CMISSFieldGaussPointBasedInterpolation,Err)
    CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldU2VariableType,5,CMISSFieldGaussPointBasedInterpolation,Err)
    CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldU2VariableType,6,CMISSFieldGaussPointBasedInterpolation,Err)

    ! set the corresponding mesh component
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,3,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,4,LinearMeshComponentNumber,Err)

    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,2,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,3,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,4,LinearMeshComponentNumber,Err)

    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldU1VariableType,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldU1VariableType,2,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldU1VariableType,3,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldU1VariableType,4,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldU1VariableType,5,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldU1VariableType,6,QuadraticMeshComponentNumber,Err)

    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldU2VariableType,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldU2VariableType,2,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldU2VariableType,3,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldU2VariableType,4,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldU2VariableType,5,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldU2VariableType,6,QuadraticMeshComponentNumber,Err)

!  CALL CMISSFieldScalingTypeSet(DependentField,CMISSFieldUnitScaling,Err)
    CALL CMISSFieldCreateFinish(DependentField,Err)
  ENDIF !DependentFieldAutoCreate

  !Create the equations_set
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSEquationsSetCreateStart(EquationSetUserNumber,Region,FibreField,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetConstitutiveLawInCellMLEvaluateSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  !Create the CellML environment
  CALL CMISSCellMLTypeInitialise(CellML,Err)
  CALL CMISSCellMLCreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import a Mooney-Rivlin material law from a file
  CALL CMISSCellMLModelImport(CellML,"mooney_rivlin.xml",MooneyRivlinModelIndex,Err)
!  CALL CMISSCellMLModelImport(CellML,"n98.xml",MooneyRivlinModelIndex,Err)
  ! Now we have imported the model we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL CMISSCellMLVariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/E11",Err)
  CALL CMISSCellMLVariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/E12",Err)
  CALL CMISSCellMLVariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/E13",Err)
  CALL CMISSCellMLVariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/E22",Err)
  CALL CMISSCellMLVariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/E23",Err)
  CALL CMISSCellMLVariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/E33",Err)
! doesn't work yet, since gauss_point_based_interpolation is not implemented TODO
!  CALL CMISSCellMLVariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/c1",Err)
!  CALL CMISSCellMLVariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/c2",Err)
  !   - to get from the CellML side
  CALL CMISSCellMLVariableSetAsWanted(CellML,MooneyRivlinModelIndex,"equations/Tdev11",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,MooneyRivlinModelIndex,"equations/Tdev12",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,MooneyRivlinModelIndex,"equations/Tdev13",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,MooneyRivlinModelIndex,"equations/Tdev22",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,MooneyRivlinModelIndex,"equations/Tdev23",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,MooneyRivlinModelIndex,"equations/Tdev33",Err)
  !Finish the CellML environment
  CALL CMISSCellMLCreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellMLFieldMapsCreateStart(CellML,Err)
  !Now we can set up the field variable component <--> CellML model variable mappings.
  !Map the strain components
  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentField,CMISSFieldU1VariableType,1,CMISSFieldValuesSetType, &
    & MooneyRivlinModelIndex,"equations/E11",CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentField,CMISSFieldU1VariableType,2,CMISSFieldValuesSetType, &
    & MooneyRivlinModelIndex,"equations/E12",CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentField,CMISSFieldU1VariableType,3,CMISSFieldValuesSetType, &
    & MooneyRivlinModelIndex,"equations/E13",CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentField,CMISSFieldU1VariableType,4,CMISSFieldValuesSetType, &
    & MooneyRivlinModelIndex,"equations/E22",CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentField,CMISSFieldU1VariableType,5,CMISSFieldValuesSetType, &
    & MooneyRivlinModelIndex,"equations/E23",CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentField,CMISSFieldU1VariableType,6,CMISSFieldValuesSetType, &
    & MooneyRivlinModelIndex,"equations/E33",CMISSFieldValuesSetType,Err)
  !Map the material parameters
! doesn't work yet, since gauss_point_based_interpolation is not implemented TODO
!  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,MaterialField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType, &
!    & MooneyRivlinModelIndex,"equations/c1",CMISSFieldValuesSetType,Err)
!  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,MaterialField,CMISSFieldUVariableType,2,CMISSFieldValuesSetType, &
!    & MooneyRivlinModelIndex,"equations/c2",CMISSFieldValuesSetType,Err)
  !Map the stress components
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,MooneyRivlinModelIndex,"equations/Tdev11",CMISSFieldValuesSetType, &
    & DependentField,CMISSFieldU2VariableType,1,CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,MooneyRivlinModelIndex,"equations/Tdev12",CMISSFieldValuesSetType, &
    & DependentField,CMISSFieldU2VariableType,2,CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,MooneyRivlinModelIndex,"equations/Tdev13",CMISSFieldValuesSetType, &
    & DependentField,CMISSFieldU2VariableType,3,CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,MooneyRivlinModelIndex,"equations/Tdev22",CMISSFieldValuesSetType, &
    & DependentField,CMISSFieldU2VariableType,4,CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,MooneyRivlinModelIndex,"equations/Tdev23",CMISSFieldValuesSetType, &
    & DependentField,CMISSFieldU2VariableType,5,CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,MooneyRivlinModelIndex,"equations/Tdev33",CMISSFieldValuesSetType, &
    & DependentField,CMISSFieldU2VariableType,6,CMISSFieldValuesSetType,Err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellMLFieldMapsCreateFinish(CellML,Err)

!Actually, don't need to create this at all - they're automatically created
! CALL CMISSFieldParameterSetCreate(DependentField,CMISSFieldDelUDelNVariableType,CMISSFieldPressureValuesSetType,Err)

  !Create the CellML models field
  CALL CMISSFieldTypeInitialise(CellMLModelsField,Err)
  CALL CMISSCellMLModelsFieldCreateStart(CellMLModelsFieldUserNumber,CellML,CellMLModelsField,Err)
  CALL CMISSCellMLModelsFieldCreateFinish(CellML,Err)

! no state field is required for a simple constitutive equation evaluation
  !Create the CellML state field
!  CALL CMISSFieldTypeInitialise(CellMLStateField,Err)
!  CALL CMISSCellMLStateFieldCreateStart(CellMLStateFieldUserNumber,CellML,CellMLStateField,Err)
!  CALL CMISSCellMLStateFieldCreateFinish(CellML,Err)

  !Create the CellML parameters field --- will be the strain field
  CALL CMISSFieldTypeInitialise(CellMLParametersField,Err)
  CALL CMISSCellMLParametersFieldCreateStart(CellMLParametersFieldUserNumber,CellML,CellMLParametersField,Err)
  CALL CMISSCellMLParametersFieldCreateFinish(CellML,Err)

  !Create the CellML intermediate field --- will be the stress field
  CALL CMISSFieldTypeInitialise(CellMLIntermediateField,Err)
  CALL CMISSCellMLIntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber,CellML,CellMLIntermediateField,Err)
  CALL CMISSCellMLIntermediateFieldCreateFinish(CellML,Err)

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)   

  !Initialise dependent field from undeformed geometry and set hydrostatic pressure
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,-8.0_CMISSDP,Err)

  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,1,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,2,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,3,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,4,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,5,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,6,0.0_CMISSDP,Err)

  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU2VariableType,CMISSFieldValuesSetType,1,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU2VariableType,CMISSFieldValuesSetType,2,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU2VariableType,CMISSFieldValuesSetType,3,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU2VariableType,CMISSFieldValuesSetType,4,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU2VariableType,CMISSFieldValuesSetType,5,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU2VariableType,CMISSFieldValuesSetType,6,0.0_CMISSDP,Err)

  !Define the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemElasticityClass,CMISSProblemFiniteElasticityType, &
    & CMISSProblemFiniteElasticityCellMLSubtype,Err)
   CALL CMISSProblemCreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  CALL CMISSControlLoopTypeSet(ControlLoop,CMISSProblemControlSimpleType,Err)
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)
  
  !Create the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(Solver,CMISSSolverNewtonJacobianFDCalculated,Err)
  CALL CMISSSolverNewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearDirectSolveType,Err)
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Create the problem solver CellML equations
  CALL CMISSSolverTypeInitialise(CellMLSolver,Err)
  CALL CMISSCellMLEquationsTypeInitialise(CellMLEquations,Err)
  CALL CMISSProblemCellMLEquationsCreateStart(Problem,Err)
  CALL CMISSSolverNewtonCellMLSolverGet(Solver,CellMLSolver,Err)
  CALL CMISSSolverCellMLEquationsGet(CellMLSolver,CellMLEquations,Err)
  CALL CMISSCellMLEquationsCellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
  CALL CMISSProblemCellMLEquationsCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)   
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularBottomSurface,BottomSurfaceNodes,BottomNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularLeftSurface,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularRightSurface,RightSurfaceNodes,RightNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularFrontSurface,FrontSurfaceNodes,FrontNormalXi,Err)

  !Set x=0 nodes to no x displacment in x
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NodeNumber,1, &
        & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  !Set x=WIDTH nodes to 10% x displacement
  DO node_idx=1,SIZE(RightSurfaceNodes,1)
    NodeNumber=RightSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NodeNumber,1, &
        & CMISSBoundaryConditionFixed,1.1_CMISSDP*WIDTH,Err)
    ENDIF
  ENDDO

  !Set y=0 nodes to no y displacement
  DO node_idx=1,SIZE(FrontSurfaceNodes,1)
    NodeNumber=FrontSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NodeNumber,2, &
        & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
    ENDIF
  ENDDO

  !Set z=0 nodes to no z displacement
  DO node_idx=1,SIZE(BottomSurfaceNodes,1)
    NodeNumber=BottomSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NodeNumber,3, &
        & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
    ENDIF
  ENDDO

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblemSolve(Problem,Err)

  !Output solution  
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"MooneyRivlinInCellML","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"MooneyRivlinInCellML","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END PROGRAM MOONEYRIVLININCELLMLEXAMPLE

