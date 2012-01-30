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

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
!  CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_FINITE_ELEMENT_CALCULATE"],Err)

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
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_TypeSet(CoordinateSystem,CMISS_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystem_OriginSet(CoordinateSystem,[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP],Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the CS to the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Define basis functions - tri-Quadratic Lagrange and tri-Linear Lagrange
  CALL CMISSBasis_Initialise(QuadraticBasis,Err)
  CALL CMISSBasis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasis_TypeSet(QuadraticBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(QuadraticBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(QuadraticBasis,[CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL CMISSBasis_CreateFinish(QuadraticBasis,Err)

  CALL CMISSBasis_Initialise(LinearBasis,Err)
  CALL CMISSBasis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasis_TypeSet(LinearBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(LinearBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(LinearBasis,[CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearBasis, &
    & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL CMISSBasis_CreateFinish(LinearBasis,Err)

  Basis(1)=QuadraticBasis
  Basis(2)=LinearBasis

  !Start the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the basis - Basis is now an array: [QuadraticBasis,LinearBasis]
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)
  !Define the mesh on the region
  CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH,WIDTH,HEIGHT],Err)
  CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
      & NumberGlobalZElements],Err)
  !Finish the creation of the generated mesh in the region
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry) - quadratic interpolation
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)  
  CALL CMISSField_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Create a fibre field and attach it to the geometric field - linear interpolation
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)  
!LinearMeshComponentNumber QuadraticMeshComponentNumber ???
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,1,LinearMeshComponentNumber,Err) 
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,2,LinearMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,3,LinearMeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(FibreField,Err)

  !Create a material field and attach it to the geometric field - quadratic interpolation
  CALL CMISSField_Initialise(MaterialField,Err)
  CALL CMISSField_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL CMISSField_TypeSet(MaterialField,CMISS_FIELD_MATERIAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(MaterialField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(MaterialField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)
  ! default is CMISS_FIELD_NODE_BASED_INTERPOLATION
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
! not implemented yet!!!
!  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
!  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
  CALL CMISSField_CreateFinish(MaterialField,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,6.0_CMISSDP,Err)

  !Create the dependent field with 4 variables and the respective number of components
  !   1   U_Var_Type            4 components: 3 displacement (quad interpol) + 1 pressure (lin interpol))
  !   2   DELUDELN_Var_Type     4 components: 3 displacement (quad interpol) + 1 pressure (lin interpol))
  !   3   U1_Var_Type           6 components: 6 independent components of the strain tensor (quad interpol) [independent]
  !   4   U2_Var_Type           6 components: 6 independent components of the stress tensor (quad interpol) [dependent]
  DependentVariableTypes = [CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE, &
    & CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_U2_VARIABLE_TYPE]
  CALL CMISSField_Initialise(DependentField,Err)
  IF(DependentFieldAutoCreate/=1) THEN
    CALL CMISSField_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
    CALL CMISSField_TypeSet(DependentField,CMISS_FIELD_GENERAL_TYPE,Err)
    CALL CMISSField_MeshDecompositionSet(DependentField,Decomposition,Err)
    CALL CMISSField_GeometricFieldSet(DependentField,GeometricField,Err)
    CALL CMISSField_DependentTypeSet(DependentField,CMISS_FIELD_DEPENDENT_TYPE,Err)
    CALL CMISSField_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
    CALL CMISSField_VariableTypesSet(DependentField,DependentVariableTypes,Err)
    CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponentsDisplPress,Err)
    CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE, &
      & FieldDependentNumberOfComponentsDisplPress,Err)
    CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE, &
      & FieldDependentNumberOfComponentsStressStrain, &
      & Err)
    CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE, &
      & FieldDependentNumberOfComponentsStressStrain, &
      & Err)

    ! set interpolation for the components of the field variables. 
    !   Default is Node Based Interpolation - so there is nothing to be done for U_Variable_Type and DelUDelN_Variable_Type
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,1, &
      & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,2, &
      & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,3, &
      & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,4, &
      & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,5, &
      & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,6, &
      & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)

    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,1, &
      & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,2, &
      & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,3, &
      & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,4, &
      & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,5, &
      & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,6, &
      & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)

    ! set the corresponding mesh component
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)

    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)

    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,4,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,5,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,6,QuadraticMeshComponentNumber,Err)

    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,4,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,5,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,6,QuadraticMeshComponentNumber,Err)

!  CALL CMISSField_ScalingTypeSet(DependentField,CMISS_FIELD_UNIT_SCALING,Err)
    CALL CMISSField_CreateFinish(DependentField,Err)
  ENDIF !DependentFieldAutoCreate

  !Create the equations_set
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSEquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
      & EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Create the CellML environment
  CALL CMISSCellML_Initialise(CellML,Err)
  CALL CMISSCellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import a Mooney-Rivlin material law from a file
  CALL CMISSCellML_ModelImport(CellML,"mooney_rivlin.xml",MooneyRivlinModelIndex,Err)
!  CALL CMISSCellML_ModelImport(CellML,"n98.xml",MooneyRivlinModelIndex,Err)
  ! Now we have imported the model we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL CMISSCellML_VariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/E11",Err)
  CALL CMISSCellML_VariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/E12",Err)
  CALL CMISSCellML_VariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/E13",Err)
  CALL CMISSCellML_VariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/E22",Err)
  CALL CMISSCellML_VariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/E23",Err)
  CALL CMISSCellML_VariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/E33",Err)
! doesn't work yet, since gauss_point_based_interpolation is not implemented TODO
!  CALL CMISSCellML_VariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/c1",Err)
!  CALL CMISSCellML_VariableSetAsKnown(CellML,MooneyRivlinModelIndex,"equations/c2",Err)
  !   - to get from the CellML side
  CALL CMISSCellML_VariableSetAsWanted(CellML,MooneyRivlinModelIndex,"equations/Tdev11",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,MooneyRivlinModelIndex,"equations/Tdev12",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,MooneyRivlinModelIndex,"equations/Tdev13",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,MooneyRivlinModelIndex,"equations/Tdev22",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,MooneyRivlinModelIndex,"equations/Tdev23",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,MooneyRivlinModelIndex,"equations/Tdev33",Err)
  !Finish the CellML environment
  CALL CMISSCellML_CreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellML_FieldMapsCreateStart(CellML,Err)
  !Now we can set up the field variable component <--> CellML model variable mappings.
  !Map the strain components
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & MooneyRivlinModelIndex,"equations/E11",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,2,CMISS_FIELD_VALUES_SET_TYPE, &
    & MooneyRivlinModelIndex,"equations/E12",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,3,CMISS_FIELD_VALUES_SET_TYPE, &
    & MooneyRivlinModelIndex,"equations/E13",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,4,CMISS_FIELD_VALUES_SET_TYPE, &
    & MooneyRivlinModelIndex,"equations/E22",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,5,CMISS_FIELD_VALUES_SET_TYPE, &
    & MooneyRivlinModelIndex,"equations/E23",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,6,CMISS_FIELD_VALUES_SET_TYPE, &
    & MooneyRivlinModelIndex,"equations/E33",CMISS_FIELD_VALUES_SET_TYPE,Err)
  !Map the material parameters
! doesn't work yet, since gauss_point_based_interpolation is not implemented TODO
!  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
!    & MooneyRivlinModelIndex,"equations/c1",CMISS_FIELD_VALUES_SET_TYPE,Err)
!  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_VALUES_SET_TYPE, &
!    & MooneyRivlinModelIndex,"equations/c2",CMISS_FIELD_VALUES_SET_TYPE,Err)
  !Map the stress components
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,MooneyRivlinModelIndex,"equations/Tdev11",CMISS_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,MooneyRivlinModelIndex,"equations/Tdev12",CMISS_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,2,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,MooneyRivlinModelIndex,"equations/Tdev13",CMISS_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,3,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,MooneyRivlinModelIndex,"equations/Tdev22",CMISS_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,4,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,MooneyRivlinModelIndex,"equations/Tdev23",CMISS_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,5,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,MooneyRivlinModelIndex,"equations/Tdev33",CMISS_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,6,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellML_FieldMapsCreateFinish(CellML,Err)

!Actually, don't need to create this at all - they're automatically created
! CALL CMISSField_ParameterSetCreate(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_PRESSURE_VALUES_SET_TYPE,Err)

  !Create the CellML models field
  CALL CMISSField_Initialise(CellMLModelsField,Err)
  CALL CMISSCellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  CALL CMISSCellML_ModelsFieldCreateFinish(CellML,Err)

! no state field is required for a simple constitutive equation evaluation
  !Create the CellML state field
!  CALL CMISSField_Initialise(CellMLStateField,Err)
!  CALL CMISSCellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,Err)
!  CALL CMISSCellML_StateFieldCreateFinish(CellML,Err)

  !Create the CellML parameters field --- will be the strain field
  CALL CMISSField_Initialise(CellMLParametersField,Err)
  CALL CMISSCellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  CALL CMISSCellML_ParametersFieldCreateFinish(CellML,Err)

  !Create the CellML intermediate field --- will be the stress field
  CALL CMISSField_Initialise(CellMLIntermediateField,Err)
  CALL CMISSCellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  CALL CMISSCellML_IntermediateFieldCreateFinish(CellML,Err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)   

  !Initialise dependent field from undeformed geometry and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,-8.0_CMISSDP, &
    & Err)

  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,5,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,6,0.0_CMISSDP, &
    & Err)

  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,5,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,6,0.0_CMISSDP, &
    & Err)

  !Define the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_ELASTICITY_CLASS,CMISS_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMISS_PROBLEM_FINITE_ELASTICITY_CELLML_SUBTYPE,Err)
   CALL CMISSProblem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL CMISSControlLoop_TypeSet(ControlLoop,CMISS_PROBLEM_CONTROL_SIMPLE_TYPE,Err)
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)
  
  !Create the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Create the problem solver CellML equations
  CALL CMISSSolver_Initialise(CellMLSolver,Err)
  CALL CMISSCellMLEquations_Initialise(CellMLEquations,Err)
  CALL CMISSProblem_CellMLEquationsCreateStart(Problem,Err)
  CALL CMISSSolver_NewtonCellMLSolverGet(Solver,CellMLSolver,Err)
  CALL CMISSSolver_CellMLEquationsGet(CellMLSolver,CellMLEquations,Err)
  CALL CMISSCellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
  CALL CMISSProblem_CellMLEquationsCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)   
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi, &
    & Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)

  !Set x=0 nodes to no x displacment in x
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  !Set x=WIDTH nodes to 10% x displacement
  DO node_idx=1,SIZE(RightSurfaceNodes,1)
    NodeNumber=RightSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED,1.1_CMISSDP*WIDTH,Err)
    ENDIF
  ENDDO

  !Set y=0 nodes to no y displacement
  DO node_idx=1,SIZE(FrontSurfaceNodes,1)
    NodeNumber=FrontSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO

  !Set z=0 nodes to no z displacement
  DO node_idx=1,SIZE(BottomSurfaceNodes,1)
    NodeNumber=BottomSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Output solution  
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"MooneyRivlinInCellML","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"MooneyRivlinInCellML","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

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

