!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a coupled Monodomain equation Finite Elasticity problem using OpenCMISS calls.
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
!> Contributor(s): Thomas Heidlauf
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

!> \example MultiPhysics/BioelectricFiniteElasticity/GudunovMonodomainElasticitySameMesh/src/GudunovMonodomainElasticitySameMeshExample.f90
!! Example program to solve a Monodomain equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/BioelectricFiniteElasticity/GudunovMonodomainElasticitySameMesh/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/BioelectricFiniteElasticity/GudunovMonodomainElasticitySameMesh/build-gnu'>Linux GNU Build</a>
!!
!<

!> Main program
PROGRAM GODUNOVMONODOMAINELASTICITYSAMEMESHEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !--------------------------------------------------------------------------------------------------------------------------------
  !Test program parameters
  REAL(CMISSDP), PARAMETER :: tol=1.0E-8_CMISSDP

  LOGICAL :: independent_field_auto_create=.FALSE.
  LOGICAL :: uniaxial_extension_bc=.FALSE.

  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_ELEMENTS=50
  !all lengths in [cm]
  REAL(CMISSDP), PARAMETER :: LENGTH= 2.5_CMISSDP ! X-direction
  REAL(CMISSDP), PARAMETER :: WIDTH= 0.05_CMISSDP ! Y-direction
  REAL(CMISSDP), PARAMETER :: HEIGHT=0.05_CMISSDP ! Z-direction
  
  !all times in [ms]
  REAL(CMISSDP), PARAMETER :: STIM_STOP=0.05_CMISSDP
  REAL(CMISSDP), PARAMETER :: TIME_STOP=50.0_CMISSDP

  REAL(CMISSDP), PARAMETER :: ODE_TIME_STEP=0.00001_CMISSDP
  REAL(CMISSDP), PARAMETER :: PDE_TIME_STEP=0.0005_CMISSDP
  REAL(CMISSDP), PARAMETER :: ELASTICITY_TIME_STEP=0.01_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: OUTPUT_FREQUENCY=2

  !--------------------------------------------------------------------------------------------------------------------------------
  !stimulation current in [uA/cm^2]
  REAL(CMISSDP), PARAMETER :: STIM_VALUE=8000.0_CMISSDP
  !condctivity in [mS/cm]
  REAL(CMISSDP), PARAMETER :: CONDUCTIVITY=3.828_CMISSDP
  !surface area to volume ratio in [cm^-1]
  REAL(CMISSDP), PARAMETER :: Am=500.0_CMISSDP
  !membrane capacitance in [uF/cm^2]
  REAL(CMISSDP), PARAMETER :: Cm=1.0_CMISSDP

  !--------------------------------------------------------------------------------------------------------------------------------
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3

  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=NumberOfSpatialCoordinates
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussPoints=3

  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2

  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumberM=2
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=NumberOfSpatialCoordinates

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumberM=4
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariablesM=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponentsM=NumberOfSpatialCoordinates+2

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumberFE=5
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariablesFE=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponentsFE1=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponentsFE2=1

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumberM=6
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariablesM=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponentsM=1

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumberFE=7
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariablesFE=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponentsFE=NumberOfSpatialCoordinates+1

  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfComponents=1

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberM=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberFE=10

  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=14

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetsUserNumberM=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetsUserNumberFE=2

  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: SolverDAEIndex=1
  INTEGER(CMISSIntg), PARAMETER :: SolverParabolicIndex=2
  INTEGER(CMISSIntg), PARAMETER :: SolverFEIndex=1

  INTEGER(CMISSIntg), PARAMETER :: ControlLoopMonodomainNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopElasticityNumber=2
  
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: EquationsSetIndexM,EquationsSetIndexFE
  INTEGER(CMISSIntg) :: CellMLIndex
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx,ComponentNumber

  INTEGER(CMISSIntg),ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg) :: BottomNormalXi,LeftNormalXi,RightNormalXi,FrontNormalXi

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: shortenModelIndex
  INTEGER(CMISSIntg) :: stimcomponent

  REAL(CMISSDP) :: YVALUE

  INTEGER(CMISSIntg) :: Err

  !CMISS variables

  TYPE(CMISSBasisType) :: QuadraticBasis,LinearBasis,Basis(2)
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsM,BoundaryConditionsFE
  TYPE(CMISSCellMLType) :: CellML
  TYPE(CMISSCellMLEquationsType) :: CellMLEquations
  TYPE(CMISSControlLoopType) :: ControlLoopMain
  TYPE(CMISSControlLoopType) :: ControlLoopM,ControlLoopFE
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: EquationsM,EquationsFE
  TYPE(CMISSEquationsSetType) :: EquationsSetM,EquationsSetFE
  TYPE(CMISSFieldType) :: EquationsSetFieldM,EquationsSetFieldFE
  TYPE(CMISSFieldType) :: GeometricFieldM,GeometricFieldFE
  TYPE(CMISSFieldType) :: DependentFieldM,DependentFieldFE
  TYPE(CMISSFieldType) :: IndependentField
  TYPE(CMISSFieldType) :: MaterialFieldM,MaterialFieldFE
  TYPE(CMISSFieldType) :: FibreField
  TYPE(CMISSFieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: SolverDAE,SolverParabolic
  TYPE(CMISSSolverType) :: SolverFE,LinearSolverFE
  TYPE(CMISSSolverEquationsType) :: SolverEquationsM,SolverEquationsFE

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
   !Generic CMISS variables

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


  !================================================================================================================================
  !  G E N E R A L   F E A T U R E S
  !================================================================================================================================

  !--------------------------------------------------------------------------------------------------------------------------------
  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Trap errors
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
  
  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  CALL CMISSOutputSetOn("GudunovMonodomainElasticitySameMesh",Err)
  
  IF(NumberOfSpatialCoordinates==1) THEN
    WRITE(*,'(A)') "1D case not implemented."
    STOP
!    NumberGlobalXElements=NUMBER_OF_ELEMENTS
!    NumberGlobalYElements=0
!    NumberGlobalZElements=0
  ELSE IF(NumberOfSpatialCoordinates==2) THEN
    WRITE(*,'(A)') "2D case not implemented."
    STOP
!    NumberGlobalXElements=NUMBER_OF_ELEMENTS
!    NumberGlobalYElements=NUMBER_OF_ELEMENTS
!    NumberGlobalZElements=0
  ELSE IF(NumberOfSpatialCoordinates==3) THEN
    NumberGlobalXElements=NUMBER_OF_ELEMENTS
    NumberGlobalYElements=1
    NumberGlobalZElements=1
  ELSE
    WRITE(*,'(A)') "The number of spatial coordinates must be 1, 2 or 3."
    STOP
  ENDIF
  
  NumberOfDomains=NumberOfComputationalNodes
  
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      !Set the coordinate system to be 1D
      CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,1,Err)
    ELSE
      !Set the coordinate system to be 2D
      CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
    ENDIF
  ELSE
    !Set the coordinate system to be 3D
    CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_LabelSet(Region,"Region",Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the bases
  !Define basis functions - tri-Quadratic Lagrange 
  CALL CMISSBasis_Initialise(QuadraticBasis,Err)
  CALL CMISSBasis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasis_TypeSet(QuadraticBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(QuadraticBasis,NumberOfXiCoordinates,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      CALL CMISSBasis_InterpolationXiSet(QuadraticBasis,[CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticBasis,[NumberOfGaussPoints],Err)
    ELSE
      CALL CMISSBasis_InterpolationXiSet(QuadraticBasis,[CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
        & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticBasis,[NumberOfGaussPoints,NumberOfGaussPoints],Err)
    ENDIF
  ELSE
    CALL CMISSBasis_InterpolationXiSet(QuadraticBasis,[CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
      & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
      & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  ENDIF
  CALL CMISSBasis_CreateFinish(QuadraticBasis,Err)

  !Define basis functions - tri-Linear Lagrange
  CALL CMISSBasis_Initialise(LinearBasis,Err)
  CALL CMISSBasis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasis_TypeSet(LinearBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(LinearBasis,NumberOfXiCoordinates,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      CALL CMISSBasis_InterpolationXiSet(LinearBasis,[CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearBasis,[NumberOfGaussPoints],Err)
    ELSE
      CALL CMISSBasis_InterpolationXiSet(LinearBasis,[CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
        & CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearBasis,[NumberOfGaussPoints,NumberOfGaussPoints],Err)
    ENDIF
  ELSE
    CALL CMISSBasis_InterpolationXiSet(LinearBasis,[CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
      & CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearBasis, &
      & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  ENDIF
  CALL CMISSBasis_CreateFinish(LinearBasis,Err)

  Basis(1)=QuadraticBasis
  Basis(2)=LinearBasis

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the basis - Basis is an array: [QuadraticBasis,LinearBasis]
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH],Err)
      CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements],Err)
    ELSE
      CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH,WIDTH],Err)
      CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
    ENDIF
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH,WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
        & NumberGlobalZElements],Err)
  ENDIF    
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)


  !================================================================================================================================
  !  F I N I T E   E L A S T C I T Y
  !================================================================================================================================
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for finite elasticity - quadratic interpolation
  CALL CMISSField_Initialise(GeometricFieldFE,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumberFE,Region,GeometricFieldFE,Err)
  CALL CMISSField_TypeSet(GeometricFieldFE,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricFieldFE,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricFieldFE,CMISS_FIELD_GEOMETRIC_TYPE,Err)  
  CALL CMISSField_NumberOfVariablesSet(GeometricFieldFE,FieldGeometryNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  IF(NumberGlobalYElements/=0) THEN
    CALL CMISSField_ComponentMeshComponentSet(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL CMISSField_ComponentMeshComponentSet(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
    ENDIF
  ENDIF
  CALL CMISSField_VariableLabelSet(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL CMISSField_CreateFinish(GeometricFieldFE,Err)
  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricFieldFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a fibre field and attach it to the geometric field - quadratic interpolation
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricFieldFE,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err) 
  IF(NumberGlobalYElements/=0) THEN
    CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
    ENDIF
  ENDIF
  CALL CMISSField_VariableLabelSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL CMISSField_CreateFinish(FibreField,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a material field for Finite Elasticity and attach it to the geometric field - quadratic interpolation
  CALL CMISSField_Initialise(MaterialFieldFE,Err)
  CALL CMISSField_CreateStart(FieldMaterialUserNumberFE,Region,MaterialFieldFE,Err)
  CALL CMISSField_TypeSet(MaterialFieldFE,CMISS_FIELD_MATERIAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(MaterialFieldFE,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(MaterialFieldFE,GeometricFieldFE,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialFieldFE,FieldMaterialNumberOfVariablesFE,Err)
  CALL CMISSField_VariableTypesSet(MaterialFieldFE,[CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_V_VARIABLE_TYPE],Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponentsFE1,Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,FieldMaterialNumberOfComponentsFE2,Err)
! default is CMISS_FIELD_NODE_BASED_INTERPOLATION
  CALL CMISSField_ComponentMeshComponentSet(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(MaterialFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_VariableLabelSet(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,"MaterialFE",Err)
  CALL CMISSField_VariableLabelSet(MaterialFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,"Gravity",Err)
  CALL CMISSField_CreateFinish(MaterialFieldFE,Err)
  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,6.0_CMISSDP, &
    & Err)
!  CALL CMISSField_ComponentValuesInitialise(MaterialFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)

  !Create the dependent field for FE with 2 variables and * components 
  !  3-d: 3 displacement (quad interpol), 1 pressure (lin interpol) --> * = 4
  !  2-d: 2 displacement (quad interpol), 1 pressure (lin interpol) --> * = 3
  !  1-d: 1 displacement (quad interpol), 1 pressure (lin interpol) --> * = 2
  CALL CMISSField_Initialise(DependentFieldFE,Err)
  CALL CMISSField_CreateStart(FieldDependentUserNumberFE,Region,DependentFieldFE,Err)
  CALL CMISSField_TypeSet(DependentFieldFE,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentFieldFE,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentFieldFE,GeometricFieldFE,Err)
  CALL CMISSField_DependentTypeSet(DependentFieldFE,CMISS_FIELD_DEPENDENT_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentFieldFE,FieldDependentNumberOfVariablesFE,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponentsFE,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponentsFE,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,2,LinearMeshComponentNumber,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1, &
        & QuadraticMeshComponentNumber, &
        & Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2,LinearMeshComponentNumber,Err)
    ELSE
      CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,3,LinearMeshComponentNumber,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1, &
        & QuadraticMeshComponentNumber, &
        & Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2, &
        & QuadraticMeshComponentNumber, &
        & Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3,LinearMeshComponentNumber,Err)
    ENDIF
  ELSE
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber, &
      & Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber, &
      & Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber, &
      & Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  ENDIF
!  CALL CMISSField_ScalingTypeSet(DependentFieldFE,CMISS_FIELD_UNIT_SCALING,Err)
  CALL CMISSField_VariableLabelSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,"DependentFE",Err)
  CALL CMISSField_VariableLabelSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,"Reaction_Force",Err)
  CALL CMISSField_CreateFinish(DependentFieldFE,Err)



  !================================================================================================================================
  !  M O N O D O M A I N
  !================================================================================================================================
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for monodomain - quadratic interpolation
  CALL CMISSField_Initialise(GeometricFieldM,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumberM,Region,GeometricFieldM,Err)
  CALL CMISSField_TypeSet(GeometricFieldM,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricFieldM,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricFieldM,CMISS_FIELD_GEOMETRIC_TYPE,Err)  
  CALL CMISSField_NumberOfVariablesSet(GeometricFieldM,FieldGeometryNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricFieldM,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  IF(NumberGlobalYElements/=0) THEN
    CALL CMISSField_ComponentMeshComponentSet(GeometricFieldM,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL CMISSField_ComponentMeshComponentSet(GeometricFieldM,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
    ENDIF
  ENDIF
  CALL CMISSField_VariableLabelSet(GeometricFieldM,CMISS_FIELD_U_VARIABLE_TYPE,"GeometryM",Err)
  CALL CMISSField_CreateFinish(GeometricFieldM,Err)
  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricFieldM,Err)
        

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a materials field for monodomain and attach it to the geometric field - constant interpolation
  CALL CMISSField_Initialise(MaterialFieldM,Err)
  CALL CMISSField_CreateStart(FieldMaterialUserNumberM,Region,MaterialFieldM,Err)
  CALL CMISSField_TypeSet(MaterialFieldM,CMISS_FIELD_MATERIAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(MaterialFieldM,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(MaterialFieldM,GeometricFieldM,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialFieldM,FieldMaterialNumberOfVariablesM,Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponentsM,Err)
! default is CMISS_FIELD_NODE_BASED_INTERPOLATION
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,3,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  IF(NumberGlobalYElements/=0) THEN
    CALL CMISSField_ComponentInterpolationSet(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL CMISSField_ComponentInterpolationSet(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,5,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    ENDIF
  ENDIF
  CALL CMISSField_VariableLabelSet(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,"MaterialM",Err)
  CALL CMISSField_CreateFinish(MaterialFieldM,Err)
  !Set Am
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Am,Err)
  !Set Cm
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Cm,Err)
  !Set Conductivity
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,CONDUCTIVITY, &
    & Err)
  IF(NumberGlobalYElements/=0) THEN
    CALL CMISSField_ComponentValuesInitialise(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4, &
      & CONDUCTIVITY,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL CMISSField_ComponentValuesInitialise(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,5, &
        & CONDUCTIVITY,Err)
    ENDIF
  ENDIF

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the dependent field for monodomain with 2 variables and 1 components 
  CALL CMISSField_Initialise(DependentFieldM,Err)
  CALL CMISSField_CreateStart(FieldDependentUserNumberM,Region,DependentFieldM,Err)
  CALL CMISSField_TypeSet(DependentFieldM,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentFieldM,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentFieldM,GeometricFieldM,Err)
  CALL CMISSField_DependentTypeSet(DependentFieldM,CMISS_FIELD_DEPENDENT_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentFieldM,FieldDependentNumberOfVariablesM,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponentsM,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldM,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponentsM,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldM,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_DimensionSet(DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMISSField_DimensionSet(DependentFieldM,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMISSField_VariableLabelSet(DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,"Vm",Err)
  CALL CMISSField_VariableLabelSet(DependentFieldM,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,"dVm/dt",Err)
  CALL CMISSField_CreateFinish(DependentFieldM,Err)


  !================================================================================================================================
  !  INDEPENDENT FIELD

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the independent field for the active stress determined in the cell model of the monodomain equation and used in finite elasticity
!  independent_field_auto_create = .TRUE.
  CALL CMISSField_Initialise(IndependentField,Err)
  IF(.NOT. independent_field_auto_create) THEN
    CALL CMISSField_CreateStart(FieldIndependentUserNumber,Region,IndependentField,Err)
    CALL CMISSField_TypeSet(IndependentField,CMISS_FIELD_GENERAL_TYPE,Err)
    CALL CMISSField_MeshDecompositionSet(IndependentField,Decomposition,Err)
    CALL CMISSField_GeometricFieldSet(IndependentField,GeometricFieldM,Err)
    CALL CMISSField_DependentTypeSet(IndependentField,CMISS_FIELD_INDEPENDENT_TYPE,Err)
    CALL CMISSField_NumberOfVariablesSet(IndependentField,FieldIndependentNumberOfVariables,Err)
    CALL CMISSField_DimensionSet(IndependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
    CALL CMISSField_NumberOfComponentsSet(IndependentField,CMISS_FIELD_U_VARIABLE_TYPE,FieldIndependentNumberOfComponents,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_NODE_BASED_INTERPOLATION, &
      & Err)
    CALL CMISSField_ComponentMeshComponentSet(IndependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_VariableLabelSet(IndependentField,CMISS_FIELD_U_VARIABLE_TYPE,"Acive_Stress",Err)
    CALL CMISSField_CreateFinish(IndependentField,Err)
  ENDIF


  !================================================================================================================================
  !  EQUATIONS SET

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for Finite Elasticity
  CALL CMISSField_Initialise(EquationsSetFieldFE,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetFE,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetsUserNumberFE,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
      & EquationsSetFieldUserNumberFE, &
    & EquationsSetFieldFE,EquationsSetFE,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSetFE,Err)

  !Create the equations set dependent field variables for Finite Elasticity
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetFE,FieldDependentUserNumberFE,DependentFieldFE,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetFE,Err)

  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL CMISSEquationsSet_IndependentCreateStart(EquationsSetFE,FieldIndependentUserNumber,IndependentField,Err)
  CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSetFE,Err)

  !Create the equations set materials field variables for Finite Elasticity
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetFE,FieldMaterialUserNumberFE,MaterialFieldFE,Err)  
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for monodomain
  CALL CMISSField_Initialise(EquationsSetFieldM,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetM,Err)
  !Set the equations set to be a Monodomain equations set
  !> \todo solve the monodomain problem on the fibre field rather than on the geometric field: GeometricField <--> FibreField
  CALL CMISSEquationsSet_CreateStart(EquationsSetsUserNumberM,Region,GeometricFieldM,CMISS_EQUATIONS_SET_BIOELECTRICS_CLASS, &
    & CMISS_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMISS_EQUATIONS_SET_NO_SUBTYPE,EquationsSetFieldUserNumberM,EquationsSetFieldM, &
    & EquationsSetM,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSetM,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set dependent field variables for monodomain
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetM,FieldDependentUserNumberM,DependentFieldM,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetM,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set materials field variables for monodomain
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetM,FieldMaterialUserNumberM,MaterialFieldM,Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetM,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set equations for monodomain
  CALL CMISSEquations_Initialise(EquationsM,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetM,EquationsM,Err)
  CALL CMISSEquations_SparsityTypeSet(EquationsM,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(EquationsM,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetM,Err)

  !Create the equations set equations for Finite Elasticity
  CALL CMISSEquations_Initialise(EquationsFE,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetFE,EquationsFE,Err)
  CALL CMISSEquations_SparsityTypeSet(EquationsFE,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(EquationsFE,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetFE,Err)   

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML environment
  CALL CMISSCellML_Initialise(CellML,Err)
  CALL CMISSCellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import the Shorten et al. 2007 model from a file
  CALL CMISSCellML_ModelImport(CellML,"shorten_mod_2011_07_04.xml",shortenModelIndex,Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL CMISSCellML_VariableSetAsKnown(CellML,shortenModelIndex,"wal_environment/I_HH",Err)
  !   - to get from the CellML side
!  CALL CMISSCellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_T",Err)
!  CALL CMISSCellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_ionic_s",Err)
!  CALL CMISSCellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_ionic_t",Err)
  !
  !NOTE: If an INTERMEDIATE (or ALGEBRAIC) variable should be used in a mapping, it has to be set as known or wanted first!
  !       --> set "razumova/stress" as wanted!
  !       --> no need to set "wal_environment/vS" since all STATE variables are automatically set as wanted! 
  CALL CMISSCellML_VariableSetAsWanted(CellML,shortenModelIndex,"razumova/stress",Err)
  !   - and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
  !Finish the CellML environment
  CALL CMISSCellML_CreateFinish(CellML,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML <--> OpenCMISS field maps
  CALL CMISSCellML_FieldMapsCreateStart(CellML,Err)
  !Map the transmembrane voltage V_m
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & shortenModelIndex,"wal_environment/vS",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"wal_environment/vS",CMISS_FIELD_VALUES_SET_TYPE, &
    & DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"razumova/stress",CMISS_FIELD_VALUES_SET_TYPE, &
    & IndependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_FieldMapsCreateFinish(CellML,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Initialise dependent field for monodomain
  !> \todo - get V_m initialial value.
  CALL CMISSField_ComponentValuesInitialise(DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & -79.974_CMISSDP,Err)
  
  !Initialise dependent field for Finite Elasticity from undeformed geometry and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  CALL CMISSField_ComponentValuesInitialise(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4, &
    & -8.0_CMISSDP,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML models field
  CALL CMISSField_Initialise(CellMLModelsField,Err)
  CALL CMISSCellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  CALL CMISSCellML_ModelsFieldCreateFinish(CellML,Err)

  !Create the CellML state field
  CALL CMISSField_Initialise(CellMLStateField,Err)
  CALL CMISSCellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,Err)
  CALL CMISSCellML_StateFieldCreateFinish(CellML,Err)

  !Create the CellML intermediate field
  CALL CMISSField_Initialise(CellMLIntermediateField,Err)
  CALL CMISSCellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  CALL CMISSCellML_IntermediateFieldCreateFinish(CellML,Err)
  
  !Create the CellML parameters field
  CALL CMISSField_Initialise(CellMLParametersField,Err)
  CALL CMISSCellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  CALL CMISSCellML_ParametersFieldCreateFinish(CellML,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Set the Stimulus for monodomain at the left side
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi, &
    & Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)

  CALL CMISSCellML_FieldComponentGet(CellML,shortenModelIndex,CMISS_CELLML_PARAMETERS_FIELD, &
    & "wal_environment/I_HH",stimcomponent,Err)

  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      SELECT CASE(NumberOfSpatialCoordinates)
      CASE(1)
        WRITE(*,'(A)') "Setting the stimulus not implemented for 1D."
        STOP
      CASE(2)
        !simulate the lower half of the nodes on the left surface
        CALL CMISSField_ParameterSetGetNode(GeometricFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
          & NodeNumber, &
          & 2, &
          & YVALUE,Err)
        IF(YVALUE .LE. WIDTH/2.0_CMISSDP) CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField, &
          & CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,STIM_VALUE,Err)
      CASE(3)
        !stimulate all nodes on the left surface
        CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,STIM_VALUE,Err)
      END SELECT
    ENDIF
  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
  !Define the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_MULTI_PHYSICS_CLASS,CMISS_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE, &
    & CMISS_PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)

  !set the main control loop (time loop type)
  CALL CMISSControlLoop_Initialise(ControlLoopMain,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoopMain,Err)
  CALL CMISSControlLoop_LabelSet(ControlLoopMain,'MAIN_TIME_LOOP',Err)
  !Loop in time for STIM_STOP with the Stimulus applied.
  CALL CMISSControlLoop_TimesSet(ControlLoopMain,0.0_CMISSDP,STIM_STOP,ELASTICITY_TIME_STEP,Err)
  CALL CMISSControlLoop_TimeOutputSet(ControlLoopMain,OUTPUT_FREQUENCY,Err)
  CALL CMISSControlLoop_OutputTypeSet(ControlLoopMain,CMISS_CONTROL_LOOP_TIMING_OUTPUT,Err)

  !set the monodomain loop (time loop type)
  CALL CMISSControlLoop_Initialise(ControlLoopM,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,[ControlLoopMonodomainNumber,CMISS_CONTROL_LOOP_NODE],ControlLoopM,Err)
  CALL CMISSControlLoop_LabelSet(ControlLoopM,'MONODOMAIN_TIME_LOOP',Err)
  CALL CMISSControlLoop_TimesSet(ControlLoopM,0.0_CMISSDP,ELASTICITY_TIME_STEP,PDE_TIME_STEP,Err)
  CALL CMISSControlLoop_OutputTypeSet(ControlLoopM,CMISS_CONTROL_LOOP_NO_OUTPUT,Err)

  !set the finite elasticity loop (simple type)
  CALL CMISSControlLoop_Initialise(ControlLoopFE,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,[ControlLoopElasticityNumber,CMISS_CONTROL_LOOP_NODE],ControlLoopFE,Err)
  CALL CMISSControlLoop_LabelSet(ControlLoopFE,'ELASTICITY_LOOP',Err)

  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solvers
  CALL CMISSProblem_SolversCreateStart(Problem,Err)

  !Create the DAE solver
  CALL CMISSSolver_Initialise(SolverDAE,Err)
  CALL CMISSProblem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMISS_CONTROL_LOOP_NODE], &
    & SolverDAEIndex,SolverDAE,Err)
  CALL CMISSSolver_DAETimeStepSet(SolverDAE,ODE_TIME_STEP,Err)
  !> \todo - solve the CellML equations on the GPU for efficiency (later)
  !CALL CMISSSolver_DAESolverTypeSet(SolverDAE,CMISS_SOLVER_DAE_EXTERNAL,Err) 
  CALL CMISSSolver_OutputTypeSet(SolverDAE,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverDAE,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverDAE,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverDAE,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverDAE,CMISS_SOLVER_MATRIX_OUTPUT,Err)

  !Create the parabolic solver
  CALL CMISSSolver_Initialise(SolverParabolic,Err)
  CALL CMISSProblem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMISS_CONTROL_LOOP_NODE], &
    & SolverParabolicIndex,SolverParabolic,Err)
  CALL CMISSSolver_DynamicSchemeSet(SolverParabolic,CMISS_SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME,Err)
  CALL CMISSSolver_OutputTypeSet(SolverParabolic,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverParabolic,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverParabolic,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverParabolic,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverParabolic,CMISS_SOLVER_MATRIX_OUTPUT,Err)

  !Create the Finte Elasticity solver
  CALL CMISSSolver_Initialise(SolverFE,Err)
  CALL CMISSSolver_Initialise(LinearSolverFE,Err)
  CALL CMISSProblem_SolverGet(Problem,[ControlLoopElasticityNumber,CMISS_CONTROL_LOOP_NODE], &
    & SolverFEIndex,SolverFE,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverFE,CMISS_SOLVER_NO_OUTPUT,Err)
  CALL CMISSSolver_OutputTypeSet(SolverFE,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverFE,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverFE,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverFE,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(SolverFE,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(SolverFE,LinearSolverFE,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolverFE,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)

  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver CellML equations
  CALL CMISSProblem_CellMLEquationsCreateStart(Problem,Err)

  CALL CMISSSolver_Initialise(SolverDAE,Err)
  CALL CMISSCellMLEquations_Initialise(CellMLEquations,Err)
  CALL CMISSProblem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMISS_CONTROL_LOOP_NODE], &
    & SolverDAEIndex,SolverDAE,Err)
  CALL CMISSSolver_CellMLEquationsGet(SolverDAE,CellMLEquations,Err)
  CALL CMISSCellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  CALL CMISSProblem_CellMLEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)

  !Create the problem solver parabolic equations (Monodomain)
  CALL CMISSSolver_Initialise(SolverParabolic,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsM,Err)
  CALL CMISSProblem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMISS_CONTROL_LOOP_NODE], &
    & SolverParabolicIndex,SolverParabolic,Err)
  CALL CMISSSolver_SolverEquationsGet(SolverParabolic,SolverEquationsM,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsM,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsM,CMISS_SOLVER_FULL_MATRICES,Err)  
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsM,EquationsSetM,EquationsSetIndexM,Err)

  !Create the problem solver Finite Elasticity equations
  CALL CMISSSolver_Initialise(SolverFE,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsFE,Err)
  CALL CMISSProblem_SolverGet(Problem,[ControlLoopElasticityNumber,CMISS_CONTROL_LOOP_NODE], &
    & SolverFEIndex,SolverFE,Err)
  CALL CMISSSolver_SolverEquationsGet(SolverFE,SolverEquationsFE,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsFE,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsFE,CMISS_SOLVER_FULL_MATRICES,Err)  
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsFE,EquationsSetFE,EquationsSetIndexFE,Err)

  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !boundary conditions

  !Prescribe boundary conditions for monodomain
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsM,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsM,BoundaryConditionsM,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsM,Err)

  !Prescribe boundary conditions for Finite Elasticity (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsFE,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsFE,BoundaryConditionsFE,Err)

  SELECT CASE(NumberOfSpatialCoordinates)
  CASE(3)

    !  uniaxial_extension_bc=.TRUE.
    IF(uniaxial_extension_bc) THEN
      !Set x=0 nodes to no x displacment in x
      DO node_idx=1,SIZE(LeftSurfaceNodes,1)
        NodeNumber=LeftSurfaceNodes(node_idx)
        CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
            & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        ENDIF
      ENDDO
      !Set x=WIDTH nodes to 10% x displacement
      DO node_idx=1,SIZE(RightSurfaceNodes,1)
        NodeNumber=RightSurfaceNodes(node_idx)
        CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
            & CMISS_BOUNDARY_CONDITION_FIXED,1.1_CMISSDP*LENGTH,Err)
        ENDIF
      ENDDO
      !Set y=0 nodes to no y displacement
      DO node_idx=1,SIZE(FrontSurfaceNodes,1)
        NodeNumber=FrontSurfaceNodes(node_idx)
        CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
            & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        ENDIF
      ENDDO
      !Set z=0 nodes to no z displacement
      DO node_idx=1,SIZE(BottomSurfaceNodes,1)
        NodeNumber=BottomSurfaceNodes(node_idx)
        CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
            & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        ENDIF
      ENDDO

    ELSE !uniaxial_extension_bc

      DO node_idx=1,SIZE(LeftSurfaceNodes,1)
        NodeNumber=LeftSurfaceNodes(node_idx)
        CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          DO ComponentNumber=2,3
            CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
              & NodeNumber, &
              & ComponentNumber,YVALUE,Err)
            IF(YVALUE<tol) THEN
              !fix all nodes (0,0,*) in x- and y-direction
              !fix all nodes (0,*,0) in x- and z-direction
              CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
                & NodeNumber, &
                & 1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
              CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
                & NodeNumber, &
                & ComponentNumber,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      !fix all (*,0,0) in y- and z-direction
      DO node_idx=1,SIZE(BottomSurfaceNodes,1)
        NodeNumber=BottomSurfaceNodes(node_idx)
        CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
            & NodeNumber, &
            & 2,YVALUE,Err)
          IF(YVALUE<tol) THEN
            CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
              & 2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
            CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
              & 3,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
          ENDIF
        ENDIF
      ENDDO
      
    ENDIF !uniaxial_extension_bc
    
  CASE DEFAULT
    WRITE(*,'(A)') "Boundary conditions not implemented for 1D and 2D."
    STOP
  END SELECT
  
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem for the first STIM_STOP
  CALL CMISSProblem_Solve(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Now turn the stimulus off
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      SELECT CASE(NumberOfSpatialCoordinates)
      CASE(1)
        WRITE(*,'(A)') "Setting the stimulus not implemented for 1D."
        STOP
      CASE(2)
        CALL CMISSField_ParameterSetGetNode(GeometricFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
          & NodeNumber, &
          & 2, &
          & YVALUE,Err)
        IF(YVALUE .LE. WIDTH/2.0_CMISSDP) CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField, &
          & CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,0.0_CMISSDP,Err)
      CASE(3)
        CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,0.0_CMISSDP,Err)
      END SELECT
    ENDIF
  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
  !Set the time loop from STIM_STOP to TIME_STOP
  CALL CMISSControlLoop_TimesSet(ControlLoopMain,STIM_STOP,TIME_STOP,ELASTICITY_TIME_STEP,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem for the rest of the time
  CALL CMISSProblem_Solve(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"GudunovMonodomainElasticitySameMeshExample","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"GudunovMonodomainElasticitySameMeshExample","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM GODUNOVMONODOMAINELASTICITYSAMEMESHEXAMPLE
