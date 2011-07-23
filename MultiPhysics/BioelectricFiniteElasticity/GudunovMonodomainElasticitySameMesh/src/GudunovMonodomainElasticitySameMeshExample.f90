!> \file
!> $Id: GudunovMonodomainElasticitySameMeshExample.f90 -1   $
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



  !Test program parameters

  REAL(CMISSDP), PARAMETER :: LENGTH= 2.5_CMISSDP !1.0_CMISSDP ! X-direction
  REAL(CMISSDP), PARAMETER :: WIDTH=  2.5_CMISSDP !1.0_CMISSDP ! Y-direction
  REAL(CMISSDP), PARAMETER :: HEIGHT=0.25_CMISSDP !1.0_CMISSDP ! Z-direction

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3

  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=NumberOfSpatialCoordinates
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussPoints=3

  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
!  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=3
!  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2

  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=NumberOfSpatialCoordinates

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumberM=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariablesM=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponentsM=NumberOfSpatialCoordinates+2

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumberFE=4
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariablesFE=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponentsFE1=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponentsFE2=1

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumberM=5
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariablesM=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponentsM=1

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumberFE=6
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariablesFE=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponentsFE=NumberOfSpatialCoordinates+1

  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfComponents=1

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberM=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberFE=9

  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=14

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberM=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberFE=2

  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumberM=1
  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumberFE=2

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
  
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain

  INTEGER(CMISSIntg),ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg) :: BottomNormalXi,LeftNormalXi,RightNormalXi,FrontNormalXi

  LOGICAL :: EXPORT_FIELD

!  INTEGER(CMISSIntg) :: N,CELL_TYPE
!  INTEGER(CMISSIntg) :: n98ModelIndex!,JRWModelIndex,LRdModelIndex
  INTEGER(CMISSIntg) :: shortenModelIndex
!  INTEGER(CMISSIntg) :: gK1component,gNacomponent,
  INTEGER(CMISSIntg) :: stimcomponent

  REAL(CMISSDP) :: YVALUE

!  REAL(CMISSDP) :: X,Y,DISTANCE,gK1_VALUE,gNa_VALUE
  
  REAL(CMISSDP), PARAMETER :: STIM_VALUE = 30000.0_CMISSDP !100.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: STIM_STOP = 0.02_CMISSDP !0.1_CMISSDP
  REAL(CMISSDP), PARAMETER :: TIME_STOP = 2.0_CMISSDP !0.2_CMISSDP
  REAL(CMISSDP), PARAMETER :: ODE_TIME_STEP = 0.00001_CMISSDP !0.00001_CMISSDP
  REAL(CMISSDP), PARAMETER :: PDE_TIME_STEP = 0.0005_CMISSDP !0.001_CMISSDP
  REAL(CMISSDP), PARAMETER :: ELASTICITY_TIME_STEP = 0.001_CMISSDP !0.001_CMISSDP
!  REAL(CMISSDP), PARAMETER :: ELASTICITY_TIME_STEP = 0.01_CMISSDP !0.001_CMISSDP
  REAL(CMISSDP), PARAMETER :: CONDUCTIVITY = 3.828_CMISSDP !3.828_CMISSDP!0.1_CMISSDP
  REAL(CMISSDP), PARAMETER :: Am = 500_CMISSDP !193.6_CMISSDP
  REAL(CMISSDP), PARAMETER :: Cm = 1.0_CMISSDP !0.014651_CMISSDP

!  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_ELEMENTS=2

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
  TYPE(CMISSFieldType) :: GeometricField
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

  !--------------------------------------------------------------------------------------------------------------------------------
  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Trap errors
  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)
  
  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  CALL CMISSOutputSetOn("GudunovMonodomainElasticitySameMesh",Err)
  
  IF(NumberOfSpatialCoordinates==1) THEN
    NumberGlobalXElements=10
    NumberGlobalYElements=0
    NumberGlobalZElements=0
  ELSE IF(NumberOfSpatialCoordinates==2) THEN
    NumberGlobalXElements=20
    NumberGlobalYElements=20
    NumberGlobalZElements=0
  ELSE IF(NumberOfSpatialCoordinates==3) THEN
    NumberGlobalXElements=10
    NumberGlobalYElements=10
    NumberGlobalZElements=1
!    NumberGlobalYElements=NUMBER_OF_ELEMENTS
!    NumberGlobalZElements=NUMBER_OF_ELEMENTS
  ELSE
    WRITE(*,'(A)') "The number of spatial coordinates must be (1,) 2 or 3."
  ENDIF
  
  NumberOfDomains=NumberOfComputationalNodes
  
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      !Set the coordinate system to be 1D
      CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,1,Err)
    ELSE
      !Set the coordinate system to be 2D
      CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,2,Err)
    ENDIF
  ELSE
    !Set the coordinate system to be 3D
    CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionLabelSet(Region,"Region",Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of a basis (default is trilinear lagrange)
  !Define basis functions - tri-Quadratic Lagrange 
  CALL CMISSBasisTypeInitialise(QuadraticBasis,Err)
  CALL CMISSBasisCreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasisTypeSet(QuadraticBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(QuadraticBasis,NumberOfXiCoordinates,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      CALL CMISSBasisInterpolationXiSet(QuadraticBasis,[CMISSBasisQuadraticLagrangeInterpolation],Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticBasis,[NumberOfGaussPoints],Err)
    ELSE
      CALL CMISSBasisInterpolationXiSet(QuadraticBasis,[CMISSBasisQuadraticLagrangeInterpolation, &
        & CMISSBasisQuadraticLagrangeInterpolation],Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticBasis,[NumberOfGaussPoints,NumberOfGaussPoints],Err)
    ENDIF
  ELSE
    CALL CMISSBasisInterpolationXiSet(QuadraticBasis,[CMISSBasisQuadraticLagrangeInterpolation, &
      & CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation],Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticBasis, &
      & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  ENDIF
  CALL CMISSBasisCreateFinish(QuadraticBasis,Err)

  !Define basis functions - tri-Linear Lagrange
  CALL CMISSBasisTypeInitialise(LinearBasis,Err)
  CALL CMISSBasisCreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasisTypeSet(LinearBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(LinearBasis,NumberOfXiCoordinates,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      CALL CMISSBasisInterpolationXiSet(LinearBasis,[CMISSBasisLinearLagrangeInterpolation],Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(LinearBasis,[NumberOfGaussPoints],Err)
    ELSE
      CALL CMISSBasisInterpolationXiSet(LinearBasis,[CMISSBasisLinearLagrangeInterpolation, &
        & CMISSBasisLinearLagrangeInterpolation],Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(LinearBasis,[NumberOfGaussPoints,NumberOfGaussPoints],Err)
    ENDIF
  ELSE
    CALL CMISSBasisInterpolationXiSet(LinearBasis,[CMISSBasisLinearLagrangeInterpolation, &
      & CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation],Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(LinearBasis, &
      & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  ENDIF
  CALL CMISSBasisCreateFinish(LinearBasis,Err)

  Basis(1)=QuadraticBasis
  Basis(2)=LinearBasis

  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the basis - Basis is now an array: [QuadraticBasis,LinearBasis]
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[LENGTH],Err)
      CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements],Err)
    ELSE
      CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[LENGTH,WIDTH],Err)
      CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
    ENDIF
  ELSE
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[LENGTH,WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
        & NumberGlobalZElements],Err)
  ENDIF    
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecompositionCalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a field to put the geometry (default is geometry) - quadratic interpolation
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)  
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
  IF(NumberGlobalYElements/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,QuadraticMeshComponentNumber,Err)
    ENDIF
  ENDIF
  CALL CMISSFieldVariableLabelSet(GeometricField,CMISSFieldUVariableType,"Geometry",Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)
        
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a fibre field and attach it to the geometric field - linear interpolation
  CALL CMISSFieldTypeInitialise(FibreField,Err)
  CALL CMISSFieldCreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSFieldTypeSet(FibreField,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreField,CMISSFieldUVariableType,FieldFibreNumberOfComponents,Err)  
!LinearMeshComponentNumber QuadraticMeshComponentNumber
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,1,LinearMeshComponentNumber,Err) 
  IF(NumberGlobalYElements/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,2,LinearMeshComponentNumber,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,3,LinearMeshComponentNumber,Err)
    ENDIF
  ENDIF
  CALL CMISSFieldVariableLabelSet(FibreField,CMISSFieldUVariableType,"Fibre",Err)
  CALL CMISSFieldCreateFinish(FibreField,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a materials field for monodomain and attach it to the geometric field - constant interpolation
  CALL CMISSFieldTypeInitialise(MaterialFieldM,Err)
  CALL CMISSFieldCreateStart(FieldMaterialUserNumberM,Region,MaterialFieldM,Err)
  CALL CMISSFieldTypeSet(MaterialFieldM,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialFieldM,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(MaterialFieldM,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialFieldM,FieldMaterialNumberOfVariablesM,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialFieldM,CMISSFieldUVariableType,FieldMaterialNumberOfComponentsM,Err)
! default is CMISSFieldNodeBasedInterpolation
  CALL CMISSFieldComponentInterpolationSet(MaterialFieldM,CMISSFieldUVariableType,1,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(MaterialFieldM,CMISSFieldUVariableType,2,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(MaterialFieldM,CMISSFieldUVariableType,3,CMISSFieldConstantInterpolation,Err)
  IF(NumberGlobalYElements/=0) THEN
    CALL CMISSFieldComponentInterpolationSet(MaterialFieldM,CMISSFieldUVariableType,4,CMISSFieldConstantInterpolation,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL CMISSFieldComponentInterpolationSet(MaterialFieldM,CMISSFieldUVariableType,5,CMISSFieldConstantInterpolation,Err)
    ENDIF
  ENDIF
!  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldM,CMISSFieldUVariableType,1,LinearMeshComponentNumber,Err) 
!  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldM,CMISSFieldUVariableType,2,LinearMeshComponentNumber,Err) 
!  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldM,CMISSFieldUVariableType,3,LinearMeshComponentNumber,Err) 
!  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldM,CMISSFieldUVariableType,4,LinearMeshComponentNumber,Err) 
!  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldM,CMISSFieldUVariableType,5,LinearMeshComponentNumber,Err) 
  CALL CMISSFieldVariableLabelSet(MaterialFieldM,CMISSFieldUVariableType,"MaterialM",Err)
  CALL CMISSFieldCreateFinish(MaterialFieldM,Err)
  !Set Am
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldM,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Am,Err)
  !Set Cm
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldM,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Cm,Err)
  !Set Conductivity
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldM,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,CONDUCTIVITY,Err)
  IF(NumberGlobalYElements/=0) THEN
    CALL CMISSFieldComponentValuesInitialise(MaterialFieldM,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,CONDUCTIVITY,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL CMISSFieldComponentValuesInitialise(MaterialFieldM,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,CONDUCTIVITY,Err)
    ENDIF
  ENDIF

  !Create a material field for Finite Elasticity and attach it to the geometric field - constant interpolation
  CALL CMISSFieldTypeInitialise(MaterialFieldFE,Err)
  CALL CMISSFieldCreateStart(FieldMaterialUserNumberFE,Region,MaterialFieldFE,Err)
  CALL CMISSFieldTypeSet(MaterialFieldFE,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialFieldFE,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(MaterialFieldFE,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialFieldFE,FieldMaterialNumberOfVariablesFE,Err)
  CALL CMISSFieldVariableTypesSet(MaterialFieldFE,[CMISSFieldUVariableType,CMISSFieldVVariableType],Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialFieldFE,CMISSFieldUVariableType,FieldMaterialNumberOfComponentsFE1,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialFieldFE,CMISSFieldVVariableType,FieldMaterialNumberOfComponentsFE2,Err)
! default is CMISSFieldNodeBasedInterpolation
  CALL CMISSFieldComponentInterpolationSet(MaterialFieldFE,CMISSFieldUVariableType,1,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(MaterialFieldFE,CMISSFieldUVariableType,2,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(MaterialFieldFE,CMISSFieldVVariableType,1,CMISSFieldConstantInterpolation,Err)
!  CALL CMISSFieldComponentInterpolationSet(MaterialFieldFE,CMISSFieldVVariableType,2,CMISSFieldConstantInterpolation,Err)
!  CALL CMISSFieldComponentInterpolationSet(MaterialFieldFE,CMISSFieldVVariableType,3,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldVariableLabelSet(MaterialFieldFE,CMISSFieldUVariableType,"MaterialFE",Err)
  CALL CMISSFieldCreateFinish(MaterialFieldFE,Err)
  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldFE,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldFE,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,6.0_CMISSDP,Err)
!  CALL CMISSFieldComponentValuesInitialise(MaterialFieldFE,CMISSFieldVVariableType,CMISSFieldValuesSetType,1,0.0_CMISSDP,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the dependent field for monodomain with 2 variables and 1 components 
  CALL CMISSFieldTypeInitialise(DependentFieldM,Err)
  CALL CMISSFieldCreateStart(FieldDependentUserNumberM,Region,DependentFieldM,Err)
  CALL CMISSFieldTypeSet(DependentFieldM,CMISSFieldGeneralType,Err)
  CALL CMISSFieldMeshDecompositionSet(DependentFieldM,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(DependentFieldM,GeometricField,Err)
  CALL CMISSFieldDependentTypeSet(DependentFieldM,CMISSFieldDependentType,Err)
  CALL CMISSFieldNumberOfVariablesSet(DependentFieldM,FieldDependentNumberOfVariablesM,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldM,CMISSFieldUVariableType,FieldDependentNumberOfComponentsM,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldM,CMISSFieldDelUDelNVariableType,FieldDependentNumberOfComponentsM,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldM,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldM,CMISSFieldDelUDelNVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldDimensionSet(DependentFieldM,CMISSFieldUVariableType,CMISSFieldScalarDimensionType,Err)
  CALL CMISSFieldDimensionSet(DependentFieldM,CMISSFieldDelUDelNVariableType,CMISSFieldScalarDimensionType,Err)
!  CALL CMISSFieldScalingTypeSet(DependentFieldM,CMISSFieldUnitScaling,Err)
  CALL CMISSFieldVariableLabelSet(DependentFieldM,CMISSFieldUVariableType,"DependentM",Err)
  CALL CMISSFieldCreateFinish(DependentFieldM,Err)

  !Create the dependent field for FE with 2 variables and * components 
  !  3-d: 3 displacement (quad interpol), 1 pressure (lin interpol) --> * = 4
  !  2-d: 2 displacement (quad interpol), 1 pressure (lin interpol) --> * = 3
  !  1-d: 1 displacement (quad interpol), 1 pressure (lin interpol) --> * = 2
  CALL CMISSFieldTypeInitialise(DependentFieldFE,Err)
  CALL CMISSFieldCreateStart(FieldDependentUserNumberFE,Region,DependentFieldFE,Err)
  CALL CMISSFieldTypeSet(DependentFieldFE,CMISSFieldGeneralType,Err)
  CALL CMISSFieldMeshDecompositionSet(DependentFieldFE,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(DependentFieldFE,GeometricField,Err)
  CALL CMISSFieldDependentTypeSet(DependentFieldFE,CMISSFieldDependentType,Err)
  CALL CMISSFieldNumberOfVariablesSet(DependentFieldFE,FieldDependentNumberOfVariablesFE,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldFE,CMISSFieldUVariableType,FieldDependentNumberOfComponentsFE,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldFE,CMISSFieldDelUDelNVariableType,FieldDependentNumberOfComponentsFE,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldUVariableType,2,LinearMeshComponentNumber,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldDelUDelNVariableType,1,QuadraticMeshComponentNumber,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldDelUDelNVariableType,2,LinearMeshComponentNumber,Err)
    ELSE
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldUVariableType,3,LinearMeshComponentNumber,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldDelUDelNVariableType,1,QuadraticMeshComponentNumber,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldDelUDelNVariableType,2,QuadraticMeshComponentNumber,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldDelUDelNVariableType,3,LinearMeshComponentNumber,Err)
    ENDIF
  ELSE
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldUVariableType,3,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldUVariableType,4,LinearMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldDelUDelNVariableType,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldDelUDelNVariableType,2,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldDelUDelNVariableType,3,QuadraticMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldFE,CMISSFieldDelUDelNVariableType,4,LinearMeshComponentNumber,Err)
  ENDIF
!  CALL CMISSFieldScalingTypeSet(DependentFieldFE,CMISSFieldUnitScaling,Err)
  CALL CMISSFieldVariableLabelSet(DependentFieldFE,CMISSFieldUVariableType,"DependentFE",Err)
  CALL CMISSFieldCreateFinish(DependentFieldFE,Err)

  CALL CMISSFieldTypeInitialise(IndependentField,Err)
  CALL CMISSFieldCreateStart(FieldIndependentUserNumber,Region,IndependentField,Err)
  CALL CMISSFieldTypeSet(IndependentField,CMISSFieldGeneralType,Err)
  CALL CMISSFieldMeshDecompositionSet(IndependentField,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(IndependentField,GeometricField,Err)
  CALL CMISSFieldDependentTypeSet(IndependentField,CMISSFieldIndependentType,Err)
  CALL CMISSFieldNumberOfVariablesSet(IndependentField,FieldIndependentNumberOfVariables,Err)
  CALL CMISSFieldDimensionSet(IndependentField,CMISSFieldUVariableType,CMISSFieldScalarDimensionType,Err)
  CALL CMISSFieldNumberOfComponentsSet(IndependentField,CMISSFieldUVariableType,FieldIndependentNumberOfComponents,Err)
  CALL CMISSFieldComponentInterpolationSet(IndependentField,CMISSFieldUVariableType,1,CMISSFieldNodeBasedInterpolation,Err)
  CALL CMISSFieldComponentMeshComponentSet(IndependentField,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
!  CALL CMISSFieldScalingTypeSet(IndependentField,CMISSFieldUnitScaling,Err)
  CALL CMISSFieldVariableLabelSet(IndependentField,CMISSFieldUVariableType,"Independent (acive stress)",Err)
  CALL CMISSFieldCreateFinish(IndependentField,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for monodomain
  CALL CMISSFieldTypeInitialise(EquationsSetFieldM,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSetM,Err)
  !Set the equations set to be a Monodomain equations set
  !GeometricField   <--->   FibreField ! TODO later
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberM,Region,GeometricField,CMISSEquationsSetBioelectricsClass, &
    & CMISSEquationsSetMonodomainEquationType,CMISSEquationsSetNoSubtype,EquationsSetFieldUserNumberM,EquationsSetFieldM, &
    & EquationsSetM,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSetM,Err)

  !Create the equations_set for Finite Elasticity
  CALL CMISSFieldTypeInitialise(EquationsSetFieldFE,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSetFE,Err)
  CALL CMISSEquationsSetCreateStart(EquationSetUserNumberFE,Region,FibreField,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetStandardMonodomainElasticitySubtype,EquationsSetFieldUserNumberFE, &
    & EquationsSetFieldFE,EquationsSetFE,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSetFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set dependent field variables for monodomain
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetM,FieldDependentUserNumberM,DependentFieldM,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetM,Err)

  !Create the equations set dependent field variables for Finite Elasticity
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetFE,FieldDependentUserNumberFE,DependentFieldFE,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL CMISSEquationsSetIndependentCreateStart(EquationsSetFE,FieldIndependentUserNumber,IndependentField,Err)
  CALL CMISSEquationsSetIndependentCreateFinish(EquationsSetFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set materials field variables for monodomain
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetM,FieldMaterialUserNumberM,MaterialFieldM,Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetM,Err)

  !Create the equations set materials field variables for Finite Elasticity
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetFE,FieldMaterialUserNumberFE,MaterialFieldFE,Err)  
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set equations for monodomain
  CALL CMISSEquationsTypeInitialise(EquationsM,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetM,EquationsM,Err)
  CALL CMISSEquationsSparsityTypeSet(EquationsM,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(EquationsM,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetM,Err)

  !Create the equations set equations for Finite Elasticity
  CALL CMISSEquationsTypeInitialise(EquationsFE,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetFE,EquationsFE,Err)
  CALL CMISSEquationsSparsityTypeSet(EquationsFE,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(EquationsFE,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetFE,Err)   

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML environment
  CALL CMISSCellMLTypeInitialise(CellML,Err)
  CALL CMISSCellMLCreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import the Shorten et al. 2007 model from a file
!  CALL CMISSCellMLModelImport(CellML,"shorten_ocallaghan_davidson_soboleva_2007.xml",shortenModelIndex,Err)
  CALL CMISSCellMLModelImport(CellML,"shorten_mod_2011_07_04.xml",shortenModelIndex,Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL CMISSCellMLVariableSetAsKnown(CellML,shortenModelIndex,"wal_environment/I_HH",Err)
  !   - to get from the CellML side
!  CALL CMISSCellMLVariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_T",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_ionic_s",Err)
!  CALL CMISSCellMLVariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_ionic_t",Err)
  !   - and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
  !Finish the CellML environment
  CALL CMISSCellMLCreateFinish(CellML,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML <--> OpenCMISS field maps
  CALL CMISSCellMLFieldMapsCreateStart(CellML,Err)
  !Map Vm
  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentFieldM,CMISSFieldUVariableType,1,CMISSFieldValuesSetType, &
    & shortenModelIndex,"wal_environment/vS",CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,shortenModelIndex,"wal_environment/vS",CMISSFieldValuesSetType, &
    & DependentFieldM,CMISSFieldUVariableType,1,CMISSFieldValuesSetType,Err)
! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,shortenModelIndex,"razumova/A_1",CMISSFieldValuesSetType, &
    & IndependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType,Err)
! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  CALL CMISSCellMLFieldMapsCreateFinish(CellML,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Initialise dependent field for monodomain
  !todo - get V_m initialial value.
  CALL CMISSFieldComponentValuesInitialise(DependentFieldM,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,-79.974_CMISSDP,Err)
  
  !Initialise dependent field for Finite Elasticity from undeformed geometry and set hydrostatic pressure
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentFieldFE,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentFieldFE,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentFieldFE,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentFieldFE,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,-8.0_CMISSDP,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML models field
  CALL CMISSFieldTypeInitialise(CellMLModelsField,Err)
  CALL CMISSCellMLModelsFieldCreateStart(CellMLModelsFieldUserNumber,CellML,CellMLModelsField,Err)
  CALL CMISSCellMLModelsFieldCreateFinish(CellML,Err)

  !Create the CellML state field
  CALL CMISSFieldTypeInitialise(CellMLStateField,Err)
  CALL CMISSCellMLStateFieldCreateStart(CellMLStateFieldUserNumber,CellML,CellMLStateField,Err)
  CALL CMISSCellMLStateFieldCreateFinish(CellML,Err)

  !Create the CellML intermediate field
  CALL CMISSFieldTypeInitialise(CellMLIntermediateField,Err)
  CALL CMISSCellMLIntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber,CellML,CellMLIntermediateField,Err)
  CALL CMISSCellMLIntermediateFieldCreateFinish(CellML,Err)
  
  !Create the CellML parameters field
  CALL CMISSFieldTypeInitialise(CellMLParametersField,Err)
  CALL CMISSCellMLParametersFieldCreateStart(CellMLParametersFieldUserNumber,CellML,CellMLParametersField,Err)
  CALL CMISSCellMLParametersFieldCreateFinish(CellML,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Set the Stimulus for monodomain at the left side
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularBottomSurface,BottomSurfaceNodes,BottomNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularLeftSurface,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularRightSurface,RightSurfaceNodes,RightNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularFrontSurface,FrontSurfaceNodes,FrontNormalXi,Err)

  CALL CMISSCellMLFieldComponentGet(CellML,shortenModelIndex,CMISSCellMLParametersFieldType, &
    & "wal_environment/I_HH",stimcomponent,Err)

  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NodeNumber,2,YVALUE, &
        & Err)
        IF(YVALUE .LE. WIDTH/2.0_CMISSDP) CALL CMISSFieldParameterSetUpdateNode(CellMLParametersField,CMISSFieldUVariableType, &
          & CMISSFieldValuesSetType,1,1,NodeNumber,stimcomponent,STIM_VALUE,Err)
    ENDIF
  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
  !Define the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemMultiPhysicsClass,CMISSProblemBioelectricFiniteElasticityType, &
    & CMISSProblemGudunovMonodomainSimpleElasticitySubtype,Err)
  CALL CMISSProblemCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem control loop
  !Loop in time for STIM_STOP with the Stimulus applied.
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)

  !set the main control loop (time loop type)
  CALL CMISSControlLoopTypeInitialise(ControlLoopMain,Err)
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoopMain,Err)
  CALL CMISSControlLoopLabelSet(ControlLoopMain,'MAIN_TIME_LOOP',Err)
  CALL CMISSControlLoopTimesSet(ControlLoopMain,0.0_CMISSDP,STIM_STOP,ELASTICITY_TIME_STEP,Err)
  CALL CMISSControlLoopTimeOutputSet(ControlLoopMain,1,Err)
!  CALL CMISSControlLoopOutputTypeSet(ControlLoopMain,CMISSControlLoopTimingOutput,Err)

  !set the monodomain loop (time loop type)
  CALL CMISSControlLoopTypeInitialise(ControlLoopM,Err)
  CALL CMISSProblemControlLoopGet(Problem,[ControlLoopMonodomainNumber,CMISSControlLoopNode],ControlLoopM,Err)
  CALL CMISSControlLoopLabelSet(ControlLoopM,'MONODOMAIN_TIME_LOOP',Err)
  CALL CMISSControlLoopTimesSet(ControlLoopM,0.0_CMISSDP,ELASTICITY_TIME_STEP,PDE_TIME_STEP,Err)
!  CALL CMISSControlLoopTimesSet(ControlLoopM,0.0_CMISSDP,STIM_STOP,PDE_TIME_STEP,Err)
  CALL CMISSControlLoopTimeOutputSet(ControlLoopM,1,Err)
!  CALL CMISSControlLoopOutputTypeSet(ControlLoopM,CMISSControlLoopTimingOutput,Err)

  !set the finite elasticity loop (simple type)
  CALL CMISSControlLoopTypeInitialise(ControlLoopFE,Err)
  CALL CMISSProblemControlLoopGet(Problem,[ControlLoopElasticityNumber,CMISSControlLoopNode],ControlLoopFE,Err)
  CALL CMISSControlLoopLabelSet(ControlLoopFE,'ELASTICITY_SIMPLE_LOOP',Err)

  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

 
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solvers
  CALL CMISSProblemSolversCreateStart(Problem,Err)

  !Create the DAE solver
  CALL CMISSSolverTypeInitialise(SolverDAE,Err)
  CALL CMISSProblemSolverGet(Problem,[ControlLoopMonodomainNumber,CMISSControlLoopNode], &
    & SolverDAEIndex,SolverDAE,Err)
  !Set the DAE time step to by 10 us
  CALL CMISSSolverDAETimeStepSet(SolverDAE,ODE_TIME_STEP,Err)
  !CALL CMISSSolverDAESolverTypeSet(SolverDAE,CMISSSolverDAEExternal,Err)
  CALL CMISSSolverOutputTypeSet(SolverDAE,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(SolverDAE,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(SolverDAE,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(SolverDAE,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(SolverDAE,CMISSSolverSolverMatrixOutput,Err)

  !Create the parabolic solver
  CALL CMISSSolverTypeInitialise(SolverParabolic,Err)
  CALL CMISSProblemSolverGet(Problem,[ControlLoopMonodomainNumber,CMISSControlLoopNode], &
    & SolverParabolicIndex,SolverParabolic,Err)
  CALL CMISSSolverOutputTypeSet(SolverParabolic,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(SolverParabolic,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(SolverParabolic,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(SolverParabolic,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(SolverParabolic,CMISSSolverSolverMatrixOutput,Err)

  !Create the Finte Elasticity solver
  CALL CMISSSolverTypeInitialise(SolverFE,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverFE,Err)
  CALL CMISSProblemSolverGet(Problem,[ControlLoopElasticityNumber,CMISSControlLoopNode], &
    & SolverFEIndex,SolverFE,Err)
  !CALL CMISSSolverOutputTypeSet(SolverFE,CMISSSolverNoOutput,Err)
  CALL CMISSSolverOutputTypeSet(SolverFE,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(SolverFE,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(SolverFE,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(SolverFE,CMISSSolverSolverMatrixOutput,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(SolverFE,CMISSSolverNewtonJacobianFDCalculated,Err)
  CALL CMISSSolverNewtonLinearSolverGet(SolverFE,LinearSolverFE,Err)
  CALL CMISSSolverLinearTypeSet(LinearSolverFE,CMISSSolverLinearDirectSolveType,Err)

  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver CellML equations
  CALL CMISSProblemCellMLEquationsCreateStart(Problem,Err)

  CALL CMISSSolverTypeInitialise(SolverDAE,Err)
  CALL CMISSCellMLEquationsTypeInitialise(CellMLEquations,Err)
  CALL CMISSProblemSolverGet(Problem,[ControlLoopMonodomainNumber,CMISSControlLoopNode], &
    & SolverDAEIndex,SolverDAE,Err)
  CALL CMISSSolverCellMLEquationsGet(SolverDAE,CellMLEquations,Err)
  CALL CMISSCellMLEquationsCellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  CALL CMISSProblemCellMLEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver equations
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)

  !Create the problem solver parabolic equations (Monodomain)
  CALL CMISSSolverTypeInitialise(SolverParabolic,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsM,Err)
  CALL CMISSProblemSolverGet(Problem,[ControlLoopMonodomainNumber,CMISSControlLoopNode], &
    & SolverParabolicIndex,SolverParabolic,Err)
  CALL CMISSSolverSolverEquationsGet(SolverParabolic,SolverEquationsM,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsM,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsM,CMISSSolverEquationsFullMatrices,Err)  
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsM,EquationsSetM,EquationsSetIndexM,Err)

  !Create the problem solver Finite Elasticity equations
  CALL CMISSSolverTypeInitialise(SolverFE,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsFE,Err)
  CALL CMISSProblemSolverGet(Problem,[ControlLoopElasticityNumber,CMISSControlLoopNode], &
    & SolverFEIndex,SolverFE,Err)
  CALL CMISSSolverSolverEquationsGet(SolverFE,SolverEquationsFE,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsFE,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsFE,CMISSSolverEquationsFullMatrices,Err)  
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsFE,EquationsSetFE,EquationsSetIndexFE,Err)

  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !boundary conditions

  !Prescribe boundary conditions for monodomain
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsM,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsM,BoundaryConditionsM,Err)
  FirstNodeNumber=1
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      LastNodeNumber=(NumberGlobalXElements+1)
    ELSE
      LastNodeNumber=(NumberGlobalXElements+1)*(NumberGlobalYElements+1)
    ENDIF
  ELSE
    LastNodeNumber=(NumberGlobalXElements+1)*(NumberGlobalYElements+1)*(NumberGlobalZElements+1)
  ENDIF
  CALL CMISSDecompositionNodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL CMISSDecompositionNodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  !IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,1,FirstNodeNumber,1, &
    !  & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  !ENDIF
  !IF(LastNodeDomain==ComputationalNodeNumber) THEN
    !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,1,LastNodeNumber,1, &
    !  & CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  !ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsM,Err)

  !Prescribe boundary conditions for Finite Elasticity (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsFE,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsFE,BoundaryConditionsFE,Err)
  !Set x=0 nodes to no x displacment in x
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsFE,DependentFieldFE,CMISSFieldUVariableType,1,1,NodeNumber,1, &
        & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  !Set x=WIDTH nodes to 10% x displacement
  DO node_idx=1,SIZE(RightSurfaceNodes,1)
    NodeNumber=RightSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsFE,DependentFieldFE,CMISSFieldUVariableType,1,1,NodeNumber,1, &
        & CMISSBoundaryConditionFixed,1.0_CMISSDP*LENGTH,Err)
    ENDIF
  ENDDO
  !Set y=0 nodes to no y displacement
  DO node_idx=1,SIZE(FrontSurfaceNodes,1)
    NodeNumber=FrontSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsFE,DependentFieldFE,CMISSFieldUVariableType,1,1,NodeNumber,2, &
        & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  !Set z=0 nodes to no z displacement
  DO node_idx=1,SIZE(BottomSurfaceNodes,1)
    NodeNumber=BottomSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsFE,DependentFieldFE,CMISSFieldUVariableType,1,1,NodeNumber,3, &
        & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem for the first STIM_STOP
  CALL CMISSProblemSolve(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Now turn the stimulus off
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NodeNumber,2,YVALUE, &
        & Err)
        IF(YVALUE .LE. WIDTH/2.0) CALL CMISSFieldParameterSetUpdateNode(CellMLParametersField,CMISSFieldUVariableType, &
          & CMISSFieldValuesSetType,1,1,NodeNumber,stimcomponent,0.0_CMISSDP,Err)
    ENDIF
  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
  !Set the time loop from STIM_STOP to TIME_STOP
  CALL CMISSControlLoopTimesSet(ControlLoopMain,STIM_STOP,TIME_STOP,PDE_TIME_STEP,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem for the rest of the time
!  CALL CMISSProblemSolve(Problem,Err)
! TODO TODO TODO TODO TODO TODO TODO 
  
  !--------------------------------------------------------------------------------------------------------------------------------
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"GudunovMonodomainElasticitySameMeshExample","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"GudunovMonodomainElasticitySameMeshExample","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
  ENDIF
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM GODUNOVMONODOMAINELASTICITYSAMEMESHEXAMPLE
