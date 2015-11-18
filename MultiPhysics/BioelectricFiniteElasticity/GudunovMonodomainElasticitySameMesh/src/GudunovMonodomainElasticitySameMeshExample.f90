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

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !--------------------------------------------------------------------------------------------------------------------------------
  !Test program parameters
  REAL(CMISSRP), PARAMETER :: tol=1.0E-8_CMISSRP

  LOGICAL :: independent_field_auto_create=.FALSE.
  LOGICAL :: uniaxial_extension_bc=.FALSE.

  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_ELEMENTS=50
  !all lengths in [cm]
  REAL(CMISSRP), PARAMETER :: LENGTH= 2.5_CMISSRP ! X-direction
  REAL(CMISSRP), PARAMETER :: WIDTH= 0.05_CMISSRP ! Y-direction
  REAL(CMISSRP), PARAMETER :: HEIGHT=0.05_CMISSRP ! Z-direction
  
  !all times in [ms]
  REAL(CMISSRP), PARAMETER :: STIM_STOP=0.05_CMISSRP
  REAL(CMISSRP), PARAMETER :: TIME_STOP=50.0_CMISSRP

  REAL(CMISSRP), PARAMETER :: ODE_TIME_STEP=0.00001_CMISSRP
  REAL(CMISSRP), PARAMETER :: PDE_TIME_STEP=0.0005_CMISSRP
  REAL(CMISSRP), PARAMETER :: ELASTICITY_TIME_STEP=0.01_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: OUTPUT_FREQUENCY=2

  !--------------------------------------------------------------------------------------------------------------------------------
  !stimulation current in [uA/cm^2]
  REAL(CMISSRP), PARAMETER :: STIM_VALUE=8000.0_CMISSRP
  !condctivity in [mS/cm]
  REAL(CMISSRP), PARAMETER :: CONDUCTIVITY=3.828_CMISSRP
  !surface area to volume ratio in [cm^-1]
  REAL(CMISSRP), PARAMETER :: Am=500.0_CMISSRP
  !membrane capacitance in [uF/cm^2]
  REAL(CMISSRP), PARAMETER :: Cm=1.0_CMISSRP

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

  REAL(CMISSRP) :: YVALUE

  INTEGER(CMISSIntg) :: Err

  !CMISS variables

  TYPE(cmfe_BasisType) :: QuadraticBasis,LinearBasis,Basis(2)
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsM,BoundaryConditionsFE
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoopMain
  TYPE(cmfe_ControlLoopType) :: ControlLoopM,ControlLoopFE
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: EquationsM,EquationsFE
  TYPE(cmfe_EquationsSetType) :: EquationsSetM,EquationsSetFE
  TYPE(cmfe_FieldType) :: EquationsSetFieldM,EquationsSetFieldFE
  TYPE(cmfe_FieldType) :: GeometricFieldM,GeometricFieldFE
  TYPE(cmfe_FieldType) :: DependentFieldM,DependentFieldFE
  TYPE(cmfe_FieldType) :: IndependentField
  TYPE(cmfe_FieldType) :: MaterialFieldM,MaterialFieldFE
  TYPE(cmfe_FieldType) :: FibreField
  TYPE(cmfe_FieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh  
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: SolverDAE,SolverParabolic
  TYPE(cmfe_SolverType) :: SolverFE,LinearSolverFE
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsM,SolverEquationsFE

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
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  !Trap errors
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  
  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  CALL cmfe_OutputSetOn("GudunovMonodomainElasticitySameMesh",Err)
  
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
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      !Set the coordinate system to be 1D
      CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,1,Err)
    ELSE
      !Set the coordinate system to be 2D
      CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
    ENDIF
  ELSE
    !Set the coordinate system to be 3D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_LabelSet(Region,"Region",Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the bases
  !Define basis functions - tri-Quadratic Lagrange 
  CALL cmfe_Basis_Initialise(QuadraticBasis,Err)
  CALL cmfe_Basis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL cmfe_Basis_TypeSet(QuadraticBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(QuadraticBasis,NumberOfXiCoordinates,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      CALL cmfe_Basis_InterpolationXiSet(QuadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis,[NumberOfGaussPoints],Err)
    ELSE
      CALL cmfe_Basis_InterpolationXiSet(QuadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
        & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis,[NumberOfGaussPoints,NumberOfGaussPoints],Err)
    ENDIF
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(QuadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
      & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
      & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  ENDIF
  CALL cmfe_Basis_CreateFinish(QuadraticBasis,Err)

  !Define basis functions - tri-Linear Lagrange
  CALL cmfe_Basis_Initialise(LinearBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL cmfe_Basis_TypeSet(LinearBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(LinearBasis,NumberOfXiCoordinates,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      CALL cmfe_Basis_InterpolationXiSet(LinearBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis,[NumberOfGaussPoints],Err)
    ELSE
      CALL cmfe_Basis_InterpolationXiSet(LinearBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
        & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis,[NumberOfGaussPoints,NumberOfGaussPoints],Err)
    ENDIF
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(LinearBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
      & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis, &
      & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  ENDIF
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)

  Basis(1)=QuadraticBasis
  Basis(2)=LinearBasis

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the basis - Basis is an array: [QuadraticBasis,LinearBasis]
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements],Err)
    ELSE
      CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH,WIDTH],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
    ENDIF
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH,WIDTH,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
        & NumberGlobalZElements],Err)
  ENDIF    
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)


  !================================================================================================================================
  !  F I N I T E   E L A S T C I T Y
  !================================================================================================================================
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for finite elasticity - quadratic interpolation
  CALL cmfe_Field_Initialise(GeometricFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumberFE,Region,GeometricFieldFE,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldFE,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricFieldFE,Decomposition,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldFE,CMFE_FIELD_GEOMETRIC_TYPE,Err)  
  CALL cmfe_Field_NumberOfVariablesSet(GeometricFieldFE,FieldGeometryNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  IF(NumberGlobalYElements/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
    ENDIF
  ENDIF
  CALL cmfe_Field_VariableLabelSet(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_CreateFinish(GeometricFieldFE,Err)
  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricFieldFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a fibre field and attach it to the geometric field - quadratic interpolation
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricFieldFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err) 
  IF(NumberGlobalYElements/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
    ENDIF
  ENDIF
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a material field for Finite Elasticity and attach it to the geometric field - quadratic interpolation
  CALL cmfe_Field_Initialise(MaterialFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumberFE,Region,MaterialFieldFE,Err)
  CALL cmfe_Field_TypeSet(MaterialFieldFE,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialFieldFE,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialFieldFE,GeometricFieldFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialFieldFE,FieldMaterialNumberOfVariablesFE,Err)
  CALL cmfe_Field_VariableTypesSet(MaterialFieldFE,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_V_VARIABLE_TYPE],Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponentsFE1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,FieldMaterialNumberOfComponentsFE2,Err)
! default is CMFE_FIELD_NODE_BASED_INTERPOLATION
  CALL cmfe_Field_ComponentMeshComponentSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"MaterialFE",Err)
  CALL cmfe_Field_VariableLabelSet(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,"Gravity",Err)
  CALL cmfe_Field_CreateFinish(MaterialFieldFE,Err)
  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2.0_CMISSRP, &
    & Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,6.0_CMISSRP, &
    & Err)
!  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldFE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP,Err)

  !Create the dependent field for FE with 2 variables and * components 
  !  3-d: 3 displacement (quad interpol), 1 pressure (lin interpol) --> * = 4
  !  2-d: 2 displacement (quad interpol), 1 pressure (lin interpol) --> * = 3
  !  1-d: 1 displacement (quad interpol), 1 pressure (lin interpol) --> * = 2
  CALL cmfe_Field_Initialise(DependentFieldFE,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumberFE,Region,DependentFieldFE,Err)
  CALL cmfe_Field_TypeSet(DependentFieldFE,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentFieldFE,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentFieldFE,GeometricFieldFE,Err)
  CALL cmfe_Field_DependentTypeSet(DependentFieldFE,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentFieldFE,FieldDependentNumberOfVariablesFE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponentsFE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponentsFE,Err)
  IF(NumberGlobalZElements==0) THEN
    IF(NumberGlobalYElements==0) THEN
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,LinearMeshComponentNumber,Err)
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1, &
        & QuadraticMeshComponentNumber, &
        & Err)
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,LinearMeshComponentNumber,Err)
    ELSE
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,LinearMeshComponentNumber,Err)
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1, &
        & QuadraticMeshComponentNumber, &
        & Err)
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2, &
        & QuadraticMeshComponentNumber, &
        & Err)
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,LinearMeshComponentNumber,Err)
    ENDIF
  ELSE
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber, &
      & Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber, &
      & Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber, &
      & Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  ENDIF
!  CALL cmfe_Field_ScalingTypeSet(DependentFieldFE,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,"DependentFE",Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"Reaction_Force",Err)
  CALL cmfe_Field_CreateFinish(DependentFieldFE,Err)



  !================================================================================================================================
  !  M O N O D O M A I N
  !================================================================================================================================
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for monodomain - quadratic interpolation
  CALL cmfe_Field_Initialise(GeometricFieldM,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumberM,Region,GeometricFieldM,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldM,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricFieldM,Decomposition,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldM,CMFE_FIELD_GEOMETRIC_TYPE,Err)  
  CALL cmfe_Field_NumberOfVariablesSet(GeometricFieldM,FieldGeometryNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  IF(NumberGlobalYElements/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
    ENDIF
  ENDIF
  CALL cmfe_Field_VariableLabelSet(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"GeometryM",Err)
  CALL cmfe_Field_CreateFinish(GeometricFieldM,Err)
  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricFieldM,Err)
        

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a materials field for monodomain and attach it to the geometric field - constant interpolation
  CALL cmfe_Field_Initialise(MaterialFieldM,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumberM,Region,MaterialFieldM,Err)
  CALL cmfe_Field_TypeSet(MaterialFieldM,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialFieldM,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialFieldM,GeometricFieldM,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialFieldM,FieldMaterialNumberOfVariablesM,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponentsM,Err)
! default is CMFE_FIELD_NODE_BASED_INTERPOLATION
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  IF(NumberGlobalYElements/=0) THEN
    CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,5,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    ENDIF
  ENDIF
  CALL cmfe_Field_VariableLabelSet(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"MaterialM",Err)
  CALL cmfe_Field_CreateFinish(MaterialFieldM,Err)
  !Set Am
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Am,Err)
  !Set Cm
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Cm,Err)
  !Set Conductivity
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,CONDUCTIVITY, &
    & Err)
  IF(NumberGlobalYElements/=0) THEN
    CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
      & CONDUCTIVITY,Err)
    IF(NumberGlobalZElements/=0) THEN
      CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5, &
        & CONDUCTIVITY,Err)
    ENDIF
  ENDIF

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the dependent field for monodomain with 2 variables and 1 components 
  CALL cmfe_Field_Initialise(DependentFieldM,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumberM,Region,DependentFieldM,Err)
  CALL cmfe_Field_TypeSet(DependentFieldM,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentFieldM,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentFieldM,GeometricFieldM,Err)
  CALL cmfe_Field_DependentTypeSet(DependentFieldM,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentFieldM,FieldDependentNumberOfVariablesM,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponentsM,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponentsM,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_DimensionSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_DimensionSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,"Vm",Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldM,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"dVm/dt",Err)
  CALL cmfe_Field_CreateFinish(DependentFieldM,Err)


  !================================================================================================================================
  !  INDEPENDENT FIELD

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the independent field for the active stress determined in the cell model of the monodomain equation and used in finite elasticity
!  independent_field_auto_create = .TRUE.
  CALL cmfe_Field_Initialise(IndependentField,Err)
  IF(.NOT. independent_field_auto_create) THEN
    CALL cmfe_Field_CreateStart(FieldIndependentUserNumber,Region,IndependentField,Err)
    CALL cmfe_Field_TypeSet(IndependentField,CMFE_FIELD_GENERAL_TYPE,Err)
    CALL cmfe_Field_MeshDecompositionSet(IndependentField,Decomposition,Err)
    CALL cmfe_Field_GeometricFieldSet(IndependentField,GeometricFieldM,Err)
    CALL cmfe_Field_DependentTypeSet(IndependentField,CMFE_FIELD_INDEPENDENT_TYPE,Err)
    CALL cmfe_Field_NumberOfVariablesSet(IndependentField,FieldIndependentNumberOfVariables,Err)
    CALL cmfe_Field_DimensionSet(IndependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(IndependentField,CMFE_FIELD_U_VARIABLE_TYPE,FieldIndependentNumberOfComponents,Err)
    CALL cmfe_Field_ComponentInterpolationSet(IndependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_NODE_BASED_INTERPOLATION, &
      & Err)
    CALL cmfe_Field_ComponentMeshComponentSet(IndependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
    CALL cmfe_Field_VariableLabelSet(IndependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Acive_Stress",Err)
    CALL cmfe_Field_CreateFinish(IndependentField,Err)
  ENDIF


  !================================================================================================================================
  !  EQUATIONS SET

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for Finite Elasticity
  CALL cmfe_Field_Initialise(EquationsSetFieldFE,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetFE,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberFE,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE], &
    & EquationsSetFieldUserNumberFE,EquationsSetFieldFE,EquationsSetFE,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetFE,Err)

  !Create the equations set dependent field variables for Finite Elasticity
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetFE,FieldDependentUserNumberFE,DependentFieldFE,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetFE,Err)

  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetFE,FieldIndependentUserNumber,IndependentField,Err)
  CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetFE,Err)

  !Create the equations set materials field variables for Finite Elasticity
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetFE,FieldMaterialUserNumberFE,MaterialFieldFE,Err)  
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for monodomain
  CALL cmfe_Field_Initialise(EquationsSetFieldM,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetM,Err)
  !Set the equations set to be a Monodomain equations set
  !> \todo solve the monodomain problem on the fibre field rather than on the geometric field: GeometricField <--> FibreField
  CALL cmfe_EquationsSet_CreateStart(EquationsSetsUserNumberM,Region,GeometricFieldM,[CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
    & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_NO_SUBTYPE],EquationsSetFieldUserNumberM, &
    & EquationsSetFieldM,EquationsSetM,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetM,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set dependent field variables for monodomain
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetM,FieldDependentUserNumberM,DependentFieldM,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetM,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set materials field variables for monodomain
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetM,FieldMaterialUserNumberM,MaterialFieldM,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetM,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set equations for monodomain
  CALL cmfe_Equations_Initialise(EquationsM,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetM,EquationsM,Err)
  CALL cmfe_Equations_SparsityTypeSet(EquationsM,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(EquationsM,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetM,Err)

  !Create the equations set equations for Finite Elasticity
  CALL cmfe_Equations_Initialise(EquationsFE,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetFE,EquationsFE,Err)
  CALL cmfe_Equations_SparsityTypeSet(EquationsFE,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(EquationsFE,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetFE,Err)   

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML environment
  CALL cmfe_CellML_Initialise(CellML,Err)
  CALL cmfe_CellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import the Shorten et al. 2007 model from a file
  CALL cmfe_CellML_ModelImport(CellML,"shorten_mod_2011_07_04.xml",shortenModelIndex,Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL cmfe_CellML_VariableSetAsKnown(CellML,shortenModelIndex,"wal_environment/I_HH",Err)
  !   - to get from the CellML side
!  CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_T",Err)
!  CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_ionic_s",Err)
!  CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_ionic_t",Err)
  !
  !NOTE: If an INTERMEDIATE (or ALGEBRAIC) variable should be used in a mapping, it has to be set as known or wanted first!
  !       --> set "razumova/stress" as wanted!
  !       --> no need to set "wal_environment/vS" since all STATE variables are automatically set as wanted! 
  CALL cmfe_CellML_VariableSetAsWanted(CellML,shortenModelIndex,"razumova/stress",Err)
  !   - and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(CellML,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(CellML,Err)
  !Map the transmembrane voltage V_m
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & shortenModelIndex,"wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"wal_environment/vS",CMFE_FIELD_VALUES_SET_TYPE, &
    & DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"razumova/stress",CMFE_FIELD_VALUES_SET_TYPE, &
    & IndependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_FieldMapsCreateFinish(CellML,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Initialise dependent field for monodomain
  !> \todo - get V_m initialial value.
  CALL cmfe_Field_ComponentValuesInitialise(DependentFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & -79.974_CMISSRP,Err)
  
  !Initialise dependent field for Finite Elasticity from undeformed geometry and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
    & -8.0_CMISSRP,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML models field
  CALL cmfe_Field_Initialise(CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateFinish(CellML,Err)

  !Create the CellML state field
  CALL cmfe_Field_Initialise(CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateFinish(CellML,Err)

  !Create the CellML intermediate field
  CALL cmfe_Field_Initialise(CellMLIntermediateField,Err)
  CALL cmfe_CellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  CALL cmfe_CellML_IntermediateFieldCreateFinish(CellML,Err)
  
  !Create the CellML parameters field
  CALL cmfe_Field_Initialise(CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateFinish(CellML,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Set the Stimulus for monodomain at the left side
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi, &
    & Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)

  CALL cmfe_CellML_FieldComponentGet(CellML,shortenModelIndex,CMFE_CELLML_PARAMETERS_FIELD, &
    & "wal_environment/I_HH",stimcomponent,Err)

  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      SELECT CASE(NumberOfSpatialCoordinates)
      CASE(1)
        WRITE(*,'(A)') "Setting the stimulus not implemented for 1D."
        STOP
      CASE(2)
        !simulate the lower half of the nodes on the left surface
        CALL cmfe_Field_ParameterSetGetNode(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
          & NodeNumber, &
          & 2, &
          & YVALUE,Err)
        IF(YVALUE .LE. WIDTH/2.0_CMISSRP) CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField, &
          & CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,STIM_VALUE,Err)
      CASE(3)
        !stimulate all nodes on the left surface
        CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,STIM_VALUE,Err)
      END SELECT
    ENDIF
  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
    & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)

  !set the main control loop (time loop type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopMain,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoopMain,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopMain,'MAIN_TIME_LOOP',Err)
  !Loop in time for STIM_STOP with the Stimulus applied.
  CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,0.0_CMISSRP,STIM_STOP,ELASTICITY_TIME_STEP,Err)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoopMain,OUTPUT_FREQUENCY,Err)
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,CMFE_CONTROL_LOOP_TIMING_OUTPUT,Err)

  !set the monodomain loop (time loop type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopM,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopM,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopM,'MONODOMAIN_TIME_LOOP',Err)
  CALL cmfe_ControlLoop_TimesSet(ControlLoopM,0.0_CMISSRP,ELASTICITY_TIME_STEP,PDE_TIME_STEP,Err)
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopM,CMFE_CONTROL_LOOP_NO_OUTPUT,Err)

  !set the finite elasticity loop (simple type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopFE,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE],ControlLoopFE,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopFE,'ELASTICITY_LOOP',Err)

  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solvers
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)

  !Create the DAE solver
  CALL cmfe_Solver_Initialise(SolverDAE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
    & SolverDAEIndex,SolverDAE,Err)
  CALL cmfe_Solver_DAETimeStepSet(SolverDAE,ODE_TIME_STEP,Err)
  !> \todo - solve the CellML equations on the GPU for efficiency (later)
  !CALL cmfe_Solver_DAESolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_EXTERNAL,Err) 
  CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_MATRIX_OUTPUT,Err)

  !Create the parabolic solver
  CALL cmfe_Solver_Initialise(SolverParabolic,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
    & SolverParabolicIndex,SolverParabolic,Err)
  CALL cmfe_Solver_DynamicSchemeSet(SolverParabolic,CMFE_SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME,Err)
  CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_MATRIX_OUTPUT,Err)

  !Create the Finte Elasticity solver
  CALL cmfe_Solver_Initialise(SolverFE,Err)
  CALL cmfe_Solver_Initialise(LinearSolverFE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE], &
    & SolverFEIndex,SolverFE,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_NO_OUTPUT,Err)
  CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(SolverFE,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(SolverFE,LinearSolverFE,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolverFE,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)

  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)

  CALL cmfe_Solver_Initialise(SolverDAE,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
    & SolverDAEIndex,SolverDAE,Err)
  CALL cmfe_Solver_CellMLEquationsGet(SolverDAE,CellMLEquations,Err)
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)

  !Create the problem solver parabolic equations (Monodomain)
  CALL cmfe_Solver_Initialise(SolverParabolic,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsM,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
    & SolverParabolicIndex,SolverParabolic,Err)
  CALL cmfe_Solver_SolverEquationsGet(SolverParabolic,SolverEquationsM,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsM,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsM,CMFE_SOLVER_FULL_MATRICES,Err)  
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsM,EquationsSetM,EquationsSetIndexM,Err)

  !Create the problem solver Finite Elasticity equations
  CALL cmfe_Solver_Initialise(SolverFE,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsFE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE], &
    & SolverFEIndex,SolverFE,Err)
  CALL cmfe_Solver_SolverEquationsGet(SolverFE,SolverEquationsFE,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsFE,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsFE,CMFE_SOLVER_FULL_MATRICES,Err)  
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsFE,EquationsSetFE,EquationsSetIndexFE,Err)

  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !boundary conditions

  !Prescribe boundary conditions for monodomain
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsM,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsM,BoundaryConditionsM,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsM,Err)

  !Prescribe boundary conditions for Finite Elasticity (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsFE,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsFE,BoundaryConditionsFE,Err)

  SELECT CASE(NumberOfSpatialCoordinates)
  CASE(3)

    !  uniaxial_extension_bc=.TRUE.
    IF(uniaxial_extension_bc) THEN
      !Set x=0 nodes to no x displacment in x
      DO node_idx=1,SIZE(LeftSurfaceNodes,1)
        NodeNumber=LeftSurfaceNodes(node_idx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
            & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
        ENDIF
      ENDDO
      !Set x=WIDTH nodes to 10% x displacement
      DO node_idx=1,SIZE(RightSurfaceNodes,1)
        NodeNumber=RightSurfaceNodes(node_idx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
            & CMFE_BOUNDARY_CONDITION_FIXED,1.1_CMISSRP*LENGTH,Err)
        ENDIF
      ENDDO
      !Set y=0 nodes to no y displacement
      DO node_idx=1,SIZE(FrontSurfaceNodes,1)
        NodeNumber=FrontSurfaceNodes(node_idx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
            & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
        ENDIF
      ENDDO
      !Set z=0 nodes to no z displacement
      DO node_idx=1,SIZE(BottomSurfaceNodes,1)
        NodeNumber=BottomSurfaceNodes(node_idx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
            & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
        ENDIF
      ENDDO

    ELSE !uniaxial_extension_bc

      DO node_idx=1,SIZE(LeftSurfaceNodes,1)
        NodeNumber=LeftSurfaceNodes(node_idx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          DO ComponentNumber=2,3
            CALL cmfe_Field_ParameterSetGetNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
              & NodeNumber, &
              & ComponentNumber,YVALUE,Err)
            IF(YVALUE<tol) THEN
              !fix all nodes (0,0,*) in x- and y-direction
              !fix all nodes (0,*,0) in x- and z-direction
              CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
                & NodeNumber, &
                & 1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
              CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
                & NodeNumber, &
                & ComponentNumber,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      !fix all (*,0,0) in y- and z-direction
      DO node_idx=1,SIZE(BottomSurfaceNodes,1)
        NodeNumber=BottomSurfaceNodes(node_idx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL cmfe_Field_ParameterSetGetNode(GeometricFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
            & NodeNumber, &
            & 2,YVALUE,Err)
          IF(YVALUE<tol) THEN
            CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
              & 2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
            CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
              & 3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
          ENDIF
        ENDIF
      ENDDO
      
    ENDIF !uniaxial_extension_bc
    
  CASE DEFAULT
    WRITE(*,'(A)') "Boundary conditions not implemented for 1D and 2D."
    STOP
  END SELECT
  
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem for the first STIM_STOP
  CALL cmfe_Problem_Solve(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Now turn the stimulus off
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      SELECT CASE(NumberOfSpatialCoordinates)
      CASE(1)
        WRITE(*,'(A)') "Setting the stimulus not implemented for 1D."
        STOP
      CASE(2)
        CALL cmfe_Field_ParameterSetGetNode(GeometricFieldM,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
          & NodeNumber, &
          & 2, &
          & YVALUE,Err)
        IF(YVALUE .LE. WIDTH/2.0_CMISSRP) CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField, &
          & CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,0.0_CMISSRP,Err)
      CASE(3)
        CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,0.0_CMISSRP,Err)
      END SELECT
    ENDIF
  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
  !Set the time loop from STIM_STOP to TIME_STOP
  CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,STIM_STOP,TIME_STOP,ELASTICITY_TIME_STEP,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem for the rest of the time
  CALL cmfe_Problem_Solve(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"GudunovMonodomainElasticitySameMeshExample","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"GudunovMonodomainElasticitySameMeshExample","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM GODUNOVMONODOMAINELASTICITYSAMEMESHEXAMPLE
