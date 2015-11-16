!> \file
!> \author Ishani Roy
!> \brief This is an example program to solve an embedded meshing problem on a simple cube
!>        Currently some of the pre processing information is assumed, which will be automated
!>        in future. This example is only to check whether data is being transferred correctly 
!>        between meshes. Trying to use data projection routines to compute Gauss point maps. 
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


!> Main program
PROGRAM EMBEDDEDMESHEXAMPLE

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Test program parameters

  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP
 
  
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionOneUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: RegionTwoUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: BasisOneUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: BasisTwoUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshOneUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshTwoUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: MeshOneUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MeshTwoUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: DecompositionOneUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: DecompositionTwoUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldOneUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldTwoUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldOneUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldTwoUserNumber=16
  !INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  !INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponents=3
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
    & INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS_2,NUMBER_GLOBAL_Y_ELEMENTS_2,NUMBER_GLOBAL_Z_ELEMENTS_2, &
    & INTERPOLATION_TYPE_2,NUMBER_OF_GAUSS_XI_2
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT,Filename
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: Err

   !Other variables
  INTEGER(CMISSIntg) :: NumberofNodes1,NumberofNodes2,NumberOfElements,dim_idx,node_idx, &
    & NumberOfComponents,elem_idx,elem_idx2,NodeCount,gauss_idx,ineach,counts
  INTEGER(CMISSIntg) :: node_idx2,NumberOfElements2,GaussPointNumber,NGP

  !Allocatable
   INTEGER(CMISSIntg), DIMENSION(:), ALLOCATABLE :: NodeNumbers, NodeNumbers1,NodeNumbers2,ElemArray
   INTEGER(CMISSIntg), DIMENSION(:,:), ALLOCATABLE :: NodeElem
   REAL(CMISSRP), DIMENSION(:), ALLOCATABLE :: C1,C2, GC
   REAL(CMISSRP), ALLOCATABLE :: ChildXiCoords(:,:),ParentXiCoords(:,:), Coords(:,:)
   REAL(CMISSRP) :: X, Y, Z, FieldValue,value
  
!====================================================================================
  !Test variables for data projection
  REAL(CMISSRP) :: AbsoluteToleranceSet=1.0E-10_CMISSRP !default is 1.0E-8
  REAL(CMISSRP) :: RelativeToleranceSet=1.0E-6_CMISSRP !default is 1.0E-8
  INTEGER(CMISSIntg) :: MaximumNumberOfIterationsSet=30 !default is 25
  REAL(CMISSRP) :: MaximumIterationUpdateSet=0.4_CMISSRP !default is 0.5
  INTEGER(CMISSIntg) :: NumberOfClosestElementsSet=3 !default is 2/4/8 for 1/2/3 dimensional projection 
  INTEGER(CMISSIntg) :: ProjectionTypeSet=CMFE_DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE !same as default
  REAL(CMISSRP) :: StartingXiSet(3)=[0.5_CMISSRP,0.5_CMISSRP,0.5_CMISSRP] !default is 0.5
  REAL(CMISSRP) :: AbsoluteToleranceGet
  REAL(CMISSRP) :: RelativeToleranceGet
  INTEGER(CMISSIntg) :: MaximumNumberOfIterationsGet
  REAL(CMISSRP) :: MaximumIterationUpdateGet
  INTEGER(CMISSIntg) :: NumberOfClosestElementsGet
  INTEGER(CMISSIntg) :: ProjectionTypeGet
  REAL(CMISSRP), ALLOCATABLE :: StartingXiGet(:)

  REAL(CMISSRP), DIMENSION(2,3) :: DataPointValues!(number_of_data_points,dimension)
  REAL(CMISSRP), DIMENSION(8) :: ElementNumbers!(number_of_data_points,dimension)
  INTEGER(CMISSIntg) :: np

  !CMISS variables

  TYPE(cmfe_BasisType) :: Basis1,Basis2
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition1,Decomposition2
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField1,GeometricField2,DependentField1,DependentField2
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh  
  TYPE(cmfe_MeshType) :: FirstMesh,SecondMesh
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region1,Region2,WorldRegion
  TYPE(cmfe_SolverType) :: Solver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_MeshEmbeddingType) :: MeshEmbedding
  !TYPE(cmfe_EmbeddingXiType) :: EmbeddingXi

  NUMBER_GLOBAL_X_ELEMENTS=1
  NUMBER_GLOBAL_Y_ELEMENTS=1
  NUMBER_GLOBAL_Z_ELEMENTS=1
  INTERPOLATION_TYPE=1
  
  NUMBER_GLOBAL_X_ELEMENTS_2=2
  NUMBER_GLOBAL_Y_ELEMENTS_2=2
  NUMBER_GLOBAL_Z_ELEMENTS_2=2
  INTERPOLATION_TYPE_2=1

  NumberOfComponents=FieldDependentNumberOfComponents
 !======================================
  !Intialise data points

!   DataPointValues(1,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP]
!   DataPointValues(2,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP]  
!   DataPointValues(3,:)=[1.0_CMISSRP,1.0_CMISSRP,1.0_CMISSRP]
!   DataPointValues(4,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP]  
!   DataPointValues(5,:)=[1.0_CMISSRP,1.0_CMISSRP,1.0_CMISSRP]
!   DataPointValues(6,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP]  
!   DataPointValues(7,:)=[0.5_CMISSRP,0.5_CMISSRP,0.5_CMISSRP]
!   DataPointValues(8,:)=[0.5_CMISSRP,0.5_CMISSRP,0.0_CMISSRP]  
 DataPointValues(1,:)= [0.5_CMISSRP,0.5_CMISSRP,0.5_CMISSRP]
 DataPointValues(2,:)= [0.8_CMISSRP,0.6_CMISSRP,0.7_CMISSRP]
!================================================================================  
  

  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)
    
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

!!!!!!!!!!!!!!!!!
! First mesh
!!!!!!!!!!!!!!!!!
  
  
  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region1,Err)
  CALL cmfe_Region_CreateStart(RegionOneUserNumber,WorldRegion,Region1,Err)
  CALL cmfe_Region_LabelSet(Region1,"FirstMeshRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region1,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region1,Err)
  
  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_Initialise(Basis1,Err)
  CALL cmfe_Basis_CreateStart(BasisOneUserNumber,Basis1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL cmfe_Basis_NumberOfXiSet(Basis1,2,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis1,[1,1],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis1,[2,2],Err) 
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis1,3,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis1,[1,1,1],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis1,[2,2,2],Err) 
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis1,Err)

  !Start the creation of a generated mesh in the region - the parent mesh
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshOneUserNumber,Region1,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis1,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(FirstMesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshOneUserNumber,FirstMesh,Err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition1,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionOneUserNumber,FirstMesh,Decomposition1,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition1,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition1,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition1,Err)
  
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField1,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldOneUserNumber,Region1,GeometricField1,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField1,Decomposition1,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField1,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField1,Err)

  !Create the dependent field with 1 variables and 1 component
  CALL cmfe_Field_Initialise(DependentField1,Err)
  CALL cmfe_Field_CreateStart(DependentFieldOneUserNumber,Region1,DependentField1,Err)
  CALL cmfe_Field_TypeSet(DependentField1,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentField1,Decomposition1,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField1,GeometricField1,Err)
  CALL cmfe_Field_DependentTypeSet(DependentField1,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentField1,FieldDependentNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !CALL cmfe_Field_ScalingTypeSet(DependentField,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_CreateFinish(DependentField1,Err)

  !Initialise the field with an initial guess
  CALL cmfe_Field_ComponentValuesInitialise(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.5_CMISSRP, &
    & Err)

  !Field value with a function x^2 + y^2 + z^2
  DO node_idx=1,(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  CALL cmfe_Field_ParameterSetGetNode(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,1,x,Err)
  CALL cmfe_Field_ParameterSetGetNode(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,2,y,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
  CALL cmfe_Field_ParameterSetGetNode(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,3,z,Err)
  ELSE
  Z = 0
  ENDIF
  FieldValue = x**2+y**2+z**2
  CALL cmfe_Field_ParameterSetUpdateNode(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
    & fieldvalue,Err)

  WRITE(*,*) 'The values at the nodes for field 1',fieldvalue
  ENDDO

   CALL cmfe_Nodes_NumberOfNodesGet(RegionOneUserNumber,NumberOfNodes1,Err)
   WRITe(*,*) 'nn',NumberOfNodes1
  
   DO node_idx=1,NumberofNodes1
   CALL cmfe_Field_ParameterSetGetNode(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
     & value,Err)
   WRITE(*,*) 'The values at the nodes for field 1 before',value
   ENDDO

!!!!!!!!!!!!!!!!!!
! The second mesh
!!!!!!!!!!!!!!!!!!

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region2,Err)
  CALL cmfe_Region_CreateStart(RegionTwoUserNumber,WorldRegion,Region2,Err)
  CALL cmfe_Region_LabelSet(Region2,"SecondMeshRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region2,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region2,Err)
  
  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_Initialise(Basis2,Err)
  CALL cmfe_Basis_CreateStart(BasisTwoUserNumber,Basis2,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS_2==0) THEN
    !Set the basis to be a bilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis2,2,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis2,[1,1],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis2,[3,3],Err) 
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis2,3,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis2,[1,1,1],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis2,[3,3,3],Err) 
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis2,Err)

  !Start the creation of a generated mesh in the region - the child mesh
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshTwoUserNumber,Region2,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis2,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS_2==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS_2,NUMBER_GLOBAL_Y_ELEMENTS_2],Err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS_2,NUMBER_GLOBAL_Y_ELEMENTS_2, &
      & NUMBER_GLOBAL_Z_ELEMENTS_2],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(SecondMesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshTwoUserNumber,SecondMesh,Err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition2,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionTwoUserNumber,SecondMesh,Decomposition2,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition2,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition2,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition2,Err)
  
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField2,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldTwoUserNumber,Region2,GeometricField2,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField2,Decomposition2,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS_2/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField2,Err)
  
  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField2,Err)

 !Create the dependent field with 1 variables and 1 component
  CALL cmfe_Field_Initialise(DependentField2,Err)
  CALL cmfe_Field_CreateStart(DependentFieldTwoUserNumber,Region2,DependentField2,Err)
  CALL cmfe_Field_TypeSet(DependentField2,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentField2,Decomposition2,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField2,GeometricField2,Err)
  CALL cmfe_Field_DependentTypeSet(DependentField2,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentField2,FieldDependentNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS_2/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
   !CALL cmfe_Field_ScalingTypeSet(DependentField2,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_CreateFinish(DependentField2,Err)

  !Initialise the field with an initial guess
  CALL cmfe_Field_ComponentValuesInitialise(DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.1_CMISSRP, &
    & Err)

  
   CALL cmfe_Nodes_NumberOfNodesGet(RegionTwoUserNumber,NumberOfNodes2,Err)
   !Before
   DO node_idx=1,NumberOfNodes2
   CALL cmfe_Field_ParameterSetGetNode(DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
     & value,Err)
   WRITE(*,*) 'The values at the nodes for field 2 before',VALUE
 ENDDO
 
 !Start creating an embedded mesh
  CALL cmfe_MeshEmbedding_Initialise(MeshEmbedding,Err)
  
  CALL cmfe_MeshEmbedding_Create(MeshEmbedding, FirstMesh, SecondMesh, Err)
  
  ALLOCATE(ParentXiCoords(NumberofComponents,NumberOfNodes1)) 
  DO node_idx = 1,NumberofNodes1
    DO dim_idx = 1,NumberOfComponents
      CALL cmfe_Field_ParameterSetGetNode(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & node_idx,dim_idx,VALUE,Err)
      ParentXiCoords(dim_idx,node_idx) = VALUE
    ENDDO
  ENDDO
  
  ! Get the positions of the child nodes
  ALLOCATE(ChildXiCoords(NumberofComponents,NumberOfNodes2)) 
  DO node_idx = 1,NumberofNodes2
    DO dim_idx = 1,NumberOfComponents
      CALL cmfe_Field_ParameterSetGetNode(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & node_idx,dim_idx,VALUE,Err)
      ChildXiCoords(dim_idx,node_idx) = VALUE
    ENDDO
  ENDDO

    
  ALLOCATE(NodeNumbers(NumberofNodes1))
  ALLOCATE(NodeNumbers1(NumberofNodes1))
  ALLOCATE(NodeNumbers2(NumberofNodes2))
  
  CALL cmfe_Mesh_NumberOfElementsGet(RegionOneUserNumber,MeshOneUserNumber,NumberOfElements,Err)
  CALL cmfe_Mesh_NumberOfElementsGet(RegionTwoUserNumber,MeshTwoUserNumber,NumberOfElements2,Err)
  !Preprocessing step:

!   By each element 
    counts = 1
    Allocate(ElemArray(NumberOfElements2))
    Allocate(NodeElem(NumberOfElements,NumberofNodes2))
    DO elem_idx = 1,NumberOfElements
     CALL cmfe_MeshElements_NodesGet(RegionOneUserNumber,MeshOneUserNumber,1,elem_idx, &
      & NodeNumbers1,Err)
      DO elem_idx2 =1,NumberOfElements2
      CALL cmfe_MeshElements_NodesGet(RegionTwoUserNumber,MeshTwoUserNumber,1,elem_idx2, &
       & NodeNumbers,Err)      
      DO node_idx=1,size(NodeNumbers)
      DO dim_idx = 1,NumberOfComponents   
      !For testing only: Mappings needs to be put in
      IF(ChildXiCoords(dim_idx,NodeNumbers(node_idx)).le.ParentXiCoords(dim_idx,NodeNumbers1(node_idx)))THEN
!     Get the element number for the parent elem_idx and child elem_idx2       
      NodeElem(elem_idx,NodeNumbers(node_idx))=NodeNumbers(node_idx)    
      ElemArray(elem_idx2)=elem_idx 
      ENDIF
      ENDDO      
      ENDDO     
    ENDDO
  ENDDO

   DO elem_idx = 1,NumberOfElements
    DO node_idx = 1,NumberofNodes2
    NodeNumbers2(node_idx) = NodeElem(elem_idx,node_idx)
    ENDDO
    CALL cmfe_MeshEmbedding_SetChildNodePosition(MeshEmbedding,elem_idx, NodeNumbers2,ChildXiCoords, Err)
   ENDDO


   !Data field from one mesh to the other
    CALL cmfe_MeshEmbedding_PushData(MeshEmbedding, DependentField1,1,DependentField2, 1, Err)
 
   !After
   DO node_idx=1,NumberofNodes2
   CALL cmfe_Field_ParameterSetGetNode(DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
     & value,Err)
   WRITE(*,*) 'The values at the nodes for field 2 after',value
   ENDDO

   !Setting position of Gauss point of parent element wrt child coordinates:Mapping needs to be put in
   Allocate(C1(NumberofNodes1))
   ineach=NumberofNodes1/NumberOfElements2
   Allocate(C2(NumberofNodes1/NumberOfElements2))

   DO elem_idx2 = 1,NumberOfElements2
   elem_idx = ElemArray(elem_idx2)    
      C2=C1((elem_idx2-1)*ineach+1:(elem_idx2)*ineach)
      DO gauss_idx = (elem_idx2-1)*ineach+1,(elem_idx2)*ineach
      !DO gauss_idx =1,NGP
      CALL cmfe_MeshEmbedding_SetGaussPointData(MeshEmbedding,elem_idx,gauss_idx, &
        & C1,elem_idx2,C2, Err)
      ENDDO
   ENDDO
   
   !Move field data

    !CALL cmfe_MeshEmbedding_PullGaussPointData(MeshEmbedding, DependentField1, 1,DependentField2, 1, Err)


   !After
   !DO gauss_idx=1,1
   !CALL cmfe_Field_ParameterSetGetGaussPoint(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,gauss_idx,1,1,value,Err)
   !WRITE(*,*) 'The values at the nodes for field 1 after',value
   !ENDDO

  !================================================Gauss point map==================================
    DO dim_idx = 1,3
       CALL cmfe_MeshEmbedding_GetGaussPointCoord(MeshEmbedding,dim_idx,NGP,C1,Err)        
    DO gauss_idx =1,NGP
        WRITE(*,*) 'data point', datapointvalues(gauss_idx,dim_idx),gauss_idx
        datapointvalues(gauss_idx,dim_idx) = C1(gauss_idx)
        WRITE(*,*) 'data points after', datapointvalues(gauss_idx,dim_idx)
      ENDDO   
    ENDDO
!    
  !SUBROUTINE gauss_point_mapping(datapointvalues,elementnumbers)
   
!   !=========================================================================================================================
!   !Create Data Points and set the values
!   CALL cmfe_DataPoints_CreateStart(RegionTwoUserNumber,SIZE(DataPointValues,1),Err)
!   DO np=1,SIZE(DataPointValues,1)
!     CALL cmfe_DataPoints_ValuesSet(RegionTwoUserNumber,np,DataPointValues(np,:),Err)     
!   ENDDO
!   CALL cmfe_DataPoints_CreateFinish(RegionTwoUserNumber,Err)  
!   !=========================================================================================================================
!   !=========================================================================================================================
!   !Create a data projection
!   CALL cmfe_DataProjection_CreateStart(RegionTwoUserNumber,DependentFieldTwoUserNumber,RegionTwoUserNumber,Err)
!   !=========================================================================================================================
!   !Test parameter set functions
!   CALL cmfe_DataProjection_AbsoluteToleranceSet(RegionTwoUserNumber,AbsoluteToleranceSet,Err) !test
!   CALL cmfe_DataProjection_MaximumIterationUpdateSet(RegionTwoUserNumber,MaximumIterationUpdateSet,Err) !test
!   CALL cmfe_DataProjection_MaximumNumberOfIterationsSet(RegionTwoUserNumber,MaximumNumberOfIterationsSet,Err) !test
!   CALL cmfe_DataProjection_NumberOfClosestElementsSet(RegionTwoUserNumber,NumberOfClosestElementsSet,Err) !test
!   CALL cmfe_DataProjection_ProjectionTypeSet(RegionTwoUserNumber,ProjectionTypeSet,Err)
!   CALL cmfe_DataProjection_RelativeToleranceSet(RegionTwoUserNumber,RelativeToleranceSet,Err) !test
!   CALL cmfe_DataProjection_StartingXiSet(RegionTwoUserNumber,StartingXiSet,Err) !test
!   !=========================================================================================================================
!   !Finish data projection  
!   CALL cmfe_DataProjection_CreateFinish(RegionTwoUserNumber,Err)
!   !=========================================================================================================================
!   !Test parameter get functions
!   CALL cmfe_DataProjection_AbsoluteToleranceGet(RegionTwoUserNumber,AbsoluteToleranceGet,Err) !test
!   CALL cmfe_DataProjection_MaximumIterationUpdateGet(RegionTwoUserNumber,MaximumIterationUpdateGet,Err) !test
!   CALL cmfe_DataProjection_MaximumNumberOfIterationsGet(RegionTwoUserNumber,MaximumNumberOfIterationsGet,Err) !test
!   CALL cmfe_DataProjection_NumberOfClosestElementsGet(RegionTwoUserNumber,NumberOfClosestElementsGet,Err) !test
!   CALL cmfe_DataProjection_ProjectionTypeGet(RegionTwoUserNumber,ProjectionTypeGet,Err) !test
!   CALL cmfe_DataProjection_RelativeToleranceGet(RegionTwoUserNumber,RelativeToleranceGet,Err) !test
!   CALL cmfe_DataProjection_StartingXiGet(RegionTwoUserNumber,StartingXiGet,Err) !test !! This is the GP xi  
!   
!   !=========================================================================================================================
!   !Start data projection
!   CALL cmfe_DataProjection_Evaluate(RegionTwoUserNumber,Err)
!   
!   !Get the closest elementnumber back - in an array
!   !=========================================================================================================================
!   !Destroy used types
!   CALL cmfe_DataProjection_Destroy(RegionTwoUserNumber,Err)
!   CALL cmfe_DataPoints_Destroy(RegiontwoUserNumber,Err)
!     
!   CALL cmfe_Region_Destroy(RegionTwoUserNumber,Err)
  !=========================================================================================================================
  !Create Data Points and set the values
  CALL cmfe_DataPoints_CreateStart(RegionOneUserNumber,SIZE(DataPointValues,1),Err)
  DO np=1,SIZE(DataPointValues,1)
    CALL cmfe_DataPoints_ValuesSet(RegionOneUserNumber,np,DataPointValues(np,:),Err)     
  ENDDO
  CALL cmfe_DataPoints_CreateFinish(RegionOneUserNumber,Err)  
  !=========================================================================================================================
  !=========================================================================================================================
  !Create a data projection
  CALL cmfe_DataProjection_CreateStart(RegionOneUserNumber,DependentFieldOneUserNumber,RegionOneUserNumber,Err)
  !=========================================================================================================================
  !Test parameter set functions
  CALL cmfe_DataProjection_AbsoluteToleranceSet(RegionOneUserNumber,AbsoluteToleranceSet,Err) !test
  CALL cmfe_DataProjection_MaximumIterationUpdateSet(RegionOneUserNumber,MaximumIterationUpdateSet,Err) !test
  CALL cmfe_DataProjection_MaximumNumberOfIterationsSet(RegionOneUserNumber,MaximumNumberOfIterationsSet,Err) !test
  CALL cmfe_DataProjection_NumberOfClosestElementsSet(RegionOneUserNumber,NumberOfClosestElementsSet,Err) !test
  CALL cmfe_DataProjection_ProjectionTypeSet(RegionOneUserNumber,ProjectionTypeSet,Err)
  CALL cmfe_DataProjection_RelativeToleranceSet(RegionOneUserNumber,RelativeToleranceSet,Err) !test
  CALL cmfe_DataProjection_StartingXiSet(RegionOneUserNumber,StartingXiSet,Err) !test
  !=========================================================================================================================
  !Finish data projection  
  CALL cmfe_DataProjection_CreateFinish(RegionOneUserNumber,Err)
  !=========================================================================================================================
  !Test parameter get functions
  CALL cmfe_DataProjection_AbsoluteToleranceGet(RegionOneUserNumber,AbsoluteToleranceGet,Err) !test
  CALL cmfe_DataProjection_MaximumIterationUpdateGet(RegionOneUserNumber,MaximumIterationUpdateGet,Err) !test
  CALL cmfe_DataProjection_MaximumNumberOfIterationsGet(RegionOneUserNumber,MaximumNumberOfIterationsGet,Err) !test
  CALL cmfe_DataProjection_NumberOfClosestElementsGet(RegionOneUserNumber,NumberOfClosestElementsGet,Err) !test
  CALL cmfe_DataProjection_ProjectionTypeGet(RegionOneUserNumber,ProjectionTypeGet,Err) !test
  CALL cmfe_DataProjection_RelativeToleranceGet(RegionOneUserNumber,RelativeToleranceGet,Err) !test
  CALL cmfe_DataProjection_StartingXiGet(RegionOneUserNumber,StartingXiGet,Err) !test !! This is the GP xi  
  
  !=========================================================================================================================
  !Start data projection
  CALL cmfe_DataProjection_Evaluate(RegionOneUserNumber,Err)
  
  !Get the closest elementnumber back - in an array
  !=========================================================================================================================
  !Destroy used types
  CALL cmfe_DataProjection_Destroy(RegionOneUserNumber,Err)
  CALL cmfe_DataPoints_Destroy(RegionOneUserNumber,Err)
    
  CALL cmfe_Region_Destroy(RegionOneUserNumber,Err)
  CALL cmfe_CoordinateSystem_Destroy(CoordinateSystemUserNumber,Err)  
  
  !=========================================================================================================================

  !=================================================================================================
  !END SUBROUTINE gauss_point_mapping

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM EMBEDDEDMESHEXAMPLE
