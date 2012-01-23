!> \file
!> \author Ishani Roy
!> \brief This is an example program to solve an embedded meshing problem on a simple cube
!>        Currently some of the pre processing information is assumed, which will be automated
!>        in future. This example is only to check whether data is being transferred correctly 
!>        between meshes. Currently it is a constant field value being transferred. 
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

  USE OPENCMISS
  USE MPI

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  
  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=1.0_CMISSDP
 
  
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
  !INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldOneUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldTwoUserNumber=16
  !INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  !INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponents=2
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
    & NumberOfComponents,elem_idx,elem_idx2,NodeCount,gauss_idx,ineach,counts!,NodeNumbers(8),NodeNumbers1(8)
!  INTEGER(CMISSIntg) :: NodeElem(2,8)!,ElemArray(2)
!  REAL(CMISSDP) :: ChildXiCoords(3,12),ParentXiCoords(3,8)
!  REAL(CMISSDP) :: C1(8),C2(4)
  INTEGER(CMISSIntg) :: node_idx2,NumberOfElements2,GaussPointNumber,NGP

  !Allocatable
   INTEGER(CMISSIntg), DIMENSION(:), ALLOCATABLE :: NodeNumbers, NodeNumbers1,NodeNumbers2,ElemArray
   INTEGER(CMISSIntg), DIMENSION(:,:), ALLOCATABLE :: NodeElem
   REAL(CMISSDP), DIMENSION(:), ALLOCATABLE :: C1,C2
   REAL(CMISSDP), ALLOCATABLE :: ChildXiCoords(:,:),ParentXiCoords(:,:), Coords(:,:)
!   REAL(CMISSDP), POINTER :: GEOMETRIC_PARAMETERS(:)
   REAL(CMISSDP) :: X, Y, Z, FieldValue,value
  !CMISS variables

  TYPE(CMISSBasisType) :: Basis1,Basis2
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition1,Decomposition2
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField1,GeometricField2,DependentField1,DependentField2
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: FirstMesh,SecondMesh
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region1,Region2,WorldRegion
  TYPE(CMISSSolverType) :: Solver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSMeshEmbeddingType) :: MeshEmbedding
  !TYPE(CMISSEmbeddingXiType) :: EmbeddingXi

  NUMBER_GLOBAL_X_ELEMENTS=1
  NUMBER_GLOBAL_Y_ELEMENTS=1
  NUMBER_GLOBAL_Z_ELEMENTS=0
  INTERPOLATION_TYPE=1
  
  NUMBER_GLOBAL_X_ELEMENTS_2=1
  NUMBER_GLOBAL_Y_ELEMENTS_2=2
  NUMBER_GLOBAL_Z_ELEMENTS_2=0
  INTERPOLATION_TYPE_2=1

  NumberOfComponents=FieldDependentNumberOfComponents
 
  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)
    
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

!!!!!!!!!!!!!!!!!
! First mesh
!!!!!!!!!!!!!!!!!
  
  
  !Start the creation of the region
  CALL CMISSRegion_Initialise(Region1,Err)
  CALL CMISSRegion_CreateStart(RegionOneUserNumber,WorldRegion,Region1,Err)
  CALL CMISSRegion_LabelSet(Region1,"FirstMeshRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region1,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region1,Err)
  
  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasis_Initialise(Basis1,Err)
  CALL CMISSBasis_CreateStart(BasisOneUserNumber,Basis1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSBasis_NumberOfXiSet(Basis1,2,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis1,(/1,1/),Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis1,(/3,3/),Err) 
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL CMISSBasis_NumberOfXiSet(Basis1,3,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis1,(/1,1,1/),Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis1,(/3,3,3/),Err) 
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(Basis1,Err)

  !Start the creation of a generated mesh in the region - the parent mesh
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshOneUserNumber,Region1,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Basis1,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(FirstMesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshOneUserNumber,FirstMesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition1,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionOneUserNumber,FirstMesh,Decomposition1,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition1,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition1,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition1,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField1,Err)
  CALL CMISSField_CreateStart(GeometricFieldOneUserNumber,Region1,GeometricField1,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField1,Decomposition1,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSField_ComponentMeshComponentSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField1,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField1,Err)

  !Create the dependent field with 1 variables and 1 component
  CALL CMISSField_Initialise(DependentField1,Err)
  CALL CMISSField_CreateStart(DependentFieldOneUserNumber,Region1,DependentField1,Err)
  CALL CMISSField_TypeSet(DependentField1,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentField1,Decomposition1,Err)
  CALL CMISSField_GeometricFieldSet(DependentField1,GeometricField1,Err)
  CALL CMISSField_DependentTypeSet(DependentField1,CMISS_FIELD_DEPENDENT_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentField1,FieldDependentNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !CALL CMISSField_ScalingTypeSet(DependentField,CMISS_FIELD_UNIT_SCALING,Err)
  CALL CMISSField_CreateFinish(DependentField1,Err)

  !Initialise the field with an initial guess
  CALL CMISSField_ComponentValuesInitialise(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.5_CMISSDP, &
    & Err)

  !Field value with a function x^2 + y^2 + z^2
  DO node_idx=1,(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  CALL CMISSField_ParameterSetGetNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,1,x,Err)
  CALL CMISSField_ParameterSetGetNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,2,y,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
  CALL CMISSField_ParameterSetGetNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,3,z,Err)
  ELSE
  Z = 0
  ENDIF
  FieldValue = x**2+y**2+z**2
  CALL CMISSField_ParameterSetUpdateNode(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
    & fieldvalue,Err)

  WRITE(*,*) 'The values at the nodes for field 1',fieldvalue
  ENDDO

   CALL CMISSNodes_NumberOfNodesGet(RegionOneUserNumber,NumberOfNodes1,Err)
   WRITe(*,*) 'nn',NumberOfNodes1
  
   DO node_idx=1,NumberofNodes1
   CALL CMISSField_ParameterSetGetNode(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
     & value,Err)
   WRITE(*,*) 'The values at the nodes for field 1 before',value
   ENDDO

!!!!!!!!!!!!!!!!!!
! The second mesh
!!!!!!!!!!!!!!!!!!

  !Start the creation of the region
  CALL CMISSRegion_Initialise(Region2,Err)
  CALL CMISSRegion_CreateStart(RegionTwoUserNumber,WorldRegion,Region2,Err)
  CALL CMISSRegion_LabelSet(Region2,"SecondMeshRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region2,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region2,Err)
  
  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasis_Initialise(Basis2,Err)
  CALL CMISSBasis_CreateStart(BasisTwoUserNumber,Basis2,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS_2==0) THEN
    !Set the basis to be a bilinear Lagrange basis
    CALL CMISSBasis_NumberOfXiSet(Basis2,2,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis2,(/1,1/),Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis2,(/3,3/),Err) 
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL CMISSBasis_NumberOfXiSet(Basis2,3,Err)
        CALL CMISSBasis_InterpolationXiSet(Basis2,(/1,1,1/),Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis2,(/3,3,3/),Err) 
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(Basis2,Err)

  !Start the creation of a generated mesh in the region - the child mesh
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshTwoUserNumber,Region2,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Basis2,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS_2==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS_2,NUMBER_GLOBAL_Y_ELEMENTS_2],Err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS_2,NUMBER_GLOBAL_Y_ELEMENTS_2, &
      & NUMBER_GLOBAL_Z_ELEMENTS_2],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(SecondMesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshTwoUserNumber,SecondMesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition2,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionTwoUserNumber,SecondMesh,Decomposition2,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition2,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition2,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition2,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField2,Err)
  CALL CMISSField_CreateStart(GeometricFieldTwoUserNumber,Region2,GeometricField2,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField2,Decomposition2,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS_2/=0) THEN
    CALL CMISSField_ComponentMeshComponentSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField2,Err)
  
  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField2,Err)

 !Create the dependent field with 1 variables and 1 component
  CALL CMISSField_Initialise(DependentField2,Err)
  CALL CMISSField_CreateStart(DependentFieldTwoUserNumber,Region2,DependentField2,Err)
  CALL CMISSField_TypeSet(DependentField2,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentField2,Decomposition2,Err)
  CALL CMISSField_GeometricFieldSet(DependentField2,GeometricField2,Err)
  CALL CMISSField_DependentTypeSet(DependentField2,CMISS_FIELD_DEPENDENT_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentField2,FieldDependentNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS_2/=0) THEN
    CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
   !CALL CMISSField_ScalingTypeSet(DependentField2,CMISS_FIELD_UNIT_SCALING,Err)
  CALL CMISSField_CreateFinish(DependentField2,Err)

  !Initialise the field with an initial guess
  CALL CMISSField_ComponentValuesInitialise(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.1_CMISSDP, &
    & Err)
 
   CALL CMISSNodes_NumberOfNodesGet(RegionTwoUserNumber,NumberOfNodes2,Err)
   !Before
   DO node_idx=1,NumberOfNodes2
   CALL CMISSField_ParameterSetGetNode(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
     & value,Err)
   WRITE(*,*) 'The values at the nodes for field 2 before',value
   ENDDO

  !Start creating an embedded mesh
  CALL CMISSMeshEmbedding_Initialise(MeshEmbedding,Err)
  
  CALL CMISSMeshEmbedding_Create(MeshEmbedding, FirstMesh, SecondMesh, Err)
  
   Allocate(ParentXiCoords(NumberofComponents,NumberOfNodes1)) 
   DO node_idx = 1,NumberofNodes1
   DO dim_idx = 1,NumberOfComponents
   CALL CMISSField_ParameterSetGetNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
     & node_idx,dim_idx,value,Err)
   ParentXiCoords(dim_idx,node_idx) = value
   ENDDO
   ENDDO

  ! Get the positions of the child nodes
  Allocate(ChildXiCoords(NumberofComponents,NumberOfNodes2)) 
    DO node_idx = 1,NumberofNodes2
    DO dim_idx = 1,NumberOfComponents
    CALL CMISSField_ParameterSetGetNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
      & node_idx,dim_idx,value,Err)
     ChildXiCoords(dim_idx,node_idx) = value
    ENDDO
    ENDDO

    
    Allocate(NodeNumbers(NumberofNodes1))
    Allocate(NodeNumbers1(NumberofNodes1))
    Allocate(NodeNumbers2(NumberofNodes2))

   CALL CMISSMesh_NumberOfElementsGet(RegionOneUserNumber,MeshOneUserNumber,NumberOfElements,Err)
   CALL CMISSMesh_NumberOfElementsGet(RegionTwoUserNumber,MeshTwoUserNumber,NumberOfElements2,Err)
   !Preprocessing step:

!   By each element 
    counts = 1
    Allocate(ElemArray(NumberOfElements2))
    Allocate(NodeElem(NumberOfElements,NumberofNodes2))
    DO elem_idx = 1,NumberOfElements
     CALL CMISSMeshElements_NodesGet(RegionOneUserNumber,MeshOneUserNumber,1,elem_idx, &
      & NodeNumbers1,Err)
      DO elem_idx2 =1,NumberOfElements2
      CALL CMISSMeshElements_NodesGet(RegionTwoUserNumber,MeshTwoUserNumber,1,elem_idx2, &
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
    CALL CMISSMeshEmbedding_SetChildNodePosition(MeshEmbedding,elem_idx, NodeNumbers2,ChildXiCoords, Err)
   ENDDO


   !Data field from one mesh to the other
    CALL CMISSMeshEmbedding_PushData(MeshEmbedding, DependentField1,1,DependentField2, 1, Err)
 
   !After
   DO node_idx=1,NumberofNodes2
   CALL CMISSField_ParameterSetGetNode(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
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
      CALL CMISSMeshEmbedding_SetGaussPointData(MeshEmbedding,elem_idx,gauss_idx, &
        & C1,elem_idx2,C2, Err)
      ENDDO
   ENDDO
   
   !Move field data

    CALL CMISSMeshEmbedding_PullGaussPointData(MeshEmbedding, DependentField1, 1,DependentField2, 1, Err)

   !After
   !DO gauss_idx=1,1
   !CALL CMISSField_ParameterSetGetGaussPoint(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,gauss_idx,1,value,Err)
   !WRITE(*,*) 'The values at the nodes for field 1 after',value
   !ENDDO
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM EMBEDDEDMESHEXAMPLE
