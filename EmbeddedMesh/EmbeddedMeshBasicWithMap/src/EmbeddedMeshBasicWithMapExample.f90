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
   REAL(CMISSDP), DIMENSION(:), ALLOCATABLE :: C1,C2, GC
   REAL(CMISSDP), ALLOCATABLE :: ChildXiCoords(:,:),ParentXiCoords(:,:), Coords(:,:)
   REAL(CMISSDP) :: X, Y, Z, FieldValue,value
  
!====================================================================================
  !Test variables for data projection
  REAL(CMISSDP) :: AbsoluteToleranceSet=1.0E-10_CMISSDP !default is 1.0E-8
  REAL(CMISSDP) :: RelativeToleranceSet=1.0E-6_CMISSDP !default is 1.0E-8
  INTEGER(CMISSIntg) :: MaximumNumberOfIterationsSet=30 !default is 25
  REAL(CMISSDP) :: MaximumIterationUpdateSet=0.4_CMISSDP !default is 0.5
  INTEGER(CMISSIntg) :: NumberOfClosestElementsSet=3 !default is 2/4/8 for 1/2/3 dimensional projection 
  INTEGER(CMISSIntg) :: ProjectionTypeSet=CMISSDataProjectionAllElementsProjectionType !same as default
  REAL(CMISSDP) :: StartingXiSet(3)=(/0.5_CMISSDP,0.5_CMISSDP,0.5_CMISSDP/) !default is 0.5
  REAL(CMISSDP) :: AbsoluteToleranceGet
  REAL(CMISSDP) :: RelativeToleranceGet
  INTEGER(CMISSIntg) :: MaximumNumberOfIterationsGet
  REAL(CMISSDP) :: MaximumIterationUpdateGet
  INTEGER(CMISSIntg) :: NumberOfClosestElementsGet
  INTEGER(CMISSIntg) :: ProjectionTypeGet
  REAL(CMISSDP), ALLOCATABLE :: StartingXiGet(:)

  REAL(CMISSDP), DIMENSION(2,3) :: DataPointValues!(number_of_data_points,dimension)
  REAL(CMISSDP), DIMENSION(8) :: ElementNumbers!(number_of_data_points,dimension)
  INTEGER(CMISSIntg) :: np

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
  NUMBER_GLOBAL_Z_ELEMENTS=1
  INTERPOLATION_TYPE=1
  
  NUMBER_GLOBAL_X_ELEMENTS_2=2
  NUMBER_GLOBAL_Y_ELEMENTS_2=2
  NUMBER_GLOBAL_Z_ELEMENTS_2=2
  INTERPOLATION_TYPE_2=1

  NumberOfComponents=FieldDependentNumberOfComponents
 !======================================
  !Intialise data points

!   DataPointValues(1,:)=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
!   DataPointValues(2,:)=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)  
!   DataPointValues(3,:)=(/1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP/)
!   DataPointValues(4,:)=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)  
!   DataPointValues(5,:)=(/1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP/)
!   DataPointValues(6,:)=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)  
!   DataPointValues(7,:)=(/0.5_CMISSDP,0.5_CMISSDP,0.5_CMISSDP/)
!   DataPointValues(8,:)=(/0.5_CMISSDP,0.5_CMISSDP,0.0_CMISSDP/)  
 DataPointValues(1,:)= (/0.5_CMISSDP,0.5_CMISSDP,0.5_CMISSDP/)
 DataPointValues(2,:)= (/0.8_CMISSDP,0.6_CMISSDP,0.7_CMISSDP/)
!================================================================================  
  

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

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

!!!!!!!!!!!!!!!!!
! First mesh
!!!!!!!!!!!!!!!!!
  
  
  !Start the creation of the region
  CALL CMISSRegionTypeInitialise(Region1,Err)
  CALL CMISSRegionCreateStart(RegionOneUserNumber,WorldRegion,Region1,Err)
  CALL CMISSRegionLabelSet(Region1,"FirstMeshRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region1,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region1,Err)
  
  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasisTypeInitialise(Basis1,Err)
  CALL CMISSBasisCreateStart(BasisOneUserNumber,Basis1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSBasisNumberOfXiSet(Basis1,2,Err)
    CALL CMISSBasisInterpolationXiSet(Basis1,(/1,1/),Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis1,(/2,2/),Err) 
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL CMISSBasisNumberOfXiSet(Basis1,3,Err)
    CALL CMISSBasisInterpolationXiSet(Basis1,(/1,1,1/),Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis1,(/2,2,2/),Err) 
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(Basis1,Err)

  !Start the creation of a generated mesh in the region - the parent mesh
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshOneUserNumber,Region1,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis1,Err)   
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
  CALL CMISSMeshTypeInitialise(FirstMesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshOneUserNumber,FirstMesh,Err)

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition1,Err)
  CALL CMISSDecompositionCreateStart(DecompositionOneUserNumber,FirstMesh,Decomposition1,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition1,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition1,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition1,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField1,Err)
  CALL CMISSFieldCreateStart(GeometricFieldOneUserNumber,Region1,GeometricField1,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField1,Decomposition1,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField1,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField1,CMISSFieldUVariableType,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(GeometricField1,CMISSFieldUVariableType,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField1,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField1,GeneratedMesh,Err)

  !Create the dependent field with 1 variables and 1 component
  CALL CMISSFieldTypeInitialise(DependentField1,Err)
  CALL CMISSFieldCreateStart(DependentFieldOneUserNumber,Region1,DependentField1,Err)
  CALL CMISSFieldTypeSet(DependentField1,CMISSFieldGeneralType,Err)
  CALL CMISSFieldMeshDecompositionSet(DependentField1,Decomposition1,Err)
  CALL CMISSFieldGeometricFieldSet(DependentField1,GeometricField1,Err)
  CALL CMISSFieldDependentTypeSet(DependentField1,CMISSFieldDependentType,Err)
  CALL CMISSFieldNumberOfVariablesSet(DependentField1,FieldDependentNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField1,CMISSFieldUVariableType,FieldDependentNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField1,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField1,CMISSFieldUVariableType,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(DependentField1,CMISSFieldUVariableType,3,1,Err)
  ENDIF
  !CALL CMISSFieldScalingTypeSet(DependentField,CMISSFieldUnitScaling,Err)
  CALL CMISSFieldCreateFinish(DependentField1,Err)

  !Initialise the field with an initial guess
  CALL CMISSFieldComponentValuesInitialise(DependentField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,0.5_CMISSDP,Err)

  !Field value with a function x^2 + y^2 + z^2
  DO node_idx=1,(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  CALL CMISSFieldParameterSetGetNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,node_idx,1,x,Err)
  CALL CMISSFieldParameterSetGetNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,node_idx,2,y,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
  CALL CMISSFieldParameterSetGetNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,node_idx,3,z,Err)
  ELSE
  Z = 0
  ENDIF
  FieldValue = x**2+y**2+z**2
  CALL CMISSFieldParameterSetUpdateNode(DependentField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,node_idx,1,fieldvalue,Err)

  WRITE(*,*) 'The values at the nodes for field 1',fieldvalue
  ENDDO

   CALL CMISSNodesNumberOfNodesGet(RegionOneUserNumber,NumberOfNodes1,Err)
   WRITe(*,*) 'nn',NumberOfNodes1
  
   DO node_idx=1,NumberofNodes1
   CALL CMISSFieldParameterSetGetNode(DependentField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,node_idx,1,value,Err)
   WRITE(*,*) 'The values at the nodes for field 1 before',value
   ENDDO

!!!!!!!!!!!!!!!!!!
! The second mesh
!!!!!!!!!!!!!!!!!!

  !Start the creation of the region
  CALL CMISSRegionTypeInitialise(Region2,Err)
  CALL CMISSRegionCreateStart(RegionTwoUserNumber,WorldRegion,Region2,Err)
  CALL CMISSRegionLabelSet(Region2,"SecondMeshRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region2,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region2,Err)
  
  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasisTypeInitialise(Basis2,Err)
  CALL CMISSBasisCreateStart(BasisTwoUserNumber,Basis2,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS_2==0) THEN
    !Set the basis to be a bilinear Lagrange basis
    CALL CMISSBasisNumberOfXiSet(Basis2,2,Err)
    CALL CMISSBasisInterpolationXiSet(Basis2,(/1,1/),Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis2,(/3,3/),Err) 
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL CMISSBasisNumberOfXiSet(Basis2,3,Err)
        CALL CMISSBasisInterpolationXiSet(Basis2,(/1,1,1/),Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis2,(/3,3,3/),Err) 
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(Basis2,Err)

  !Start the creation of a generated mesh in the region - the child mesh
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshTwoUserNumber,Region2,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis2,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS_2==0) THEN
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS_2,NUMBER_GLOBAL_Y_ELEMENTS_2],Err)
  ELSE
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS_2,NUMBER_GLOBAL_Y_ELEMENTS_2, &
      & NUMBER_GLOBAL_Z_ELEMENTS_2],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(SecondMesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshTwoUserNumber,SecondMesh,Err)

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition2,Err)
  CALL CMISSDecompositionCreateStart(DecompositionTwoUserNumber,SecondMesh,Decomposition2,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition2,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition2,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition2,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField2,Err)
  CALL CMISSFieldCreateStart(GeometricFieldTwoUserNumber,Region2,GeometricField2,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField2,Decomposition2,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField2,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField2,CMISSFieldUVariableType,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS_2/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(GeometricField2,CMISSFieldUVariableType,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField2,Err)
  
  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField2,GeneratedMesh,Err)

 !Create the dependent field with 1 variables and 1 component
  CALL CMISSFieldTypeInitialise(DependentField2,Err)
  CALL CMISSFieldCreateStart(DependentFieldTwoUserNumber,Region2,DependentField2,Err)
  CALL CMISSFieldTypeSet(DependentField2,CMISSFieldGeneralType,Err)
  CALL CMISSFieldMeshDecompositionSet(DependentField2,Decomposition2,Err)
  CALL CMISSFieldGeometricFieldSet(DependentField2,GeometricField2,Err)
  CALL CMISSFieldDependentTypeSet(DependentField2,CMISSFieldDependentType,Err)
  CALL CMISSFieldNumberOfVariablesSet(DependentField2,FieldDependentNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField2,CMISSFieldUVariableType,FieldDependentNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField2,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField2,CMISSFieldUVariableType,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS_2/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(DependentField2,CMISSFieldUVariableType,3,1,Err)
  ENDIF
   !CALL CMISSFieldScalingTypeSet(DependentField2,CMISSFieldUnitScaling,Err)
  CALL CMISSFieldCreateFinish(DependentField2,Err)

  !Initialise the field with an initial guess
  CALL CMISSFieldComponentValuesInitialise(DependentField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,0.1_CMISSDP,Err)

  
   CALL CMISSNodesNumberOfNodesGet(RegionTwoUserNumber,NumberOfNodes2,Err)
   !Before
   DO node_idx=1,NumberOfNodes2
   CALL CMISSFieldParameterSetGetNode(DependentField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,node_idx,1,value,Err)
   WRITE(*,*) 'The values at the nodes for field 2 before',value
   ENDDO

  !Start creating an embedded mesh
  CALL CMISSEmbeddedMeshTypeInitialise(MeshEmbedding,Err)
  
  CALL CMISSMeshEmbeddingCreate(MeshEmbedding, FirstMesh, SecondMesh, Err)
  
   Allocate(ParentXiCoords(NumberofComponents,NumberOfNodes1)) 
   DO node_idx = 1,NumberofNodes1
   DO dim_idx = 1,NumberOfComponents
   CALL CMISSFieldParameterSetGetNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1, &
     & node_idx,dim_idx,value,Err)
   ParentXiCoords(dim_idx,node_idx) = value
   ENDDO
   ENDDO

  ! Get the positions of the child nodes
  Allocate(ChildXiCoords(NumberofComponents,NumberOfNodes2)) 
    DO node_idx = 1,NumberofNodes2
    DO dim_idx = 1,NumberOfComponents
    CALL CMISSFieldParameterSetGetNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1, &
      & node_idx,dim_idx,value,Err)
     ChildXiCoords(dim_idx,node_idx) = value
    ENDDO
    ENDDO

    
    Allocate(NodeNumbers(NumberofNodes1))
    Allocate(NodeNumbers1(NumberofNodes1))
    Allocate(NodeNumbers2(NumberofNodes2))

   CALL CMISSMeshNumberOfElementsGet(RegionOneUserNumber,MeshOneUserNumber,NumberOfElements,Err)
   CALL CMISSMeshNumberOfElementsGet(RegionTwoUserNumber,MeshTwoUserNumber,NumberOfElements2,Err)
   !Preprocessing step:

!   By each element 
    counts = 1
    Allocate(ElemArray(NumberOfElements2))
    Allocate(NodeElem(NumberOfElements,NumberofNodes2))
    DO elem_idx = 1,NumberOfElements
     CALL CMISSMeshElementsNodesGet(RegionOneUserNumber,MeshOneUserNumber,1,elem_idx, &
      & NodeNumbers1,Err)
      DO elem_idx2 =1,NumberOfElements2
      CALL CMISSMeshElementsNodesGet(RegionTwoUserNumber,MeshTwoUserNumber,1,elem_idx2, &
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
    CALL CMISSMeshEmbeddingSetChildNodePosition(MeshEmbedding,elem_idx, NodeNumbers2,ChildXiCoords, Err)
   ENDDO


   !Data field from one mesh to the other
    CALL CMISSMeshEmbeddingPushData(MeshEmbedding, DependentField1,1,DependentField2, 1, Err)
 
   !After
   DO node_idx=1,NumberofNodes2
   CALL CMISSFieldParameterSetGetNode(DependentField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,node_idx,1,value,Err)
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
      CALL CMISSMeshEmbeddingSetGaussPointData(MeshEmbedding,elem_idx,gauss_idx, &
        & C1,elem_idx2,C2, Err)
      ENDDO
   ENDDO
   
   !Move field data

    !CALL CMISSMeshEmbeddingPullGaussPointData(MeshEmbedding, DependentField1, 1,DependentField2, 1, Err)


   !After
   !DO gauss_idx=1,1
   !CALL CMISSFieldParameterSetGetGaussPoint(DependentField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,gauss_idx,1,value,Err)
   !WRITE(*,*) 'The values at the nodes for field 1 after',value
   !ENDDO

  !================================================Gauss point map==================================
    DO dim_idx = 1,3
       CALL CMISSFieldParameterSetGetGaussPointCoord(MeshEmbedding,dim_idx,NGP,C1,Err)        
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
!   CALL CMISSDataPointsCreateStart(RegionTwoUserNumber,SIZE(DataPointValues,1),Err)
!   DO np=1,SIZE(DataPointValues,1)
!     CALL CMISSDataPointsValuesSet(RegionTwoUserNumber,np,DataPointValues(np,:),Err)     
!   ENDDO
!   CALL CMISSDataPointsCreateFinish(RegionTwoUserNumber,Err)  
!   !=========================================================================================================================
!   !=========================================================================================================================
!   !Create a data projection
!   CALL CMISSDataProjectionCreateStart(RegionTwoUserNumber,DependentFieldTwoUserNumber,RegionTwoUserNumber,Err)
!   !=========================================================================================================================
!   !Test parameter set functions
!   CALL CMISSDataProjectionAbsoluteToleranceSet(RegionTwoUserNumber,AbsoluteToleranceSet,Err) !test
!   CALL CMISSDataProjectionMaximumIterationUpdateSet(RegionTwoUserNumber,MaximumIterationUpdateSet,Err) !test
!   CALL CMISSDataProjectionMaximumNumberOfIterationsSet(RegionTwoUserNumber,MaximumNumberOfIterationsSet,Err) !test
!   CALL CMISSDataProjectionNumberOfClosestElementsSet(RegionTwoUserNumber,NumberOfClosestElementsSet,Err) !test
!   CALL CMISSDataProjectionProjectionTypeSet(RegionTwoUserNumber,ProjectionTypeSet,Err)
!   CALL CMISSDataProjectionRelativeToleranceSet(RegionTwoUserNumber,RelativeToleranceSet,Err) !test
!   CALL CMISSDataProjectionStartingXiSet(RegionTwoUserNumber,StartingXiSet,Err) !test
!   !=========================================================================================================================
!   !Finish data projection  
!   CALL CMISSDataProjectionCreateFinish(RegionTwoUserNumber,Err)
!   !=========================================================================================================================
!   !Test parameter get functions
!   CALL CMISSDataProjectionAbsoluteToleranceGet(RegionTwoUserNumber,AbsoluteToleranceGet,Err) !test
!   CALL CMISSDataProjectionMaximumIterationUpdateGet(RegionTwoUserNumber,MaximumIterationUpdateGet,Err) !test
!   CALL CMISSDataProjectionMaximumNumberOfIterationsGet(RegionTwoUserNumber,MaximumNumberOfIterationsGet,Err) !test
!   CALL CMISSDataProjectionNumberOfClosestElementsGet(RegionTwoUserNumber,NumberOfClosestElementsGet,Err) !test
!   CALL CMISSDataProjectionProjectionTypeGet(RegionTwoUserNumber,ProjectionTypeGet,Err) !test
!   CALL CMISSDataProjectionRelativeToleranceGet(RegionTwoUserNumber,RelativeToleranceGet,Err) !test
!   CALL CMISSDataProjectionStartingXiGet(RegionTwoUserNumber,StartingXiGet,Err) !test !! This is the GP xi  
!   
!   !=========================================================================================================================
!   !Start data projection
!   CALL CMISSDataProjectionEvaluate(RegionTwoUserNumber,Err)
!   
!   !Get the closest elementnumber back - in an array
!   !=========================================================================================================================
!   !Destroy used types
!   CALL CMISSDataProjectionDestroy(RegionTwoUserNumber,Err)
!   CALL CMISSDataPointsDestroy(RegiontwoUserNumber,Err)
!     
!   CALL CMISSRegionDestroy(RegionTwoUserNumber,Err)
  !=========================================================================================================================
  !Create Data Points and set the values
  CALL CMISSDataPointsCreateStart(RegionOneUserNumber,SIZE(DataPointValues,1),Err)
  DO np=1,SIZE(DataPointValues,1)
    CALL CMISSDataPointsValuesSet(RegionOneUserNumber,np,DataPointValues(np,:),Err)     
  ENDDO
  CALL CMISSDataPointsCreateFinish(RegionOneUserNumber,Err)  
  !=========================================================================================================================
  !=========================================================================================================================
  !Create a data projection
  CALL CMISSDataProjectionCreateStart(RegionOneUserNumber,DependentFieldOneUserNumber,RegionOneUserNumber,Err)
  !=========================================================================================================================
  !Test parameter set functions
  CALL CMISSDataProjectionAbsoluteToleranceSet(RegionOneUserNumber,AbsoluteToleranceSet,Err) !test
  CALL CMISSDataProjectionMaximumIterationUpdateSet(RegionOneUserNumber,MaximumIterationUpdateSet,Err) !test
  CALL CMISSDataProjectionMaximumNumberOfIterationsSet(RegionOneUserNumber,MaximumNumberOfIterationsSet,Err) !test
  CALL CMISSDataProjectionNumberOfClosestElementsSet(RegionOneUserNumber,NumberOfClosestElementsSet,Err) !test
  CALL CMISSDataProjectionProjectionTypeSet(RegionOneUserNumber,ProjectionTypeSet,Err)
  CALL CMISSDataProjectionRelativeToleranceSet(RegionOneUserNumber,RelativeToleranceSet,Err) !test
  CALL CMISSDataProjectionStartingXiSet(RegionOneUserNumber,StartingXiSet,Err) !test
  !=========================================================================================================================
  !Finish data projection  
  CALL CMISSDataProjectionCreateFinish(RegionOneUserNumber,Err)
  !=========================================================================================================================
  !Test parameter get functions
  CALL CMISSDataProjectionAbsoluteToleranceGet(RegionOneUserNumber,AbsoluteToleranceGet,Err) !test
  CALL CMISSDataProjectionMaximumIterationUpdateGet(RegionOneUserNumber,MaximumIterationUpdateGet,Err) !test
  CALL CMISSDataProjectionMaximumNumberOfIterationsGet(RegionOneUserNumber,MaximumNumberOfIterationsGet,Err) !test
  CALL CMISSDataProjectionNumberOfClosestElementsGet(RegionOneUserNumber,NumberOfClosestElementsGet,Err) !test
  CALL CMISSDataProjectionProjectionTypeGet(RegionOneUserNumber,ProjectionTypeGet,Err) !test
  CALL CMISSDataProjectionRelativeToleranceGet(RegionOneUserNumber,RelativeToleranceGet,Err) !test
  CALL CMISSDataProjectionStartingXiGet(RegionOneUserNumber,StartingXiGet,Err) !test !! This is the GP xi  
  
  !=========================================================================================================================
  !Start data projection
  CALL CMISSDataProjectionEvaluate(RegionOneUserNumber,Err)
  
  !Get the closest elementnumber back - in an array
  !=========================================================================================================================
  !Destroy used types
  CALL CMISSDataProjectionDestroy(RegionOneUserNumber,Err)
  CALL CMISSDataPointsDestroy(RegionOneUserNumber,Err)
    
  CALL CMISSRegionDestroy(RegionOneUserNumber,Err)
  CALL CMISSCoordinateSystemDestroy(CoordinateSystemUserNumber,Err)  
  
  !=========================================================================================================================

  !=================================================================================================
  !END SUBROUTINE gauss_point_mapping

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM EMBEDDEDMESHEXAMPLE
