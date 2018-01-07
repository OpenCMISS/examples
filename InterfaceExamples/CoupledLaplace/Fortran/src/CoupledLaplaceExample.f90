!> \file
!> \author Chris Bradley
!> \brief This is an example program which solves a weakly coupled Laplace equation in two regions using OpenCMISS calls.
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

!> \example InterfaceExamples/CoupledLaplace/src/CoupledLaplaceExample.f90
!! Example program which sets up a field in two regions using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/InterfaceExamples/CoupledLaplace/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/InterfaceExamples/CoupledLaplace/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM COUPLEDLAPLACE

  USE OpenCMISS
  USE OpenCMISS_Iron
  
#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=3.0_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystem1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystem2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemInterfaceUserNumber=42
  INTEGER(CMISSIntg), PARAMETER :: Region1UserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: Region2UserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: Basis1UserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: Basis2UserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: InterfaceBasisUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMesh1UserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMesh2UserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: InterfaceGeneratedMeshUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: Mesh1UserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: Mesh2UserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: InterfaceMeshUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: Decomposition1UserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: Decomposition2UserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: InterfaceDecompositionUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: GeometricField1UserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: GeometricField2UserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: InterfaceGeometricFieldUserNumber=19
  INTEGER(CMISSIntg), PARAMETER :: EquationsSet1UserNumber=20
  INTEGER(CMISSIntg), PARAMETER :: EquationsSet2UserNumber=21
  INTEGER(CMISSIntg), PARAMETER :: DependentField1UserNumber=22
  INTEGER(CMISSIntg), PARAMETER :: DependentField2UserNumber=23
  INTEGER(CMISSIntg), PARAMETER :: InterfaceUserNumber=24
  INTEGER(CMISSIntg), PARAMETER :: InterfaceConditionUserNumber=25
  INTEGER(CMISSIntg), PARAMETER :: LagrangeFieldUserNumber=26
  INTEGER(CMISSIntg), PARAMETER :: CoupledProblemUserNumber=27
  INTEGER(CMISSIntg), PARAMETER :: InterfaceMappingBasisUserNumber=28
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetField1UserNumber=40
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetField2UserNumber=41
 
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
    & INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI,NUMBER_OF_NODE_XI
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT

  INTEGER(CMISSIntg) :: EquationsSet1Index,EquationsSet2Index
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain
  INTEGER(CMISSIntg) :: InterfaceConditionIndex
  INTEGER(CMISSIntg) :: Mesh1Index,Mesh2Index
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: y_element_idx,z_element_idx,mesh_local_y_node,mesh_local_z_node
  REAL(CMISSRP) :: XI2(2),XI3(3)

  !CMISS variables

  TYPE(cmfe_BasisType) :: Basis1,Basis2,InterfaceBasis,InterfaceMappingBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem1,CoordinateSystem2,CoordinateSystemInterface,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition1,Decomposition2,InterfaceDecomposition
  TYPE(cmfe_EquationsType) :: Equations1,Equations2
  TYPE(cmfe_EquationsSetType) :: EquationsSet1,EquationsSet2
  TYPE(cmfe_FieldType) :: GeometricField1,GeometricField2,InterfaceGeometricField,DependentField1, &
    & DependentField2,LagrangeField,EquationsSetField1,EquationsSetField2
  TYPE(cmfe_FieldsType) :: Fields1,Fields2,InterfaceFields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh1,GeneratedMesh2,InterfaceGeneratedMesh
  TYPE(cmfe_InterfaceType) :: Interface
  TYPE(cmfe_InterfaceConditionType) :: InterfaceCondition
  TYPE(cmfe_InterfaceEquationsType) :: InterfaceEquations
  TYPE(cmfe_InterfaceMeshConnectivityType) :: InterfaceMeshConnectivity
  TYPE(cmfe_MeshType) :: Mesh1,Mesh2,InterfaceMesh
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_ProblemType) :: CoupledProblem
  TYPE(cmfe_RegionType) :: Region1,Region2,WorldRegion
  TYPE(cmfe_SolverType) :: CoupledSolver
  TYPE(cmfe_SolverEquationsType) :: CoupledSolverEquations
  
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
    IF(INTERPOLATION_TYPE<=0.OR.INTERPOLATION_TYPE>4) CALL HANDLE_ERROR("Invalid Interpolation specification.")
  ELSE IF(NUMBER_OF_ARGUMENTS == 0) THEN
    NUMBER_GLOBAL_X_ELEMENTS=2
    NUMBER_GLOBAL_Y_ELEMENTS=2
    NUMBER_GLOBAL_Z_ELEMENTS=0
    INTERPOLATION_TYPE=1
  ELSE
    CALL HANDLE_ERROR("Invalid number of arguments.")
  ENDIF

  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  !Set error handling mode
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
 
  !Set diganostics for testing
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["SOLVER_MAPPING_CALCULATE         ", &
  !  & "SOLVER_MATRIX_STRUCTURE_CALCULATE"],Err)
  
  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  
  !Start the creation of a new RC coordinate system for the first region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM(1) << == '
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem1,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystem1UserNumber,CoordinateSystem1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem1,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem1,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem1,Err)

  !Start the creation of a new RC coordinate system for the second region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM(2) << == '
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem2,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystem2UserNumber,CoordinateSystem2,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem2,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem2,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem2,Err)
  
  !Start the creation of a new RC coordinate system for the second region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM(Interface) << == '
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystemInterface,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemInterfaceUserNumber,CoordinateSystemInterface,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystemInterface,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystemInterface,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystemInterface,Err)
  
  !Start the creation of the first region
  PRINT *, ' == >> CREATING REGION(1) << == '
  CALL cmfe_Region_Initialise(Region1,Err)
  CALL cmfe_Region_CreateStart(Region1UserNumber,WorldRegion,Region1,Err)
  CALL cmfe_Region_LabelSet(Region1,"Region1",Err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region1,CoordinateSystem1,Err)
  !Finish the creation of the first region
  CALL cmfe_Region_CreateFinish(Region1,Err)

  !Start the creation of the second region
  PRINT *, ' == >> CREATING REGION(2) << == '
  CALL cmfe_Region_Initialise(Region2,Err)
  CALL cmfe_Region_CreateStart(Region2UserNumber,WorldRegion,Region2,Err)
  CALL cmfe_Region_LabelSet(Region2,"Region2",Err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region2,CoordinateSystem2,Err)
  !Finish the creation of the second region
  CALL cmfe_Region_CreateFinish(Region2,Err)

  !Start the creation of a bI/tri-linear-Lagrange basis
  PRINT *, ' == >> CREATING BASIS(1) << == '
  CALL cmfe_Basis_Initialise(Basis1,Err)
  CALL cmfe_Basis_CreateStart(Basis1UserNumber,Basis1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis1,2,Err)
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis1,3,Err)
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis1,Err)
   
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1)
    NUMBER_OF_GAUSS_XI=2
  CASE(2)
    NUMBER_OF_GAUSS_XI=3
  CASE(3,4)
    NUMBER_OF_GAUSS_XI=4
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  !Start the creation of a bI/tri-XXX-Lagrange basis
  PRINT *, ' == >> CREATING BASIS(2) << == '
  CALL cmfe_Basis_Initialise(Basis2,Err)
  CALL cmfe_Basis_CreateStart(Basis2UserNumber,Basis2,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bi-XXX Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis2,2,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis2,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis2,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
  ELSE
    !Set the basis to be a tri-XXX Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis2,3,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis2,[INTERPOLATION_TYPE,INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis2,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI, &
      & NUMBER_OF_GAUSS_XI],Err)
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis2,Err)
  
  !Start the creation of a generated mesh in the first region
  PRINT *, ' == >> CREATING GENERATED MESH(1) << == '
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh1,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMesh1UserNumber,Region1,GeneratedMesh1,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh1,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh1,Basis1,Err)   
  !Define the mesh on the first region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh1,[WIDTH,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh1,[NUMBER_GLOBAL_X_ELEMENTS, &
      & NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh1,[WIDTH,HEIGHT,LENGTH],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh1,[NUMBER_GLOBAL_X_ELEMENTS, &
      & NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the first region
  CALL cmfe_Mesh_Initialise(Mesh1,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh1,Mesh1UserNumber,Mesh1,Err)

  !Start the creation of a generated mesh in the second region
  PRINT *, ' == >> CREATING GENERATED MESH(2) << == '
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh2,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMesh2UserNumber,Region2,GeneratedMesh2,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh2,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh2,Basis2,Err)   
  !Define the mesh on the second region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL cmfe_GeneratedMesh_OriginSet(GeneratedMesh2,[WIDTH,0.0_CMISSRP],Err)
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh2,[WIDTH,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh2,[NUMBER_GLOBAL_X_ELEMENTS, &
      & NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL cmfe_GeneratedMesh_OriginSet(GeneratedMesh2,[WIDTH,0.0_CMISSRP,0.0_CMISSRP],Err)
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh2,[WIDTH,HEIGHT,LENGTH],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh2,[NUMBER_GLOBAL_X_ELEMENTS, &
      & NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the second region
  CALL cmfe_Mesh_Initialise(Mesh2,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh2,Mesh2UserNumber,Mesh2,Err)

  !Create an interface between the two meshes
  PRINT *, ' == >> CREATING INTERFACE << == '
  CALL cmfe_Interface_Initialise(Interface,Err)
  CALL cmfe_Interface_CreateStart(InterfaceUserNumber,WorldRegion,Interface,Err)
  CALL cmfe_Interface_LabelSet(Interface,"Interface",Err)
  !Add in the two meshes
  CALL cmfe_Interface_MeshAdd(Interface,Mesh1,Mesh1Index,Err)
  CALL cmfe_Interface_MeshAdd(Interface,Mesh2,Mesh2Index,Err)
  CALL cmfe_Interface_CoordinateSystemSet(Interface,CoordinateSystemInterface,Err)
  !Finish creating the interface
  CALL cmfe_Interface_CreateFinish(Interface,Err)

  !Start the creation of a (bi)-linear-Lagrange basis
  PRINT *, ' == >> CREATING INTERFACE BASIS << == '
  CALL cmfe_Basis_Initialise(InterfaceBasis,Err)
  CALL cmfe_Basis_CreateStart(InterfaceBasisUserNumber,InterfaceBasis,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN  
    !Set the basis to be a linear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(InterfaceBasis,1,Err)
    CALL cmfe_Basis_InterpolationXiSet(InterfaceBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  ELSE
    !Set the basis to be a bilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(InterfaceBasis,2,Err)
    CALL cmfe_Basis_InterpolationXiSet(InterfaceBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
      & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(InterfaceBasis,Err)

  !Start the creation of a (bi)-linear-Lagrange basis
  PRINT *, ' == >> CREATING INTERFACE MAPPING BASIS << == '
  CALL cmfe_Basis_Initialise(InterfaceMappingBasis,Err)
  CALL cmfe_Basis_CreateStart(InterfaceMappingBasisUserNumber,InterfaceMappingBasis,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN  
    !Set the basis to be a linear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(InterfaceMappingBasis,1,Err)
    CALL cmfe_Basis_InterpolationXiSet(InterfaceMappingBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  ELSE
    !Set the basis to be a bilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(InterfaceMappingBasis,2,Err)
    CALL cmfe_Basis_InterpolationXiSet(InterfaceMappingBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
      & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(InterfaceMappingBasis,Err)
  
  !Start the creation of a generated mesh for the interface
  PRINT *, ' == >> CREATING INTERFACE GENERATED MESH << == '
  CALL cmfe_GeneratedMesh_Initialise(InterfaceGeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(InterfaceGeneratedMeshUserNumber,Interface,InterfaceGeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(InterfaceGeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(InterfaceGeneratedMesh,InterfaceBasis,Err)   
  !Define the mesh on the interface
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL cmfe_GeneratedMesh_OriginSet(InterfaceGeneratedMesh,[WIDTH,0.0_CMISSRP],Err)
    CALL cmfe_GeneratedMesh_ExtentSet(InterfaceGeneratedMesh,[0.0_CMISSRP,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(InterfaceGeneratedMesh,[NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL cmfe_GeneratedMesh_OriginSet(InterfaceGeneratedMesh,[WIDTH,0.0_CMISSRP,0.0_CMISSRP],Err)
    CALL cmfe_GeneratedMesh_ExtentSet(InterfaceGeneratedMesh,[0.0_CMISSRP,HEIGHT,LENGTH],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(InterfaceGeneratedMesh,[NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF    
  !Finish the creation of a generated mesh in interface
  CALL cmfe_Mesh_Initialise(InterfaceMesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(InterfaceGeneratedMesh,InterfaceMeshUserNumber,InterfaceMesh,Err)

  !Couple the interface meshes
  PRINT *, ' == >> CREATING INTERFACE MESHES CONNECTIVITY << == '
  CALL cmfe_InterfaceMeshConnectivity_Initialise(InterfaceMeshConnectivity,Err)
  CALL cmfe_InterfaceMeshConnectivity_CreateStart(Interface,InterfaceMesh,InterfaceMeshConnectivity,Err)
  CALL cmfe_InterfaceMeshConnectivity_BasisSet(InterfaceMeshConnectivity,InterfaceMappingBasis,Err)
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1,4)
    NUMBER_OF_NODE_XI=2
  CASE(2)
    NUMBER_OF_NODE_XI=3
  CASE(3)
    NUMBER_OF_NODE_XI=4
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN   
    DO y_element_idx=1,NUMBER_GLOBAL_Y_ELEMENTS
      !Map the interface element to the elements in mesh 1
      CALL cmfe_InterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity,y_element_idx,Mesh1Index, &
        y_element_idx*NUMBER_GLOBAL_X_ELEMENTS,Err)
      XI2 = [ 1.0_CMISSRP, 0.0_CMISSRP ]
      CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx,Mesh1Index, &
        & y_element_idx*NUMBER_GLOBAL_X_ELEMENTS,1,1,XI2,Err)
      XI2 = [ 1.0_CMISSRP, 1.0_CMISSRP ]
      CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx,Mesh1Index, &
        & y_element_idx*NUMBER_GLOBAL_X_ELEMENTS,2,1,XI2,Err)      
      !Map the interface element to the elements in mesh 2
      CALL cmfe_InterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity,y_element_idx,Mesh2Index, &
        & 1+(y_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS,Err)
      DO mesh_local_y_node = 1,NUMBER_OF_NODE_XI-1
        XI2 = [ 0.0_CMISSRP, REAL(mesh_local_y_node-1,CMISSRP)/REAL(NUMBER_OF_NODE_XI-1,CMISSRP) ]
        CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx,Mesh2Index, &
          & 1+(y_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS,1,1,XI2,Err)
        XI2 = [ 0.0_CMISSRP, REAL(mesh_local_y_node,CMISSRP)/REAL(NUMBER_OF_NODE_XI-1,CMISSRP) ]
        CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx,Mesh2Index, &
          & 1+(y_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS,2,1,XI2,Err)
      ENDDO !mesh_local_y_node
    ENDDO !y_element_idx
  ELSE
    DO y_element_idx=1,NUMBER_GLOBAL_Y_ELEMENTS
      DO z_element_idx=1,NUMBER_GLOBAL_Z_ELEMENTS
        !Map the interface element to the elements in mesh 1
        CALL cmfe_InterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity, &
          & y_element_idx+(z_element_idx-1)*NUMBER_GLOBAL_Y_ELEMENTS,Mesh1Index, &
          y_element_idx*NUMBER_GLOBAL_X_ELEMENTS+(z_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS* &
          & NUMBER_GLOBAL_Y_ELEMENTS,Err)
        XI3 = [ 1.0_CMISSRP, 0.0_CMISSRP, 0.0_CMISSRP ]
        CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
          & (z_element_idx-1)*NUMBER_GLOBAL_Y_ELEMENTS,Mesh1Index,y_element_idx* &
          & NUMBER_GLOBAL_X_ELEMENTS+(z_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS* &
          & NUMBER_GLOBAL_Y_ELEMENTS,1,1,XI3,Err)
        XI3 = [ 1.0_CMISSRP, 1.0_CMISSRP, 0.0_CMISSRP ]
        CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
          & (z_element_idx-1)*NUMBER_GLOBAL_Y_ELEMENTS,Mesh1Index,y_element_idx* &
          & NUMBER_GLOBAL_X_ELEMENTS+(z_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS* &
          & NUMBER_GLOBAL_Y_ELEMENTS,2,1,XI3,Err)
        XI3 = [ 1.0_CMISSRP, 0.0_CMISSRP, 1.0_CMISSRP ]
        CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
          & (z_element_idx-1)*NUMBER_GLOBAL_Y_ELEMENTS,Mesh1Index,y_element_idx* &
          & NUMBER_GLOBAL_X_ELEMENTS+(z_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS* &
          & NUMBER_GLOBAL_Y_ELEMENTS,3,1,XI3,Err)
        XI3 = [ 1.0_CMISSRP, 1.0_CMISSRP, 1.0_CMISSRP ]
        CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
          & (z_element_idx-1)*NUMBER_GLOBAL_Y_ELEMENTS,Mesh1Index,y_element_idx* &
          & NUMBER_GLOBAL_X_ELEMENTS+(z_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS* &
          & NUMBER_GLOBAL_Y_ELEMENTS,4,1,XI3,Err)
        !Map the interface element to the elements in mesh 2
        CALL cmfe_InterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity,y_element_idx+ &
          & (z_element_idx-1)*NUMBER_GLOBAL_Y_ELEMENTS,Mesh2Index,1+(y_element_idx-1)* &
          & NUMBER_GLOBAL_X_ELEMENTS+(z_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS* &
          & NUMBER_GLOBAL_Y_ELEMENTS,Err)
        DO mesh_local_y_node = 1,NUMBER_OF_NODE_XI-1
          DO mesh_local_z_node = 1,NUMBER_OF_NODE_XI-1
            XI3 = [ 0.0_CMISSRP,  &
              & REAL(mesh_local_y_node-1,CMISSRP)/REAL(NUMBER_OF_NODE_XI-1,CMISSRP), &
              & REAL(mesh_local_z_node-1,CMISSRP)/REAL(NUMBER_OF_NODE_XI-1,CMISSRP) ]
            CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
              & (z_element_idx-1)*NUMBER_GLOBAL_Y_ELEMENTS,Mesh2Index,1+(y_element_idx-1)* &
              & NUMBER_GLOBAL_X_ELEMENTS+(z_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS* &
              & NUMBER_GLOBAL_Y_ELEMENTS,1,1,XI3,Err)
            XI3 = [ 0.0_CMISSRP,  &
              & REAL(mesh_local_y_node,CMISSRP)/REAL(NUMBER_OF_NODE_XI-1,CMISSRP), &
              & REAL(mesh_local_z_node-1,CMISSRP)/REAL(NUMBER_OF_NODE_XI-1,CMISSRP) ]
            CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
              & (z_element_idx-1)*NUMBER_GLOBAL_Y_ELEMENTS,Mesh2Index,1+(y_element_idx-1)* &
              & NUMBER_GLOBAL_X_ELEMENTS+(z_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS* &
              & NUMBER_GLOBAL_Y_ELEMENTS,2,1,XI3,Err)
            XI3 = [ 0.0_CMISSRP,  &
              & REAL(mesh_local_y_node-1,CMISSRP)/REAL(NUMBER_OF_NODE_XI-1,CMISSRP), &
              & REAL(mesh_local_z_node,CMISSRP)/REAL(NUMBER_OF_NODE_XI-1,CMISSRP) ]
            CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
              & (z_element_idx-1)*NUMBER_GLOBAL_Y_ELEMENTS,Mesh2Index,1+(y_element_idx-1)* &
              & NUMBER_GLOBAL_X_ELEMENTS+(z_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS* &
              & NUMBER_GLOBAL_Y_ELEMENTS,3,1,XI3,Err)
            XI3 = [ 0.0_CMISSRP,  &
              & REAL(mesh_local_y_node,CMISSRP)/REAL(NUMBER_OF_NODE_XI-1,CMISSRP), &
              & REAL(mesh_local_z_node,CMISSRP)/REAL(NUMBER_OF_NODE_XI-1,CMISSRP) ]
            CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,y_element_idx+ &
              & (z_element_idx-1)*NUMBER_GLOBAL_Y_ELEMENTS,Mesh2Index,1+(y_element_idx-1)* &
              & NUMBER_GLOBAL_X_ELEMENTS+(z_element_idx-1)*NUMBER_GLOBAL_X_ELEMENTS* &
              & NUMBER_GLOBAL_Y_ELEMENTS,4,1,XI3,Err)
          ENDDO !mesh_local_z_node
        ENDDO !mesh_local_y_node
      ENDDO !z_element_idx
    ENDDO !y_element_idx
  ENDIF
  CALL cmfe_InterfaceMeshConnectivity_CreateFinish(InterfaceMeshConnectivity,Err)

  !Create a decomposition for mesh1
  PRINT *, ' == >> CREATING MESH(1) DECOMPOSITION << == '
  CALL cmfe_Decomposition_Initialise(Decomposition1,Err)
  CALL cmfe_Decomposition_CreateStart(Decomposition1UserNumber,Mesh1,Decomposition1,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition1,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition1,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition1,Err)

  !Create a decomposition for mesh2
  PRINT *, ' == >> CREATING MESH(2) DECOMPOSITION << == '
  CALL cmfe_Decomposition_Initialise(Decomposition2,Err)
  CALL cmfe_Decomposition_CreateStart(Decomposition2UserNumber,Mesh2,Decomposition2,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition2,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition2,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition2,Err)
  
  !Create a decomposition for the interface mesh
  PRINT *, ' == >> CREATING INTERFACE DECOMPOSITION << == '
  CALL cmfe_Decomposition_Initialise(InterfaceDecomposition,Err)
  CALL cmfe_Decomposition_CreateStart(InterfaceDecompositionUserNumber,InterfaceMesh,InterfaceDecomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(InterfaceDecomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(InterfaceDecomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(InterfaceDecomposition,Err)

  !Start to create a default (geometric) field on the first region
  PRINT *, ' == >> CREATING MESH(1) GEOMETRIC FIELD << == '
  CALL cmfe_Field_Initialise(GeometricField1,Err)
  CALL cmfe_Field_CreateStart(GeometricField1UserNumber,Region1,GeometricField1,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField1,Decomposition1,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the first field
  CALL cmfe_Field_CreateFinish(GeometricField1,Err)

  !Start to create a default (geometric) field on the second region
  PRINT *, ' == >> CREATING MESH(2) GEOMETRIC FIELD << == '
  CALL cmfe_Field_Initialise(GeometricField2,Err)
  CALL cmfe_Field_CreateStart(GeometricField2UserNumber,Region2,GeometricField2,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField2,Decomposition2,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the second field
  CALL cmfe_Field_CreateFinish(GeometricField2,Err)

  !Update the geometric field parameters for the first field
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh1,GeometricField1,Err)
  !Update the geometric field parameters for the second field
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh2,GeometricField2,Err)

  !Create the equations set for the first region
  PRINT *, ' == >> CREATING EQUATION SET(1) << == '
  CALL cmfe_Field_Initialise(EquationsSetField1,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSet1,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSet1UserNumber,Region1,GeometricField1,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMFE_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE],EquationsSetField1UserNumber, &
    & EquationsSetField1,EquationsSet1,Err)
  !Set the equations set to be a standard Laplace problem
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet1,Err)

  !Create the equations set for the second region
  PRINT *, ' == >> CREATING EQUATION SET(2) << == '
  CALL cmfe_Field_Initialise(EquationsSetField2,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSet2,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSet2UserNumber,Region2,GeometricField2,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMFE_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE],EquationsSetField2UserNumber, &
    & EquationsSetField2,EquationsSet2,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet2,Err)

  !Create the equations set dependent field variables for the first equations set
  PRINT *, ' == >> CREATING DEPENDENT FIELD(1) << == '
  CALL cmfe_Field_Initialise(DependentField1,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet1,DependentField1UserNumber,DependentField1,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet1,Err)

  !Create the equations set dependent field variables for the second equations set
  PRINT *, ' == >> CREATING DEPENDENT FIELD(2) << == '
  CALL cmfe_Field_Initialise(DependentField2,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet2,DependentField2UserNumber,DependentField2,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet2,Err)

  !Create the equations set equations for the first equations set
  PRINT *, ' == >> CREATING EQUATIONS(1) << == '
  CALL cmfe_Equations_Initialise(Equations1,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet1,Equations1,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations1,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  !CALL cmfe_Equations_OutputTypeSet(Equations1,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations1,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations1,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations1,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet1,Err)

  !Create the equations set equations for the second equations set
  PRINT *, ' == >> CREATING EQUATIONS(2) << == '
  CALL cmfe_Equations_Initialise(Equations2,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet2,Equations2,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations2,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  !CALL cmfe_Equations_OutputTypeSet(Equations2,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations2,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations2,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations2,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet2,Err)

  !Start to create a default (geometric) field on the Interface
  PRINT *, ' == >> CREATING INTERFACE GEOMETRIC FIELD << == '
  CALL cmfe_Field_Initialise(InterfaceGeometricField,Err)
  CALL cmfe_Field_CreateStart(InterfaceGeometricFieldUserNumber,Interface,InterfaceGeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(InterfaceGeometricField,InterfaceDecomposition,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(InterfaceGeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(InterfaceGeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(InterfaceGeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the first field
  CALL cmfe_Field_CreateFinish(InterfaceGeometricField,Err)

  !Update the geometric field parameters for the interface field
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(InterfaceGeneratedMesh,InterfaceGeometricField,Err)

  !Create an interface condition between the two meshes
  PRINT *, ' == >> CREATING INTERFACE CONDITIONS << == '
  CALL cmfe_InterfaceCondition_Initialise(InterfaceCondition,Err)
  CALL cmfe_InterfaceCondition_CreateStart(InterfaceConditionUserNumber,Interface,InterfaceGeometricField, &
    & InterfaceCondition,Err)
  !Specify the method for the interface condition
  CALL cmfe_InterfaceCondition_MethodSet(InterfaceCondition,CMFE_INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,Err)
  !Specify the type of interface condition operator
  CALL cmfe_InterfaceCondition_OperatorSet(InterfaceCondition,CMFE_INTERFACE_CONDITION_FIELD_CONTINUITY_OPERATOR,Err)
  !Add in the dependent variables from the equations sets
  CALL cmfe_InterfaceCondition_DependentVariableAdd(InterfaceCondition,Mesh1Index,EquationsSet1, &
    & CMFE_FIELD_U_VARIABLE_TYPE,Err)
  CALL cmfe_InterfaceCondition_DependentVariableAdd(InterfaceCondition,Mesh2Index,EquationsSet2, &
    & CMFE_FIELD_U_VARIABLE_TYPE,Err)
  !Finish creating the interface condition
  CALL cmfe_InterfaceCondition_CreateFinish(InterfaceCondition,Err)

  !Create the Lagrange multipliers field
  PRINT *, ' == >> CREATING INTERFACE LAGRANGE FIELD << == '
  CALL cmfe_Field_Initialise(LagrangeField,Err)
  CALL cmfe_InterfaceCondition_LagrangeFieldCreateStart(InterfaceCondition,LagrangeFieldUserNumber,LagrangeField,Err)
  !Finish the Lagrange multipliers field
  CALL cmfe_InterfaceCondition_LagrangeFieldCreateFinish(InterfaceCondition,Err)

  !Create the interface condition equations
  PRINT *, ' == >> CREATING INTERFACE EQUATIONS << == '
  CALL cmfe_InterfaceEquations_Initialise(InterfaceEquations,Err)
  CALL cmfe_InterfaceCondition_EquationsCreateStart(InterfaceCondition,InterfaceEquations,Err)
  !Set the interface equations sparsity
  CALL cmfe_InterfaceEquations_SparsitySet(InterfaceEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the interface equations output
  CALL cmfe_InterfaceEquations_OutputTypeSet(InterfaceEquations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !Finish creating the interface equations
  CALL cmfe_InterfaceCondition_EquationsCreateFinish(InterfaceCondition,Err)
  
  !Start the creation of a coupled problem.
  PRINT *, ' == >> CREATING PROBLEM << == '
  CALL cmfe_Problem_Initialise(CoupledProblem,Err)
  CALL cmfe_Problem_CreateStart(CoupledProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS,CMFE_PROBLEM_LAPLACE_EQUATION_TYPE, &
    & CMFE_PROBLEM_STANDARD_LAPLACE_SUBTYPE],CoupledProblem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(CoupledProblem,Err)

  !Start the creation of the problem control loop for the coupled problem
  PRINT *, ' == >> CREATING PROBLEM CONTROL LOOP << == '
  CALL cmfe_Problem_ControlLoopCreateStart(CoupledProblem,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(CoupledProblem,Err)
 
  !Start the creation of the problem solver for the coupled problem
  PRINT *, ' == >> CREATING PROBLEM SOLVERS << == '
  CALL cmfe_Solver_Initialise(CoupledSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(CoupledProblem,Err)
  CALL cmfe_Problem_SolverGet(CoupledProblem,CMFE_CONTROL_LOOP_NODE,1,CoupledSolver,Err)
  !CALL cmfe_Solver_OutputTypeSet(CoupledSolver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(CoupledSolver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(CoupledSolver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(CoupledSolver,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  CALL cmfe_Solver_OutputTypeSet(CoupledSolver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  CALL cmfe_Solver_LinearTypeSet(CoupledSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL cmfe_Solver_LibraryTypeSet(CoupledSolver,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(CoupledProblem,Err)

  !Start the creation of the problem solver equations for the coupled problem
  PRINT *, ' == >> CREATING PROBLEM SOLVER EQUATIONS << == '
  CALL cmfe_Solver_Initialise(CoupledSolver,Err)
  CALL cmfe_SolverEquations_Initialise(CoupledSolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(CoupledProblem,Err)
  !Get the solve equations
  CALL cmfe_Problem_SolverGet(CoupledProblem,CMFE_CONTROL_LOOP_NODE,1,CoupledSolver,Err)
  CALL cmfe_Solver_SolverEquationsGet(CoupledSolver,CoupledSolverEquations,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(CoupledSolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(CoupledSolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)  
  !Add in the first equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(CoupledSolverEquations,EquationsSet1,EquationsSet1Index,Err)
  !Add in the second equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(CoupledSolverEquations,EquationsSet2,EquationsSet2Index,Err)
  !Add in the interface condition
  CALL cmfe_SolverEquations_InterfaceConditionAdd(CoupledSolverEquations,InterfaceCondition,InterfaceConditionIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(CoupledProblem,Err)

  !Start the creation of the equations set boundary conditions for the first equations set
  PRINT *, ' == >> CREATING BOUNDARY CONDITIONS << == '
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(CoupledSolverEquations,BoundaryConditions,Err)
  !Set the first node to 0.0
  FirstNodeNumber=1
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition1,FirstNodeNumber,1,FirstNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
  ENDIF
  !Set boundary conditions for second dependent field
  !Set the last node to 1.0
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Region_NodesGet(Region2,Nodes,Err)
  CALL cmfe_Nodes_NumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition2,LastNodeNumber,1,LastNodeDomain,Err)
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,Err)
  ENDIF
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(CoupledSolverEquations,Err)


  !Solve the problem
  PRINT *, ' == >> SOLVING PROBLEM << == '
  CALL cmfe_Problem_Solve(CoupledProblem,Err)

  !Export the fields
  PRINT *, ' == >> EXPORTING FIELDS << == '
  CALL cmfe_Fields_Initialise(Fields1,Err)
  CALL cmfe_Fields_Create(Region1,Fields1,Err)
  CALL cmfe_Fields_NodesExport(Fields1,"CoupledLaplace_1","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields1,"CoupledLaplace_1","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields1,Err)
  CALL cmfe_Fields_Initialise(Fields2,Err)
  CALL cmfe_Fields_Create(Region2,Fields2,Err)
  CALL cmfe_Fields_NodesExport(Fields2,"CoupledLaplace_2","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields2,"CoupledLaplace_2","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields2,Err)
  CALL cmfe_Fields_Initialise(InterfaceFields,Err)
  CALL cmfe_Fields_Create(INTERFACE,InterfaceFields,Err)
  CALL cmfe_Fields_NodesExport(InterfaceFields,"CoupledLaplace_Interface","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(InterfaceFields,"CoupledLaplace_Interface","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(InterfaceFields,Err)
  
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP
 
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
     
END PROGRAM COUPLEDLAPLACE
