!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a Laplace equation using OpenCMISS calls.
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

!> \example Meshes/Versions3DLaplace/Versions3DLaplaceExample.f90
!! Example program to solve a Laplace equation over a cube using OpenCMISS calls. Two tricubic hermite elements are used to create the cube using versions to collapse a few nodes. See Example 5i1: Mesh with 2 collapsed elements on the CMISS examples website for more information about the mesh (http://cmiss.bioeng.auckland.ac.nz/development/examples/5/5i/5i1/index.html)
!! \htmlinclude Meshes/Versions3DLaplace/Versions3DLaplace/history.html
!!
!<

!> Main program
PROGRAM VERSIONS3DLAPLACE

  USE OPENCMISS
  USE MPI


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER ::   ZERO = 0.0_CMISSDP
  REAL(CMISSDP), PARAMETER ::   ONE = 1.0_CMISSDP
  REAL(CMISSDP), PARAMETER ::   A = 1.0_CMISSDP
  REAL(CMISSDP), PARAMETER ::   B = 1.0_CMISSDP
  REAL(CMISSDP), PARAMETER ::   C = 0.70710678118654746_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=3
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Basis1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Basis2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: Basis3UserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfNodes=8
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfElements=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshComponent1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshComponent2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshComponent3UserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1
 
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  
  LOGICAL :: EXPORT_FIELD

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis(3)
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,EquationsSetField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSMeshElementsType) :: Elements(3)
  TYPE(CMISSNodesType) :: Nodes

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
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

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  CALL CMISSRandomSeedsSet(9999,Err)
  

  !CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",[&
  !  & "FIELD_GEOMETRIC_PARAMETERS_LINE_LENGTHS_CALCULATE", &
  !  & "DECOMPOSITION_TOPOLOGY_LINES_CALCULATE           " &
  ! & ],Err)

  !CALL CMISSOutputSetOn(Filename,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  NUMBER_OF_DOMAINS=NumberOfComputationalNodes

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 3D
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Start the creation of a basis (Cubic Hermite)
  CALL CMISSBasis_Initialise(Basis(1),Err)
  CALL CMISSBasis_CreateStart(Basis1UserNumber,Basis(1),Err)
  CALL CMISSBasis_TypeSet(Basis(1),CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(Basis(1),NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(Basis(1),(/CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION,CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION, &
	& CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis(1),([4,4,4]),Err)
  CALL CMISSBasis_CreateFinish(Basis(1),Err)

  CALL CMISSBasis_Initialise(Basis(2),Err)
  CALL CMISSBasis_CreateStart(Basis2UserNumber,Basis(2),Err)
  CALL CMISSBasis_TypeSet(Basis(2),CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(Basis(2),NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(Basis(2),(/CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION,CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION, &
	& CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis(2),([4,4,4]),Err)
  CALL CMISSBasis_CreateFinish(Basis(2),Err)

  CALL CMISSBasis_Initialise(Basis(3),Err)
  CALL CMISSBasis_CreateStart(Basis3UserNumber,Basis(3),Err)
  CALL CMISSBasis_TypeSet(Basis(3),CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(Basis(3),NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(Basis(3),(/CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION,CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION, &
	& CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis(3),([4,4,4]),Err)
  CALL CMISSBasis_CreateFinish(Basis(3),Err)
   
  !Create a mesh
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NumberOfXiCoordinates,Mesh,Err)

  CALL CMISSMesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  CALL CMISSMesh_NumberOfElementsSet(Mesh,TotalNumberOfElements,Err)

  !Define nodes for the mesh
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(Region,TotalNumberOfNodes,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)

  CALL CMISSDiagnosticsSetOn(CMISS_ALL_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["FIELD_MAPPINGS_CALCULATE"],Err)

  !Create elements for the mesh
  !Mesh Component 1
  CALL CMISSMeshElements_Initialise(Elements(1),Err)
  CALL CMISSMeshElements_CreateStart(Mesh,MeshComponent1UserNumber,Basis(1),Elements(1),Err)
  !Element 1
  CALL CMISSMeshElements_NodesSet(Elements(1),1,[1,2,3,4,5,6,5,6],Err)
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(1),1,2,5,3,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(1),1,2,5,4,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(1),1,2,5,7,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(1),1,2,5,8,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_NodesSet(Elements(1),2,[3,4,3,4,5,6,7,8],Err)
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(1),2,2,5,1,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(1),2,2,5,2,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(1),2,1,5,3,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(1),2,1,5,4,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(1),2,2,5,5,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(1),2,2,5,6,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_CreateFinish(Elements(1),Err)

  !Mesh Component 2
  CALL CMISSMeshElements_Initialise(Elements(2),Err)
  CALL CMISSMeshElements_CreateStart(Mesh,MeshComponent2UserNumber,Basis(2),Elements(2),Err)
  !Element 1
  CALL CMISSMeshElements_NodesSet(Elements(2),1,[1,2,3,4,5,6,5,6],Err)
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(2),1,2,5,3,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(2),1,2,5,4,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(2),1,2,5,7,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(2),1,2,5,8,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  !Element 2
  CALL CMISSMeshElements_NodesSet(Elements(2),2,[3,4,3,4,5,6,7,8],Err)
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(2),2,2,5,1,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(2),2,2,5,2,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(2),2,1,5,3,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(2),2,1,5,4,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(2),2,2,5,5,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(2),2,2,5,6,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_CreateFinish(Elements(2),Err)

  !Mesh Component 3
  CALL CMISSMeshElements_Initialise(Elements(3),Err)
  CALL CMISSMeshElements_CreateStart(Mesh,MeshComponent3UserNumber,Basis(3),Elements(3),Err)
  !Element1
  CALL CMISSMeshElements_NodesSet(Elements(3),1,[1,2,3,4,5,6,5,6],Err)
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(3),1,2,5,3,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(3),1,2,5,4,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(3),1,2,5,7,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(3),1,2,5,8,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  !Element3
  CALL CMISSMeshElements_NodesSet(Elements(3),2,[3,4,3,4,5,6,7,8],Err)
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(3),2,2,5,1,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(3),2,2,5,2,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(3),2,1,5,3,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(3),2,1,5,4,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(3),2,2,5,5,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_LocalElementNodeVersionSet(Elements(3),2,2,5,6,Err) ! GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex
  CALL CMISSMeshElements_CreateFinish(Elements(3),Err)
  CALL CMISSMesh_CreateFinish(Mesh,Err) 
  
  CALL CMISSDiagnosticsSetOff(Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)
  
  !Destory the mesh now that we have decomposed it
  !CALL CMISSMesh_Destroy(Mesh,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !CALL CMISSField_ScalingTypeSet(GeometricField,CMISS_FIELD_UNIT_SCALING,Err)
  CALL CMISSField_ScalingTypeSet(GeometricField,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)

  CALL CMISSDiagnosticsSetOn(CMISS_ALL_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["FIELD_MAPPINGS_CALCULATE"],Err)

  !Node 1 
  !Geometric x component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,1,ZERO,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,1,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,1,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,1,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,1,1,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,1,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,1,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,1,1,ZERO,Err)
  !Geometric y component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,2,ZERO,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,1,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,1,2,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,1,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,1,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,1,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,1,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,1,2,ZERO,Err)
  !Geometric z component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,3,ZERO,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,1,3,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,1,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,1,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,1,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,1,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,1,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,1,3,ZERO,Err)

  !Node 2
  !Geometric x component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,1,ZERO,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,2,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,2,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,2,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,2,1,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,2,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,2,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,2,1,ZERO,Err)
  !Geometric y component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,2,ZERO,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,2,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,2,2,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,2,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,2,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,2,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,2,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,2,2,ZERO,Err)
  !Geometric z component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,3,ONE,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,2,3,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,2,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,2,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,2,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,2,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,2,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,2,3,ZERO,Err)

  !Node 3
  !Geometric x component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,1,ZERO,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,3,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,3,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,3,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,3,1,ONE,Err) !Version1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,5,3,1,C,Err) !Version2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,3,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,3,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,3,1,ZERO,Err)
  !Geometric y component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,2,ONE,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,3,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,3,2,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,3,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,3,2,ZERO,Err) !Version1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,5,3,2,-C,Err) !Version2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,3,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,3,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,3,2,ZERO,Err)
  !Geometric z component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,3,ZERO,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,3,3,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,3,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,3,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,3,3,ZERO,Err) !Version1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,5,3,3,ZERO,Err) !Version2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,3,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,3,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,3,3,ZERO,Err)

  !Node 4
  !Geometric x component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,1,ZERO,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,4,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,4,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,4,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,4,1,ONE,Err) !Version1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,5,4,1,C,Err) !Version2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,4,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,4,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,4,1,ZERO,Err)
  !Geometric y component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,2,ONE,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,4,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,4,2,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,4,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,4,2,ZERO,Err) !Version1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,5,4,2,-C,Err) !Version2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,4,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,4,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,4,2,ZERO,Err)
  !Geometric z component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,3,ONE,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,4,3,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,4,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,4,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,4,3,ZERO,Err) !Version1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,5,4,3,ZERO,Err) !Version2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,4,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,4,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,4,3,ZERO,Err)

  !Node 5
  !Geometric x component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,5,1,ONE,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,5,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,5,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,5,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,5,1,ONE,Err) !Version1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,5,5,1,C,Err) !Version2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,5,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,5,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,5,1,ZERO,Err)
  !Geometric y component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,5,2,ZERO,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,5,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,5,2,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,5,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,5,2,ZERO,Err) !Version1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,5,5,2,-C,Err) !Version2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,5,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,5,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,5,2,ZERO,Err)
  !Geometric z component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,5,3,ZERO,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,5,3,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,5,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,5,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,5,3,ZERO,Err) !Version1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,5,5,3,ZERO,Err) !Version2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,5,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,5,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,5,3,ZERO,Err)

  !Node 6
  !Geometric x component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,6,1,ONE,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,6,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,6,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,6,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,6,1,ONE,Err) !Version1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,5,6,1,C,Err) !Version2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,6,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,6,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,6,1,ZERO,Err)
  !Geometric y component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,6,2,ZERO,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,6,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,6,2,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,6,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,6,2,ZERO,Err) !Version1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,5,6,2,-C,Err) !Version2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,6,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,6,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,6,2,ZERO,Err)
  !Geometric z component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,6,3,ONE,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,6,3,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,6,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,6,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,6,3,ZERO,Err) !Version1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,5,6,3,ZERO,Err) !Version2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,6,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,6,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,6,3,ZERO,Err)

  !Node 7
  !Geometric x component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,7,1,ONE,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,7,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,7,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,7,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,7,1,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,7,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,7,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,7,1,ZERO,Err)
  !Geometric y component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,7,2,ONE,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,7,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,7,2,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,7,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,7,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,7,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,7,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,7,2,ZERO,Err)
  !Geometric z component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,7,3,ZERO,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,7,3,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,7,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,7,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,7,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,7,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,7,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,7,3,ZERO,Err)

  !Node 8
  !Geometric x component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,8,1,ONE,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,8,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,8,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,8,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,8,1,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,8,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,8,1,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,8,1,ZERO,Err)
  !Geometric y component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,8,2,ONE,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,8,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,8,2,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,8,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,8,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,8,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,8,2,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,8,2,ZERO,Err)
  !Geometric z component
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,8,3,ONE,Err) !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,8,3,ONE,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,8,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4,8,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,8,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,6,8,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,7,8,3,ZERO,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,8,8,3,ZERO,Err)

  CALL CMISSField_ParameterSetUpdateStart(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  CALL CMISSDiagnosticsSetOff(Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"InitialGeometry","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"InitialGeometry","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF

  !Create the Standard Laplace Equations set
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMISS_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Set the DOFs to be contiguous across components
  !CALL CMISSField_DOFOrderTypeSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  !CALL CMISSField_DOFOrderTypeSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Initialise the field with an initial guess
  !CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.5_CMISSDP,Err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  !CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)
  
  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS,CMISS_PROBLEM_LAPLACE_EQUATION_TYPE, &
    & CMISS_PROBLEM_STANDARD_LAPLACE_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)
 
  !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
  CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-12_CMISSDP,Err)
  CALL CMISSSolver_LinearIterativeRelativeToleranceSet(Solver,1.0E-12_CMISSDP,Err)
  !CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  !CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Start the creation of the equations set boundary conditions
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  LastNodeNumber=8
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED,1.0_CMISSDP,Err)
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblem_Solve(Problem,Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"LaplaceSolution","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"LaplaceSolution","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM VERSIONS3DLAPLACE
