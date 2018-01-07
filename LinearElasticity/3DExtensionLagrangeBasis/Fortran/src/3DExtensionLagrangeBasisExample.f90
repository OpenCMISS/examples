! \file
!> \author Chris Bradley
!> \brief This is an example program to solve a linear elasticity equation using OpenCMISS calls.
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

!> \example LinearElasticity/src/LinearElasticityExample.f90
!! Example program to solve a linear elasticity equation using OpenCMISS calls.
!<

!> Main program
PROGRAM LinearElasticity3DLagrangeBasis
#ifndef NOMPIMOD
  USE MPI
#endif
  USE OpenCMISS
  USE OpenCMISS_Iron

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  REAL(CMISSRP), PARAMETER :: LENGTH=120.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=160.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: HEIGHT=10.0_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Basis1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Basis2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: Basis3UserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfNodes=8
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=3
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfElements=1
  INTEGER(CMISSIntg), PARAMETER :: MeshComponent1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshComponent2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshComponent3UserNumber=3

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariables=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponents=6

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  REAL(CMISSRP), PARAMETER ::   ZERO = 0.0_CMISSRP

  !Program types


  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  LOGICAL :: EXPORT_FIELD

  !CMISS variables

  TYPE(cmfe_RegionType) :: WorldRegion
  TYPE(cmfe_CoordinateSystemType) :: WorldCoordinateSystem
  TYPE(cmfe_BasisType) :: Basis(3)
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField,MaterialField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region
  TYPE(cmfe_SolverType) :: Solver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_MeshElementsType) :: Elements(3)

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

  !Intialise cmiss
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_FINITE_ELEMENT_CALCULATE"],Err)

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberGlobalXElements=1
  NumberGlobalYElements=1
  NumberGlobalZElements=1
  NumberOfDomains=1

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_TypeSet(CoordinateSystem,CMFE_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL cmfe_CoordinateSystem_OriginSet(CoordinateSystem,[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP],Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the CS to the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Define 3 sets of basis functions, one describing each independent coordinate specified by InterpolationType
  !NOTE if you change interpolation you need to change Boundary Conditions
  !NOTE:: Num of Gauss points must be the same across X,Y & Z coordinates and be sufficient to accurately integrate the hightest order interpolation being used

  CALL cmfe_Basis_Initialise(Basis(1),Err)
  CALL cmfe_Basis_CreateStart(Basis1UserNumber,Basis(1),Err)
  CALL cmfe_Basis_TypeSet(Basis(1),CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(Basis(1),NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(Basis(1),[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis(1), &
    & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_CreateFinish(Basis(1),Err)

  CALL cmfe_Basis_Initialise(Basis(2),Err)
  CALL cmfe_Basis_CreateStart(Basis2UserNumber,Basis(2),Err)
  CALL cmfe_Basis_TypeSet(Basis(2),CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(Basis(2),NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(Basis(2),[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis(2), &
    & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_CreateFinish(Basis(2),Err)

  CALL cmfe_Basis_Initialise(Basis(3),Err)
  CALL cmfe_Basis_CreateStart(Basis3UserNumber,Basis(3),Err)
  CALL cmfe_Basis_TypeSet(Basis(3),CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(Basis(3),NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(Basis(3),[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis(3), &
    & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_CreateFinish(Basis(3),Err)

  !Create a mesh
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)

  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,TotalNumberOfElements,Err)

  !Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,TotalNumberOfNodes,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)

  !Create elements for the mesh
  !Mesh Component 1
  CALL cmfe_MeshElements_Initialise(Elements(1),Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,MeshComponent1UserNumber,Basis(1),Elements(1),Err)
  CALL cmfe_MeshElements_NodesSet(Elements(1),1,[1,2,3,4,5,6,7,8],Err)
  CALL cmfe_MeshElements_CreateFinish(Elements(1),Err)
  !Mesh Component 2
  CALL cmfe_MeshElements_Initialise(Elements(2),Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,MeshComponent2UserNumber,Basis(2),Elements(2),Err)
  CALL cmfe_MeshElements_NodesSet(Elements(2),1,[1,2,3,4,5,6,7,8],Err)
  CALL cmfe_MeshElements_CreateFinish(Elements(2),Err)
  !Mesh Component 3
  CALL cmfe_MeshElements_Initialise(Elements(3),Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,MeshComponent3UserNumber,Basis(3),Elements(3),Err)
  CALL cmfe_MeshElements_NodesSet(Elements(3),1,[1,2,3,4,5,6,7,8],Err)
  CALL cmfe_MeshElements_CreateFinish(Elements(3),Err)

  CALL cmfe_Mesh_CreateFinish(Mesh,Err) 

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)  
  CALL cmfe_Field_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,2,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,3,Err)
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Set geometric node coordinates (x)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,1,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,2,1,LENGTH,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,3,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4,1,LENGTH,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,5,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,6,1,LENGTH,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,7,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,8,1,LENGTH,Err)

  !Set geometric node coordinates (y)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,1,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,2,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,3,2,WIDTH,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4,2,WIDTH,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,5,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,6,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,7,2,WIDTH,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,8,2,WIDTH,Err)

  !Set geometric node coordinates (z)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,1,3, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,2,3, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,3,3, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4,3, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,5,3,HEIGHT,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,6,3,HEIGHT,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,7,3,HEIGHT,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,8,3,HEIGHT,Err)

  !Create a dependent field with two variables and three components
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL cmfe_Field_TypeSet(DependentField,CMFE_FIELD_GENERAL_TYPE,Err)  
  CALL cmfe_Field_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField,GeometricField,Err) 
  CALL cmfe_Field_DependentTypeSet(DependentField,CMFE_FIELD_DEPENDENT_TYPE,Err) 
  CALL cmfe_Field_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,2,2,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,3,3,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,2,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,3,Err)
  CALL cmfe_Field_CreateFinish(DependentField,Err)

  !Create a material field and attach it to the geometric field  
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL cmfe_Field_TypeSet(MaterialField,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,2,2,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,3,3,Err)
  CALL cmfe_Field_CreateFinish(MaterialField,Err)

  !Set isotropic elasticity material parameters - Young's Modulus & Poisson's Ratio
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & 10.0E3_CMISSRP, &
    & Err) !E1
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
    & 10.0E3_CMISSRP, &
    & Err) !E2
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
    & 10.0E3_CMISSRP, &
    & Err) !E3
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,0.3_CMISSRP, &
    & Err) !v13
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,0.3_CMISSRP, &
    & Err) !v23
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,0.3_CMISSRP, &
    & Err) !v12

  !Create a Elasticity Class, Linear Elasticity type, no subtype, EquationsSet
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_LINEAR_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(EQUATIONS,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
                                              !CMFE_EQUATIONS_SPARSE_MATRICES=1 !<Use sparse matrices for the equations.
                                              !CMFE_EQUATIONS_FULL_MATRICES=2 !<Use fully populated matrices for the equations. 
  CALL cmfe_Equations_OutputTypeSet(EQUATIONS,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
                                            !CMFE_EQUATIONS_NO_OUTPUT !<No output from the equations.
                                            !CMFE_EQUATIONS_TIMING_OUTPUT !<Timing information output.
                                            !CMFE_EQUATIONS_MATRIX_OUTPUT !<All below and equation matrices output.
                                            !CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT !<All below and Element matrices output.
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)
  
  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_LINEAR_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_NO_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !Start the creation of the Problem solvers
  !Create the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(SOLVER,CMFE_SOLVER_MATRIX_OUTPUT,Err)
                                      !CMFE_SOLVER_NO_OUTPUT !<No output from the solver routines. \see OPENCMISS_SolverOutputTypes,OPENCMISS
                                      !CMFE_SOLVER_PROGRESS_OUTPUT !<Progress output from solver routines.
                                      !CMFE_SOLVER_TIMING_OUTPUT !<Timing output from the solver routines plus below.
                                      !CMFE_SOLVER_SOLVER_OUTPUT !<Solver specific output from the solver routines plus below.
                                      !CMFE_SOLVER_MATRIX_OUTPUT !<Solver matrices output from the solver routines plus below.
  CALL cmfe_Solver_LibraryTypeSet(SOLVER,CMFE_SOLVER_PETSC_LIBRARY,Err)
                                        !CMFE_SOLVER_CMISS_LIBRARY     !<CMISS (internal) solver library.
                                        !CMFE_SOLVER_PETSC_LIBRARY     !<PETSc solver library.
                                        !CMFE_SOLVER_MUMPS_LIBRARY     !<MUMPS solver library.
                                        !CMFE_SOLVER_SUPERLU_LIBRARY   !<SuperLU solver library.
                                        !CMFE_SOLVER_SPOOLES_LIBRARY !<SPOOLES solver library.
                                        !CMFE_SOLVER_UMFPACK_LIBRARY   !<UMFPACK solver library.
                                        !CMFE_SOLVER_LUSOL_LIBRARY     !<LUSOL solver library.
                                        !CMFE_SOLVER_ESSL_LIBRARY      !<ESSL solver library.
                                        !CMFE_SOLVER_LAPACK_LIBRARY    !<LAPACK solver library.
                                        !CMFE_SOLVER_TAO_LIBRARY       !<TAO solver library.
                                        !CMFE_SOLVER_HYPRE_LIBRARY     !<Hypre solver library.
  CALL cmfe_Solver_LinearTypeSet(SOLVER,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
                                      !CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE    !<Direct linear solver type.
                                      !CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE !<Iterative linear solver type.
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)   
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
                                                          !CMFE_SOLVER_SPARSE_MATRICES !<Use sparse solver matrices.
                                                          !CMFE_SOLVER_FULL_MATRICES !<Use fully populated solver matrices.
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Fix nodes 1,3,5,7 at x=0
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED,ZERO,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,3,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED,ZERO,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,5,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED,ZERO,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,7,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED,ZERO,Err)

  !Fix nodes 1,2,5,6 at y=0
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED,ZERO,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,2,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED,ZERO,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,5,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED,ZERO,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,6,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED,ZERO,Err)

  !Fix nodes 1,2,3,4 at x=0
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED,ZERO,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,2,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED,ZERO,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,3,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED,ZERO,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,4,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED,ZERO,Err)

  !Apply force at nodes 1,2,3,4 at x=l
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,2,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & -400.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,4,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & -400.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,6,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & -400.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,8,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & -400.0_CMISSRP,Err)

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !=SOLVE Problem==================================================================================================================
  !Solve the Problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !=OUTPUT SOLUTION================================================================================================================
  !!TODO:: Output reaction forces in ipnode files
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"LinearElasticity3DLinearLagrangeBasisExample","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"LinearElasticity3DLinearLagrangeBasisExample","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM LinearElasticity3DLagrangeBasis


