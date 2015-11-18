!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a finite elasticity equation using OpenCMISS calls.
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
!> Contributor(s): Kumar Mithraratne, Adam Reeve, Thomas Heidlauf
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

!> \example FiniteElasticity/SimplexElements/LargeQuadraticTet/src/LargeQuadraticTetExample.f90
!! Example program to solve a finite elasticity equation using OpenCMISS calls.
!!
!! \htmlinclude FiniteElasticity/SimplexElements/LargeQuadraticTet/history.html
!<

!> Main program
PROGRAM LARGEQUADRATICTETEXAMPLE

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


  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: BasisQuadserNumberQuad=1
  INTEGER(CMISSIntg), PARAMETER :: BasisQuadserNumberLin=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfElements=5
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfNodes=26

  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: MeshComponentNumberQuad=1
  INTEGER(CMISSIntg), PARAMETER :: MeshComponentNumberLin=2

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
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariables=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponents=4

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: DerivativeUserNumber=1

  !Program types


  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  INTEGER(CMISSIntg) :: node_idx,component_idx

  REAL(CMISSRP) :: FibreFieldAngle(3)

  !CMISS variables

  TYPE(cmfe_BasisType) :: BasisQuad
  TYPE(cmfe_BasisType) :: BasisLin
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,FibreField,EquationsSetField,MaterialField,DependentField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_MeshElementsType) :: ElementsQuad
  TYPE(cmfe_MeshElementsType) :: ElementsLin
  TYPE(cmfe_ControlLoopType) :: ControlLoop

  !REAL(CMISSRP), POINTER :: FieldData(:)

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
  CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_FINITE_ELEMENT_CALCULATE"],Err) !CMFE_ALL_DIAG_TYPE

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

  !Define basis function - Simplex tri-quadratic (displacement)
  CALL cmfe_Basis_Initialise(BasisQuad,Err)
  CALL cmfe_Basis_CreateStart(BasisQuadserNumberQuad,BasisQuad,Err) 
  CALL cmfe_Basis_TypeSet(BasisQuad,CMFE_BASIS_SIMPLEX_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(BasisQuad,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(BasisQuad,[CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION,CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION],Err)
  CALL cmfe_Basis_CreateFinish(BasisQuad,Err)

  !Define basis function - Simplex tri-linear (pressure)
  CALL cmfe_Basis_Initialise(BasisLin,Err)
  CALL cmfe_Basis_CreateStart(BasisQuadserNumberLin,BasisLin,Err) 
  CALL cmfe_Basis_TypeSet(BasisLin,CMFE_BASIS_SIMPLEX_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(BasisLin,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(BasisLin,[CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
    & CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION,CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],Err)
  CALL cmfe_Basis_CreateFinish(BasisLin,Err)

  !Create a mesh
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)    

  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err) 
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,TotalNumberOfElements,Err)  

  !Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,TotalNumberOfNodes,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)

  !Geometry (and displacement) nodes
  CALL cmfe_MeshElements_Initialise(ElementsQuad,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,MeshComponentNumberQuad,BasisQuad,ElementsQuad,Err)
  CALL cmfe_MeshElements_NodesSet(ElementsQuad,1,[1,5,2,7,6,3,8,10,9,4],Err)
  CALL cmfe_MeshElements_NodesSet(ElementsQuad,2,[3,9,4,6,10,2,18,16,14,17],Err)
  CALL cmfe_MeshElements_NodesSet(ElementsQuad,3,[4,19,21,16,20,17,9,22,18,3],Err)
  CALL cmfe_MeshElements_NodesSet(ElementsQuad,4,[25,24,3,23,6,2,26,18,14,17],Err)
  CALL cmfe_MeshElements_NodesSet(ElementsQuad,5,[4,12,13,16,15,17,10,11,14,2],Err)
  CALL cmfe_MeshElements_CreateFinish(ElementsQuad,Err)
  !Pressure nodes
  CALL cmfe_MeshElements_Initialise(ElementsLin,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,MeshComponentNumberLin,BasisLin,ElementsLin,Err)
  CALL cmfe_MeshElements_NodesSet(ElementsLin,1,[1,2,3,4],Err)
  CALL cmfe_MeshElements_NodesSet(ElementsLin,2,[3,4,2,17],Err)
  CALL cmfe_MeshElements_NodesSet(ElementsLin,3,[4,21,17,3],Err)
  CALL cmfe_MeshElements_NodesSet(ElementsLin,4,[25,3,2,17],Err)
  CALL cmfe_MeshElements_NodesSet(ElementsLin,5,[4,13,17,2],Err)
  CALL cmfe_MeshElements_CreateFinish(ElementsLin,Err)

  CALL cmfe_Mesh_CreateFinish(Mesh,Err) 

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)  
  CALL cmfe_Field_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  !Geometry approximation equivalent to displacement approximation (MeshComponentNumberQuad)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,MeshComponentNumberQuad,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,MeshComponentNumberQuad,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,MeshComponentNumberQuad,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !node 1
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,1,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,1,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,1,3, &
    & 0.0_CMISSRP,Err)
  !node 2
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,2,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,2,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,2,3, &
    & 0.0_CMISSRP,Err)
  !node 3
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,3,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,3,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,3,3, &
    & 0.0_CMISSRP,Err)
  !node 4
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4,3, &
    & 1.0_CMISSRP,Err)
  !node 5
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,5,1, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,5,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,5,3, &
    & 0.0_CMISSRP,Err)
  !node 6
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,6,1, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,6,2, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,6,3, &
    & 0.0_CMISSRP,Err)
  !node 7
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,7,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,7,2, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,7,3, &
    & 0.0_CMISSRP,Err)
  !node 8
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,8,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,8,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,8,3, &
    & 0.5_CMISSRP,Err)
  !node 9
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,9,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,9,2, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,9,3, &
    & 0.5_CMISSRP,Err)
  !node 10
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,10,1, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,10,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,10,3, &
    & 0.5_CMISSRP,Err)

  !node 11
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,11,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,11,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,11,3, &
    & 0.5_CMISSRP,Err)
  !node 12
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,12,1, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,12,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,12,3, &
    & 1.0_CMISSRP,Err)
  !node 13
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,13,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,13,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,13,3, &
    & 1.0_CMISSRP,Err)
  !node 14
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,14,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,14,2, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,14,3, &
    & 0.5_CMISSRP,Err)
  !node 15
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,15,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,15,2, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,15,3, &
    & 1.0_CMISSRP,Err)
  !node 16
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,16,1, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,16,2, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,16,3, &
    & 1.0_CMISSRP,Err)
  !node 17
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,17,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,17,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,17,3, &
    & 1.0_CMISSRP,Err)
  !node 18
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,18,1, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,18,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,18,3, &
    & 0.5_CMISSRP,Err)
  !node 19
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,19,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,19,2, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,19,3, &
    & 1.0_CMISSRP,Err)
  !node 20
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,20,1, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,20,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,20,3, &
    & 1.0_CMISSRP,Err)

  !node 21
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,21,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,21,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,21,3, &
    & 1.0_CMISSRP,Err)
  !node 22
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,22,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,22,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,22,3, &
    & 0.5_CMISSRP,Err)
  !node 23
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,23,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,23,2, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,23,3, &
    & 0.0_CMISSRP,Err)
  !node 24
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,24,1, &
    & 0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,24,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,24,3, &
    & 0.0_CMISSRP,Err)
  !node 25
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,25,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,25,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,25,3, &
    & 0.0_CMISSRP,Err)
  !node 26
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,26,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,26,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,26,3, &
    & 0.5_CMISSRP,Err)

  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,MeshComponentNumberQuad,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,MeshComponentNumberQuad,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,MeshComponentNumberQuad,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !Rotation Angles (radians)
  FibreFieldAngle=[1.0_CMISSRP,1.0_CMISSRP,1.0_CMISSRP]
  DO node_idx=1,TotalNumberOfNodes
    DO component_idx=1,FieldFibreNumberOfComponents
      CALL cmfe_Field_ParameterSetUpdateNode(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
        & DerivativeUserNumber, &
        & node_idx,component_idx,FibreFieldAngle(component_idx),Err)
    ENDDO
  ENDDO

  !Create a material field and attach it to the geometric field  
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL cmfe_Field_TypeSet(MaterialField,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialField,Decomposition,Err)        
  CALL cmfe_Field_GeometricFieldSet(MaterialField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,1,MeshComponentNumberQuad,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,2,MeshComponentNumberQuad,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL cmfe_Field_CreateFinish(MaterialField,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,6.0_CMISSRP,Err)

  !Create a dependent field with two variables and four components
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL cmfe_Field_TypeSet(DependentField,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)  
  CALL cmfe_Field_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField,GeometricField,Err) 
  CALL cmfe_Field_DependentTypeSet(DependentField,CMFE_FIELD_DEPENDENT_TYPE,Err) 
  CALL cmfe_Field_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,MeshComponentNumberQuad,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,2,MeshComponentNumberQuad,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,3,MeshComponentNumberQuad,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,MeshComponentNumberLin,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,MeshComponentNumberQuad,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,MeshComponentNumberQuad,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,MeshComponentNumberQuad,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,MeshComponentNumberLin,Err)  
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL cmfe_Field_CreateFinish(DependentField,Err)  

  !Create the equations_set
    CALL cmfe_Field_Initialise(EquationsSetField,Err)
CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
  & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],EquationsSetFieldUserNumber, &
  & EquationsSetField,EquationsSet,Err)
  
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)   

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,-8.0_CMISSRP, &
    & Err)

  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_NO_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoop,2,Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(Solver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)   
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Fix nodes 1,3,4,7,8,9,19,21,22 at x=0
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,3,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,4,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,7,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,8,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,9,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,19,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,21,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,22,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)

  !Fix nodes 2,11,13,14,15,17,23,25,26 at x=1.1
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,2,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,11,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,13,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,14,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,15,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,17,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,23,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,25,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,26,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)

  !Fix nodes 1,2,4,5,8,10,11,12,13 at y=0
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,2,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,4,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,5,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,8,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,10,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,11,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,12,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,13,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)

  !Fix nodes 1,2,3,5,6,7,23,24,25 at z=0
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,2,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,3,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,5,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,6,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,7,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,23,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,24,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,25,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSRP,Err)

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !Output solution  
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"LargeQuadraticTet","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"LargeQuadraticTet","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM LARGEQUADRATICTETEXAMPLE

