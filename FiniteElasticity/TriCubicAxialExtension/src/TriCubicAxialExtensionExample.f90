!> \file
!> \author Kumar Mithraratne
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

!> \example FiniteElasticity/TriCubicAxialExtension/src/TriCubicAxialExtensionExample.f90
!! Example program to solve a finite elasticity equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/TriCubicAxialExtension/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/TriCubicAxialExtension/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM TRICUBICAXIALEXTENSIONEXAMPLE

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
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: CubicBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: CubicMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program types

  !Program variables

  INTEGER(CMISSIntg) :: TotalNumberElements,TotalNumberNodes,NumberOfMeshDimensions
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber

  !CMISS variables
  TYPE(cmfe_BasisType) :: CubicBasis, LinearBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,FibreField,MaterialField,DependentField,EquationsSetField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_MeshElementsType) :: CubicElements,LinearElements
  TYPE(cmfe_ControlLoopType) :: ControlLoop

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

  !Set all diganostic levels on for testing
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_RESIDUAL_EVALUATE"],Err)

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !Create a 3D rectangular cartesian coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Define basis functions - tri-linear Lagrange and tri-cubic Lagrange
  CALL cmfe_Basis_Initialise(LinearBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis,[3,3,3],Err)
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)

  CALL cmfe_Basis_Initialise(CubicBasis,Err)
  CALL cmfe_Basis_CreateStart(CubicBasisUserNumber,CubicBasis,Err)
  CALL cmfe_Basis_InterpolationXiSet(CubicBasis,[CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION, &
    & CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(CubicBasis,[3,3,3],Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(CubicBasis,.TRUE.,Err)
  CALL cmfe_Basis_CreateFinish(CubicBasis,Err)

  !Create a mesh with two components, cubic hermite for geometry and linear lagrange
  !for hydrostatic pressure, fibre angles and material properties
  TotalNumberElements=1
  NumberOfMeshDimensions=3
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,2,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,TotalNumberElements,Err)
  !define nodes for the mesh
  TotalNumberNodes=8
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,TotalNumberNodes,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  !cubic Hermite component
  CALL cmfe_MeshElements_Initialise(CubicElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,CubicMeshComponentNumber,CubicBasis,CubicElements,Err)
  CALL cmfe_MeshElements_NodesSet(CubicElements,1,[1,2,3,4,5,6,7,8],Err)
  CALL cmfe_MeshElements_CreateFinish(CubicElements,Err)
  !linear Lagrange component
  CALL cmfe_MeshElements_Initialise(LinearElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,LinearMeshComponentNumber,LinearBasis,LinearElements,Err)
  CALL cmfe_MeshElements_NodesSet(LinearElements,1,[1,2,3,4,5,6,7,8],Err)
  CALL cmfe_MeshElements_CreateFinish(LinearElements,Err)
  !finish mesh creation
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,CubicMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,CubicMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,CubicMeshComponentNumber,Err)
  CALL cmfe_Field_ScalingTypeSet(GeometricField,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Set node parameters
  !node 1
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,1,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,1,1, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,1,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,1,2, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,1,3, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,1,3, &
    & 1.0_CMISSRP,Err)

  !node 2
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,2,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,2,1, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,2,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,2,2, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,2,3, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,2,3, &
    & 1.0_CMISSRP,Err)

  !node 3
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,3,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,3,1, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,3,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,3,2, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,3,3, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,3,3, &
    & 1.0_CMISSRP,Err)

  !node 4
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,4,1, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,4,2, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4,3, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,4,3, &
    & 1.0_CMISSRP,Err)

  !node 5
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,5,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,5,1, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,5,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,5,2, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,5,3, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,5,3, &
    & 1.0_CMISSRP,Err)

  !node 6
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,6,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,6,1, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,6,2, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,6,2, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,6,3, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,6,3, &
    & 1.0_CMISSRP,Err)

  !node 7
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,7,1, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,7,1, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,7,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,7,2, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,7,3, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,7,3, &
    & 1.0_CMISSRP,Err)

  !node 8
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,8,1, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2,8,1, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,8,2, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3,8,2, &
    & 1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,8,3, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,5,8,3, &
    & 1.0_CMISSRP,Err)

  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ScalingTypeSet(FibreField,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !Create the equations_set
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field with 2 variables and 4 components (3 displacement, 1 pressure)
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL cmfe_Field_TypeSet(DependentField,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL cmfe_Field_DependentTypeSet(DependentField,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentField,2,Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CubicMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,2,CubicMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,3,CubicMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CubicMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,CubicMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,CubicMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ScalingTypeSet(DependentField,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_CreateFinish(DependentField,Err)

  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 4.0 respectively.
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,4.0_CMISSRP,Err)

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
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(Solver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(Solver,1.0E-10_CMISSRP,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,1,1,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,1,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,1,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,1,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,1,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,1,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,1,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,1,2,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,1,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,1,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,1,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,1,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,1,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,1,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,1,3,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,1,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,1,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,1,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,1,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,1,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,1,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)

  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,2,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,2,1,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,2,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,2,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,2,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,2,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,2,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,2,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,2,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,2,2,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,2,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,2,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,2,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,2,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,2,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,2,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,2,3,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,2,3,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,2,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,2,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,2,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,2,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,2,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,2,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)

  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,3,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,3,1,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,3,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,3,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,3,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,3,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,3,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,3,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,3,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,3,2,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,3,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,3,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,3,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,3,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,3,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,3,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,3,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,3,3,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,3,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,3,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,3,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,3,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,3,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,3,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)

  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,4,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,4,1,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,4,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,4,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,4,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,4,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,4,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,4,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,4,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,4,2,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,4,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,4,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,4,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,4,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,4,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,4,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,4,3,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,4,3,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,4,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,4,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,4,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,4,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,4,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,4,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)

  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,5,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,5,1,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,5,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,5,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,5,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,5,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,5,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,5,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,5,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,5,2,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,5,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,5,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,5,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,5,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,5,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,5,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,5,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,5,3,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  ! & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,5,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,5,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,5,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,5,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,5,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,5,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)

  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,6,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,6,1,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,6,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,6,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,6,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,6,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,6,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,6,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,6,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,6,2,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,6,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,6,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,6,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,6,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,6,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,6,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,6,3,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,6,3,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,6,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,6,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,6,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,6,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,6,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,6,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)

  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,7,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,7,1,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,7,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,7,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,7,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,7,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,7,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,7,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,7,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,7,2,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,7,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,7,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,7,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,7,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,7,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,7,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,7,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,7,3,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,7,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,7,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,7,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,7,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,7,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,7,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)

  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,8,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,8,1,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  &1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,8,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,8,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,8,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,8,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,8,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,8,1, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,8,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,8,2,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,8,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,8,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,8,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,8,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,8,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,8,2, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,8,3,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSRP,Err)
  !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,8,3,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,3,8,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,4,8,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,5,8,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,6,8,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,7,8,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)
  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,8,8,3, &
    & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSRP,Err)

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !Output solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"TriCubicAxialExtension","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"TriCubicAxialExtension","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program completed."

  STOP

END PROGRAM TRICUBICAXIALEXTENSIONEXAMPLE

