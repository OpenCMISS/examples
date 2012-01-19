!> \file
!> \author Kumar Mithraratne
!> \brief This is an example program to solve a finite elasticity equation using openCMISS calls.
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
!>

!> \example FiniteElasticity/TriCubicAxialExtension/src/TriCubicAxialExtensionExample.f90
!! Example program to solve a finite elasticity equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/TriCubicAxialExtension/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/TriCubicAxialExtension/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM TRICUBICAXIALEXTENSIONEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

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

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: TotalNumberElements,TotalNumberNodes,NumberOfMeshDimensions
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  !CMISS variables
  TYPE(CMISSBasisType) :: CubicBasis, LinearBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,FibreField,MaterialField,DependentField,EquationsSetField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSMeshElementsType) :: CubicElements,LinearElements
  TYPE(CMISSControlLoopType) :: ControlLoop

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
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  !Set all diganostic levels on for testing
  !CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_RESIDUAL_EVALUATE"],Err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberGlobalXElements=1
  NumberGlobalYElements=1
  NumberGlobalZElements=1
  NumberOfDomains=1

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Create a 3D rectangular cartesian coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Define basis functions - tri-linear Lagrange and tri-cubic Lagrange
  CALL CMISSBasis_Initialise(LinearBasis,Err)
  CALL CMISSBasis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearBasis, &
    & [CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME],Err)
  CALL CMISSBasis_CreateFinish(LinearBasis,Err)

  CALL CMISSBasis_Initialise(CubicBasis,Err)
  CALL CMISSBasis_CreateStart(CubicBasisUserNumber,CubicBasis,Err)
  CALL CMISSBasis_InterpolationXiSet(CubicBasis,[CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION, &
    & CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION,CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION],Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(CubicBasis, &
    & [CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME],Err)
  CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(CubicBasis,.TRUE.,Err)
  CALL CMISSBasis_CreateFinish(CubicBasis,Err)

  !Create a mesh with two components, cubic hermite for geometry and linear lagrange
  !for hydrostatic pressure, fibre angles and material properties
  TotalNumberElements=1
  NumberOfMeshDimensions=3
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,2,Err)
  CALL CMISSMesh_NumberOfElementsSet(Mesh,TotalNumberElements,Err)
  !define nodes for the mesh
  TotalNumberNodes=8
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(Region,TotalNumberNodes,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)
  !cubic Hermite component
  CALL CMISSMeshElements_Initialise(CubicElements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,CubicMeshComponentNumber,CubicBasis,CubicElements,Err)
  CALL CMISSMeshElements_NodesSet(CubicElements,1,[1,2,3,4,5,6,7,8],Err)
  CALL CMISSMeshElements_CreateFinish(CubicElements,Err)
  !linear Lagrange component
  CALL CMISSMeshElements_Initialise(LinearElements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,LinearMeshComponentNumber,LinearBasis,LinearElements,Err)
  CALL CMISSMeshElements_NodesSet(LinearElements,1,[1,2,3,4,5,6,7,8],Err)
  CALL CMISSMeshElements_CreateFinish(LinearElements,Err)
  !finish mesh creation
  CALL CMISSMesh_CreateFinish(Mesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSField_VariableLabelSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,CubicMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,CubicMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,CubicMeshComponentNumber,Err)
  CALL CMISSField_ScalingTypeSet(GeometricField,CMISS_FIELD_UNIT_SCALING,Err)
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Set node parameters
  !node 1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,1,1, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,1,2, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,3, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,1,3, &
    & 1.0_CMISSDP,Err)

  !node 2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,1, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,2,1, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,2,2, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,3, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,2,3, &
    & 1.0_CMISSDP,Err)

  !node 3
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,3,1, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,2, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,3,2, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,3, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,3,3, &
    & 1.0_CMISSDP,Err)

  !node 4
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,1, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,4,1, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,2, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,4,2, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,3, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,4,3, &
    & 1.0_CMISSDP,Err)

  !node 5
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,5,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,5,1, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,5,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,5,2, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,5,3, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,5,3, &
    & 1.0_CMISSDP,Err)

  !node 6
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,6,1, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,6,1, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,6,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,6,2, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,6,3, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,6,3, &
    & 1.0_CMISSDP,Err)

  !node 7
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,7,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,7,1, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,7,2, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,7,2, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,7,3, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,7,3, &
    & 1.0_CMISSDP,Err)

  !node 8
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,8,1, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,8,1, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,8,2, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3,8,2, &
    & 1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,8,3, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,5,8,3, &
    & 1.0_CMISSDP,Err)

  !Create a fibre field and attach it to the geometric field
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_VariableLabelSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,1,LinearMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,2,LinearMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,3,LinearMeshComponentNumber,Err)
  CALL CMISSField_ScalingTypeSet(FibreField,CMISS_FIELD_UNIT_SCALING,Err)
  CALL CMISSField_CreateFinish(FibreField,Err)

  !Create the equations_set
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSEquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EquationsSetFieldUserNumber, &
      & EquationsSetField,&
    & EquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field with 2 variables and 4 components (3 displacement, 1 pressure)
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSField_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL CMISSField_TypeSet(DependentField,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL CMISSField_DependentTypeSet(DependentField,CMISS_FIELD_DEPENDENT_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentField,2,Err)
  CALL CMISSField_VariableLabelSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,CubicMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,2,CubicMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,3,CubicMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CubicMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2,CubicMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3,CubicMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL CMISSField_ScalingTypeSet(DependentField,CMISS_FIELD_UNIT_SCALING,Err)
  CALL CMISSField_CreateFinish(DependentField,Err)

  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL CMISSField_Initialise(MaterialField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 4.0 respectively.
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,4.0_CMISSDP,Err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,-8.0_CMISSDP, &
    & Err)

  !Define the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_ELASTICITY_CLASS,CMISS_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMISS_PROBLEM_NO_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL CMISSControlLoop_MaximumIterationsSet(ControlLoop,2,Err)
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(Solver,1.0E-10_CMISSDP,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,1,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,1,1,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,1,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,1,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,1,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,1,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,1,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,1,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,1,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,1,2,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,1,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,1,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,1,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,1,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,1,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,1,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,1,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,1,3,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,1,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,1,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,1,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,1,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,1,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,1,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,2,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,2,1,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,2,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,2,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,2,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,2,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,2,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,2,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,2,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,2,2,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,2,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,2,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,2,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,2,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,2,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,2,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,2,3,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,2,3,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,2,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,2,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,2,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,2,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,2,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,2,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,3,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,3,1,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,3,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,3,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,3,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,3,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,3,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,3,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,3,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,3,2,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,3,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,3,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,3,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,3,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,3,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,3,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,3,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,3,3,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,3,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,3,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,3,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,3,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,3,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,3,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,4,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,4,1,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,4,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,4,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,4,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,4,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,4,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,4,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,4,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,4,2,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,4,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,4,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,4,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,4,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,4,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,4,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,4,3,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,4,3,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,4,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,4,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,4,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,4,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,4,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,4,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,5,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,5,1,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,5,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,5,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,5,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,5,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,5,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,5,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,5,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,5,2,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,5,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,5,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,5,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,5,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,5,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,5,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,5,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,5,3,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  ! & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,5,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,5,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,5,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,5,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,5,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,5,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,6,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,6,1,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,6,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,6,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,6,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,6,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,6,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,6,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,6,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,6,2,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,6,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,6,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,6,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,6,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,6,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,6,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,6,3,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,6,3,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,6,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,6,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,6,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,6,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,6,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,6,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,7,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,7,1,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,7,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,7,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,7,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,7,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,7,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,7,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,7,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,7,2,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,7,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,7,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,7,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,7,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,7,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,7,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,7,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,7,3,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,7,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,7,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,7,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,7,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,7,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,7,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)

  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,8,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.1_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,8,1,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  &1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,8,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,8,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,8,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,8,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,8,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,8,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,8,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,8,2,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,8,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,8,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,8,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,8,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,8,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,8,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,8,3,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 1.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,8,3,CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
  !  & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,3,8,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4,8,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,5,8,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,6,8,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,7,8,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,8,8,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED, &
    & 0.0_CMISSDP,Err)

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Output solution
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"TriCubicAxialExtension","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"TriCubicAxialExtension","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program completed."

  STOP

END PROGRAM TRICUBICAXIALEXTENSIONEXAMPLE

