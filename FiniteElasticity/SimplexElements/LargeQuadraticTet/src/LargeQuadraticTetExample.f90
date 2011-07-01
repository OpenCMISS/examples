!> \file
!> \author Chris Bradley
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
!! Example program to solve a finite elasticity equation using openCMISS calls.
!!
!! \htmlinclude FiniteElasticity/SimplexElements/LargeQuadraticTet/history.html
!<

!> Main program
PROGRAM LARGEQUADRATICTETEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField


  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=1.0_CMISSDP

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
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: DerivativeUserNumber=1

  !Program types


  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  INTEGER(CMISSIntg) :: node_idx,component_idx

  REAL(CMISSDP) :: FibreFieldAngle(3)

  !CMISS variables

  TYPE(CMISSBasisType) :: BasisQuad
  TYPE(CMISSBasisType) :: BasisLin
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,FibreField,MaterialField,DependentField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSMeshElementsType) :: ElementsQuad
  TYPE(CMISSMeshElementsType) :: ElementsLin
  TYPE(CMISSControlLoopType) :: ControlLoop

  !REAL(CMISSDP), POINTER :: FieldData(:)

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

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  CALL CMISSDiagnosticsSetOn(CMISSFromDiagType,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_FINITE_ELEMENT_CALCULATE"/),Err) !CMISSAllDiagType

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

  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemTypeSet(CoordinateSystem,CMISSCoordinateRectangularCartesianType,Err)
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystemOriginSet(CoordinateSystem,(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/),Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Create a region and assign the CS to the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !Define basis function - Simplex tri-quadratic (displacement)
  CALL CMISSBasisTypeInitialise(BasisQuad,Err)
  CALL CMISSBasisCreateStart(BasisQuadserNumberQuad,BasisQuad,Err) 
  CALL CMISSBasisTypeSet(BasisQuad,CMISSBasisSimplexType,Err)
  CALL CMISSBasisNumberOfXiSet(BasisQuad,NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(BasisQuad,(/CMISSBasisQuadraticSimplexInterpolation, &
    & CMISSBasisQuadraticSimplexInterpolation,CMISSBasisQuadraticSimplexInterpolation/),Err)
  CALL CMISSBasisCreateFinish(BasisQuad,Err)

  !Define basis function - Simplex tri-linear (pressure)
  CALL CMISSBasisTypeInitialise(BasisLin,Err)
  CALL CMISSBasisCreateStart(BasisQuadserNumberLin,BasisLin,Err) 
  CALL CMISSBasisTypeSet(BasisLin,CMISSBasisSimplexType,Err)
  CALL CMISSBasisNumberOfXiSet(BasisLin,NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(BasisLin,(/CMISSBasisLinearSimplexInterpolation, &
    & CMISSBasisLinearSimplexInterpolation,CMISSBasisLinearSimplexInterpolation/),Err)
  CALL CMISSBasisCreateFinish(BasisLin,Err)

  !Create a mesh
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)    

  CALL CMISSMeshNumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err) 
  CALL CMISSMeshNumberOfElementsSet(Mesh,TotalNumberOfElements,Err)  

  !Define nodes for the mesh
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSNodesCreateStart(Region,TotalNumberOfNodes,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)

  !Geometry (and displacement) nodes
  CALL CMISSMeshElementsTypeInitialise(ElementsQuad,Err)
  CALL CMISSMeshElementsCreateStart(Mesh,MeshComponentNumberQuad,BasisQuad,ElementsQuad,Err)
  CALL CMISSMeshElementsNodesSet(ElementsQuad,1,(/1,5,2,7,6,3,8,10,9,4/),Err)
  CALL CMISSMeshElementsNodesSet(ElementsQuad,2,(/3,9,4,6,10,2,18,16,14,17/),Err)
  CALL CMISSMeshElementsNodesSet(ElementsQuad,3,(/4,19,21,16,20,17,9,22,18,3/),Err)
  CALL CMISSMeshElementsNodesSet(ElementsQuad,4,(/25,24,3,23,6,2,26,18,14,17/),Err)
  CALL CMISSMeshElementsNodesSet(ElementsQuad,5,(/4,12,13,16,15,17,10,11,14,2/),Err)
  CALL CMISSMeshElementsCreateFinish(ElementsQuad,Err)
  !Pressure nodes
  CALL CMISSMeshElementsTypeInitialise(ElementsLin,Err)
  CALL CMISSMeshElementsCreateStart(Mesh,MeshComponentNumberLin,BasisLin,ElementsLin,Err)
  CALL CMISSMeshElementsNodesSet(ElementsLin,1,(/1,2,3,4/),Err)
  CALL CMISSMeshElementsNodesSet(ElementsLin,2,(/3,4,2,17/),Err)
  CALL CMISSMeshElementsNodesSet(ElementsLin,3,(/4,21,17,3/),Err)
  CALL CMISSMeshElementsNodesSet(ElementsLin,4,(/25,3,2,17/),Err)
  CALL CMISSMeshElementsNodesSet(ElementsLin,5,(/4,13,17,2/),Err)
  CALL CMISSMeshElementsCreateFinish(ElementsLin,Err)

  CALL CMISSMeshCreateFinish(Mesh,Err) 

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry)
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)  
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)  
  !Geometry approximation equivalent to displacement approximation (MeshComponentNumberQuad)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,MeshComponentNumberQuad,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,MeshComponentNumberQuad,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,MeshComponentNumberQuad,Err)
  CALL CMISSFieldVariableLabelSet(GeometricField,CMISSFieldUVariableType,"Geometry",Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !node 1
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,2,0.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,3,0.0_CMISSDP,Err)
  !node 2
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,2,0.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,3,0.0_CMISSDP,Err)
  !node 3
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,2,1.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,3,0.0_CMISSDP,Err)
  !node 4
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,2,0.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,3,1.0_CMISSDP,Err)
  !node 5
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,5,1,0.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,5,2,0.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,5,3,0.0_CMISSDP,Err)
  !node 6
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,6,1,0.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,6,2,0.5_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,6,3,0.0_CMISSDP,Err)
  !node 7
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,7,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,7,2,0.5_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,7,3,0.0_CMISSDP,Err)
  !node 8
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,8,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,8,2,0.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,8,3,0.5_CMISSDP,Err)
  !node 9
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,9,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,9,2,0.5_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,9,3,0.5_CMISSDP,Err)
  !node 10
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,10,1,0.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,10,2,0.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,10,3,0.5_CMISSDP,Err)

  !node 11
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,11,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,11,2,0.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,11,3,0.5_CMISSDP,Err)
  !node 12
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,12,1,0.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,12,2,0.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,12,3,1.0_CMISSDP,Err)
  !node 13
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,13,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,13,2,0.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,13,3,1.0_CMISSDP,Err)
  !node 14
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,14,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,14,2,0.5_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,14,3,0.5_CMISSDP,Err)
  !node 15
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,15,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,15,2,0.5_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,15,3,1.0_CMISSDP,Err)
  !node 16
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,16,1,0.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,16,2,0.5_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,16,3,1.0_CMISSDP,Err)
  !node 17
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,17,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,17,2,1.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,17,3,1.0_CMISSDP,Err)
  !node 18
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,18,1,0.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,18,2,1.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,18,3,0.5_CMISSDP,Err)
  !node 19
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,19,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,19,2,0.5_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,19,3,1.0_CMISSDP,Err)
  !node 20
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,20,1,0.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,20,2,1.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,20,3,1.0_CMISSDP,Err)

  !node 21
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,21,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,21,2,1.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,21,3,1.0_CMISSDP,Err)
  !node 22
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,22,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,22,2,1.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,22,3,0.5_CMISSDP,Err)
  !node 23
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,23,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,23,2,0.5_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,23,3,0.0_CMISSDP,Err)
  !node 24
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,24,1,0.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,24,2,1.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,24,3,0.0_CMISSDP,Err)
  !node 25
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,25,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,25,2,1.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,25,3,0.0_CMISSDP,Err)
  !node 26
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,26,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,26,2,1.0_CMISSDP,Err)  
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,26,3,0.5_CMISSDP,Err)

  !Create a fibre field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(FibreField,Err)
  CALL CMISSFieldCreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSFieldTypeSet(FibreField,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreField,CMISSFieldUVariableType,FieldFibreNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,1,MeshComponentNumberQuad,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,2,MeshComponentNumberQuad,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,3,MeshComponentNumberQuad,Err)
  CALL CMISSFieldVariableLabelSet(FibreField,CMISSFieldUVariableType,"Fibre",Err)
  CALL CMISSFieldCreateFinish(FibreField,Err)

  !Rotation Angles (radians)
  FibreFieldAngle=(/1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP/)
  DO node_idx=1,TotalNumberOfNodes
    DO component_idx=1,FieldFibreNumberOfComponents
      CALL CMISSFieldParameterSetUpdateNode(FibreField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,DerivativeUserNumber, &
        & node_idx,component_idx,FibreFieldAngle(component_idx),Err)
    ENDDO
  ENDDO

  !Create a material field and attach it to the geometric field  
  CALL CMISSFieldTypeInitialise(MaterialField,Err)
  CALL CMISSFieldCreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL CMISSFieldTypeSet(MaterialField,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialField,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(MaterialField,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialField,CMISSFieldUVariableType,FieldMaterialNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(MaterialField,CMISSFieldUVariableType,1,MeshComponentNumberQuad,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialField,CMISSFieldUVariableType,2,MeshComponentNumberQuad,Err)
  CALL CMISSFieldVariableLabelSet(MaterialField,CMISSFieldUVariableType,"Material",Err)
  CALL CMISSFieldCreateFinish(MaterialField,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,6.0_CMISSDP,Err)

  !Create a dependent field with two variables and four components
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSFieldCreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL CMISSFieldTypeSet(DependentField,CMISSFieldGeneralType,Err)  
  CALL CMISSFieldMeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(DependentField,GeometricField,Err) 
  CALL CMISSFieldDependentTypeSet(DependentField,CMISSFieldDependentType,Err) 
  CALL CMISSFieldNumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldUVariableType,FieldDependentNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldDelUDelNVariableType,FieldDependentNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,1,MeshComponentNumberQuad,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,2,MeshComponentNumberQuad,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,3,MeshComponentNumberQuad,Err)  
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,4,MeshComponentNumberLin,Err)  
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,1,MeshComponentNumberQuad,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,2,MeshComponentNumberQuad,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,3,MeshComponentNumberQuad,Err)  
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,4,MeshComponentNumberLin,Err)  
  CALL CMISSFieldVariableLabelSet(DependentField,CMISSFieldUVariableType,"Dependent",Err)
  CALL CMISSFieldCreateFinish(DependentField,Err)  

  !Create the equations_set
    CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
CALL CMISSEquationsSetCreateStart(EquationSetUserNumber,Region,FibreField,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetMooneyRivlinSubtype,EquationsSetFieldUserNumber,EquationsSetField, &
    & EquationsSet,Err)
  
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)   

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,-8.0_CMISSDP,Err)

  !Define the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemElasticityClass,CMISSProblemFiniteElasticityType, &
    & CMISSProblemNoSubtype,Err)
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  CALL CMISSControlLoopMaximumIterationsSet(ControlLoop,2,Err)
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(Solver,CMISSSolverNewtonJacobianFDCalculated,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSSolverNewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearDirectSolveType,Err)
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)   
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Fix nodes 1,3,4,7,8,9,19,21,22 at x=0
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,1,1, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,3,1, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,4,1, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,7,1, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,8,1, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,9,1, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,19,1, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,21,1, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,22,1, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)

  !Fix nodes 2,11,13,14,15,17,23,25,26 at x=1.1
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,2,1, &
    & CMISSBoundaryConditionFixedIncremented, &
    & 1.1_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,11,1, &
    & CMISSBoundaryConditionFixedIncremented, &
    & 1.1_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,13,1, &
    & CMISSBoundaryConditionFixedIncremented, &
    & 1.1_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,14,1, &
    & CMISSBoundaryConditionFixedIncremented, &
    & 1.1_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,15,1, &
    & CMISSBoundaryConditionFixedIncremented, &
    & 1.1_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,17,1, &
    & CMISSBoundaryConditionFixedIncremented, &
    & 1.1_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,23,1, &
    & CMISSBoundaryConditionFixedIncremented, &
    & 1.1_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,25,1, &
    & CMISSBoundaryConditionFixedIncremented, &
    & 1.1_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,26,1, &
    & CMISSBoundaryConditionFixedIncremented, &
    & 1.1_CMISSDP,Err)

  !Fix nodes 1,2,4,5,8,10,11,12,13 at y=0
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,1,2, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,2,2, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,4,2, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,5,2, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,8,2, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,10,2, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,11,2, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,12,2, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,13,2, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)

  !Fix nodes 1,2,3,5,6,7,23,24,25 at z=0
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,1,3, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,2,3, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,3,3, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,5,3, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,6,3, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,7,3, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,23,3, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,24,3, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,25,3, &
    & CMISSBoundaryConditionFixed, &
    & 0.0_CMISSDP,Err)

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblemSolve(Problem,Err)

  !Output solution  
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"LargeQuadraticTet","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"LargeQuadraticTet","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM LARGEQUADRATICTETEXAMPLE

