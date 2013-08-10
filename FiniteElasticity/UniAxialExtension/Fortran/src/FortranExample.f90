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
!> Contributor(s): Kumar Mithraratne, Adam Reeve
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

!> \example FiniteElasticity/UniAxialExtension/src/UniAxialExtensionExample.f90
!! Example program to solve a finite elasticity equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/UniAxialExtension/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/UniAxialExtension/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM UNIAXIALEXTENSIONEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=1.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfNodes=8
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=1
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfElements=1
  INTEGER(CMISSIntg), PARAMETER :: MeshComponentNumber=1

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
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program types


  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
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
  TYPE(CMISSMeshElementsType) :: Elements

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

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_FINITE_ELEMENT_CALCULATE"/),Err)

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
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_TypeSet(CoordinateSystem,CMISS_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystem_OriginSet(CoordinateSystem,(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/),Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the CS to the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_LabelSet(Region,"Region",Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Define basis function - tri-linear Lagrange  
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err) 
  CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(Basis,NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(Basis,(/CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis, &
    & (/CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME/),Err)  
  CALL CMISSBasis_CreateFinish(Basis,Err)

  !Create a mesh
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)    

  CALL CMISSMesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err) 
  CALL CMISSMesh_NumberOfElementsSet(Mesh,TotalNumberOfElements,Err)  

  !Define nodes for the mesh
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(Region,TotalNumberOfNodes,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)

  CALL CMISSMeshElements_Initialise(Elements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,MeshComponentNumber,Basis,Elements,Err)
  CALL CMISSMeshElements_NodesSet(Elements,1,(/1,2,3,4,5,6,7,8/),Err)
  CALL CMISSMeshElements_CreateFinish(Elements,Err)

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
  CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)  
  CALL CMISSField_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,MeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,MeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,MeshComponentNumber,Err)
  CALL CMISSField_VariableLabelSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !node 1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,3, &
    & 0.0_CMISSDP,Err)
  !node 2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,1, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,3, &
    & 0.0_CMISSDP,Err)
  !node 3
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,2, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,3, &
    & 0.0_CMISSDP,Err)
  !node 4
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,1, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,2, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,3, &
    & 0.0_CMISSDP,Err)
  !node 5
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,5,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,5,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,5,3, &
    & 1.0_CMISSDP,Err)
  !node 6
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,6,1, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,6,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,6,3, &
    & 1.0_CMISSDP,Err)
  !node 7
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,7,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,7,2, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,7,3, &
    & 1.0_CMISSDP,Err)
  !node 8
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,8,1, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,8,2, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,8,3, &
    & 1.0_CMISSDP,Err)

  !Create a fibre field and attach it to the geometric field  
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,1,MeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,2,MeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,3,MeshComponentNumber,Err)
  CALL CMISSField_VariableLabelSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL CMISSField_CreateFinish(FibreField,Err)

  !Create a material field and attach it to the geometric field  
  CALL CMISSField_Initialise(MaterialField,Err)
  CALL CMISSField_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL CMISSField_TypeSet(MaterialField,CMISS_FIELD_MATERIAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(MaterialField,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(MaterialField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,1,MeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,2,MeshComponentNumber,Err)
  CALL CMISSField_VariableLabelSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL CMISSField_CreateFinish(MaterialField,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,6.0_CMISSDP,Err)

  !Create a dependent field with two variables and four components
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSField_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL CMISSField_TypeSet(DependentField,CMISS_FIELD_GENERAL_TYPE,Err)  
  CALL CMISSField_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentField,GeometricField,Err) 
  CALL CMISSField_DependentTypeSet(DependentField,CMISS_FIELD_DEPENDENT_TYPE,Err) 
  CALL CMISSField_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,MeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,2,MeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,3,MeshComponentNumber,Err)  
  CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION, &
    & Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,MeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2,MeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3,MeshComponentNumber,Err)  
  CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
    & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_VariableLabelSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL CMISSField_CreateFinish(DependentField,Err)  

  !Create the equations_set
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSEquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_NO_SUBTYPE,EquationsSetFieldUserNumber,EquationsSetField,&
    & EquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

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
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)   
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Fix nodes 1,3,5,7 at x=0 and nodes 2,4,6,8 at x=1.1
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,1,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,2,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 1.1_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,3,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,4,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 1.1_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,5,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,6,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 1.1_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,7,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,8,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 1.1_CMISSDP,Err)

  !Fix nodes 1,2,5,6 at y=0
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,1,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,2,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,5,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,6,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)

  !Fix nodes 1,2,3,4 at x=0
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,1,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,2,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,3,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,4,3, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Output solution  
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"UniaxialExtension","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"UniaxialExtension","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM UNIAXIALEXTENSIONEXAMPLE

