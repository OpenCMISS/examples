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
!> Contributor(s): Jack Lee
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
PROGRAM TRIQUADRATICCUBEEXAMPLE

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
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfNodes=27
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfElements=1

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

  !Program types


  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  !CMISS variables

  TYPE(CMISSBasisType) :: QuadraticBasis, LinearBasis
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
  TYPE(CMISSMeshElementsType) :: QuadraticElements,LinearElements
  TYPE(CMISSControlLoopType) :: ControlLoop

  !Other variables
  INTEGER(CMISSIntg) :: NN,I,J,K,BC_TYPE
  REAL(CMISSDP) :: X,Y,Z

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
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange
  CALL CMISSBasis_Initialise(LinearBasis,Err)
  CALL CMISSBasis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearBasis, &
    & (/CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME/),Err)
  CALL CMISSBasis_CreateFinish(LinearBasis,Err)

  CALL CMISSBasis_Initialise(QuadraticBasis,Err)
  CALL CMISSBasis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasis_InterpolationXiSet(QuadraticBasis,(/CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & (/CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME/),Err)
  CALL CMISSBasis_CreateFinish(QuadraticBasis,Err)

  !Create a mesh with two components, Quadratic for geometry and linear lagrange
  !for hydrostatic pressure, fibre angles and material properties
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  CALL CMISSMesh_NumberOfElementsSet(Mesh,TotalNumberOfElements,Err)

  !define nodes for the mesh
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(Region,TotalNumberOfNodes,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)
  !Quadratic component
  CALL CMISSMeshElements_Initialise(QuadraticElements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,QuadraticMeshComponentNumber,QuadraticBasis,QuadraticElements,Err)
  CALL CMISSMeshElements_NodesSet(QuadraticElements,1, &
    & (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27/),Err) !umm.... 27?
  CALL CMISSMeshElements_CreateFinish(QuadraticElements,Err)
  !linear Lagrange component
  CALL CMISSMeshElements_Initialise(LinearElements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,LinearMeshComponentNumber,LinearBasis,LinearElements,Err)
  CALL CMISSMeshElements_NodesSet(LinearElements,1,(/1,3,7,9,19,21,25,27/),Err) 
  CALL CMISSMeshElements_CreateFinish(LinearElements,Err)

  !finish mesh creation
  CALL CMISSMesh_CreateFinish(Mesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry)
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)  
  CALL CMISSField_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Set node parameters
  NN=0
  DO I=1,3
    Z=(I-1)*0.5_CMISSDP
    DO J=1,3
      Y=(J-1)*0.5_CMISSDP
      DO K=1,3
        X=(K-1)*0.5_CMISSDP
        NN=NN+1
        CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NN,1,X, &
          & Err)
        CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NN,2,Y, &
          & Err)
        CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NN,3,Z, &
          & Err)
      ENDDO
    ENDDO
  ENDDO

  !Create a fibre field and attach it to the geometric field  
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,1,LinearMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,2,LinearMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,3,LinearMeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(FibreField,Err)

  !Create a material field and attach it to the geometric field  
  CALL CMISSField_Initialise(MaterialField,Err)
  CALL CMISSField_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL CMISSField_TypeSet(MaterialField,CMISS_FIELD_MATERIAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(MaterialField,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(MaterialField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)  
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_CreateFinish(MaterialField,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,6.0_CMISSDP,Err)

  !Create the equations_set
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSEquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_NO_SUBTYPE,EquationsSetFieldUserNumber,EquationsSetField,&
    & EquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field with 2 variables and 4 components (3 displacement, 1 pressure)
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSField_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL CMISSField_TypeSet(DependentField,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL CMISSField_DependentTypeSet(DependentField,CMISS_FIELD_DEPENDENT_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL CMISSField_ScalingTypeSet(DependentField,CMISS_FIELD_UNIT_SCALING,Err)
  CALL CMISSField_CreateFinish(DependentField,Err)

!Actually, don't need to create this at all - they're automatically created
! CALL CMISSField_ParameterSetCreate(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_PRESSURE_VALUES_SET_TYPE,Err)

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
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL CMISSControlLoop_MaximumIterationsSet(ControlLoop,3,Err)
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

  !Fix base of the element in z direction
  DO NN=1,9
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NN,3, &
      & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
  ENDDO

  BC_TYPE=2 ! 1=displacement 2=surface pressure 3=force
  SELECT CASE(BC_TYPE)
  CASE(1)
    !1. Move top of the element by a little?
    DO NN=19,27
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NN,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED,1.1_CMISSDP,Err)
    ENDDO
  CASE(2)
    !2. Apply pressure at the top face?
    DO NN=19,27
      !NOTE: Surface pressure goes into pressure_values_set_type of the DELUDELN type
!       CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NN,3, &
!         & CMISS_BOUNDARY_CONDITION_PRESSURE,4.0_CMISSDP,Err)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NN,3, &
        & CMISS_BOUNDARY_CONDITION_PRESSURE_INCREMENTED,4.0_CMISSDP,Err)
    ENDDO
  CASE(3)
    !3. Apply force at the top nodes?
    DO NN=19,27
      !NOTE: Surface pressure goes into pressure_values_set_type of the DELUDELN type
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NN,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED,1.0_CMISSDP,Err)
    ENDDO
  END SELECT

  !Fix two more nodes at the bottom
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,1,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,1,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)
  !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,3,1,CMISS_BOUNDARY_CONDITION_FIXED, &
  !  & 1.0_CMISSDP,Err)
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,3,2, &
    & CMISS_BOUNDARY_CONDITION_FIXED, &
    & 0.0_CMISSDP,Err)

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Output solution  
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"TriQuadraticCube","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"TriQuadraticCube","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM TRIQUADRATICCUBEEXAMPLE

