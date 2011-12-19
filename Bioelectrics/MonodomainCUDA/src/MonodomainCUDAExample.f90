!> \file
!> $Id: MonodomainExample.f90 1836 2010-12-20 17:25:14Z chrispbradley $
!> \author Chris Bradley
!> \brief This is an example program to solve a Monodomain equation using OpenCMISS calls.
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

!> \example Bioelectrics/MonodomainCUDA/src/MonodomainCUDAExample.f90
!! Example program to solve a Monodomain equation using OpenCMISS and CUDA calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/Bioelectrics/Monodomain/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/Bioelectrics/Monodomain/build-gnu'>Linux GNU Build</a>
!!
!<

!> Main program
PROGRAM MONODOMAINCUDAEXAMPLE

  USE OPENCMISS
  USE MPI
  USE ISO_C_BINDING

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField


  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=3
  INTEGER(CMISSIntg), PARAMETER :: MeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS, ERROR
  ! 4 nodes per tetrahedral element
  INTEGER(CMISSIntg), PARAMETER :: NODES_PER_ELEMENT=4
  INTEGER(CMISSIntg) :: NODE_SET(NODES_PER_ELEMENT)

  ! C Pointer
  TYPE(C_PTR) :: NODE_COORDINATES_C
  TYPE(C_PTR) :: NODE_SETS_C

  ! FORTRAN Pointers to Arrays
  REAL(C_DOUBLE), POINTER ::NODE_COORDINATES_F90(:)
  INTEGER(C_INT), POINTER :: NODE_SETS_F90(:)

  INTEGER(C_INT) :: NUMBER_OF_NODES,NUMBER_OF_ELEMENTS

  INTEGER(CMISSIntg) :: MPI_IERROR

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: N,M,INDEX

  INTEGER(CMISSIntg) :: LRdModelIndex

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCellMLType) :: CellML
  TYPE(CMISSCellMLEquationsType) :: CellMLEquations
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField,SourceField
  TYPE(CMISSFieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSMeshElementsType) :: Elements
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
   !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex,CellMLIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain
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

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  CALL CMISSOutputSetOn("MonodomainCUDA",Err)

  NUMBER_OF_DOMAINS=NumberOfComputationalNodes

  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

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
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Define basis function - tri-linear Lagrange
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(Basis,NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(Basis,(/CMISS_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
    & CMISS_BASIS_LINEAR_SIMPLEX_INTERPOLATION,CMISS_BASIS_LINEAR_SIMPLEX_INTERPOLATION/),Err)
  CALL CMISSBasis_CreateFinish(Basis,Err)


  CALL CMISSReadVTK("Atlas_Mesh.vtk", NODE_COORDINATES_C, NUMBER_OF_NODES, NODE_SETS_C, NUMBER_OF_ELEMENTS,Err)
  CALL C_F_POINTER(NODE_COORDINATES_C, NODE_COORDINATES_F90,[NUMBER_OF_NODES*NumberOfMeshDimensions])
  CALL C_F_POINTER(NODE_SETS_C, NODE_SETS_F90,[NUMBER_OF_ELEMENTS*NODES_PER_ELEMENT])

  !Create a mesh
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  CALL CMISSMesh_NumberOfElementsSet(Mesh,NUMBER_OF_ELEMENTS,Err)

  !Define nodes for the mesh
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(Region,NUMBER_OF_NODES,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)

  CALL CMISSMeshElements_Initialise(Elements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,MeshComponentNumber,Basis,Elements,Err)
  DO N=0,NUMBER_OF_ELEMENTS-1
      INDEX=(NODES_PER_ELEMENT*N)+1
      DO M=1, NODES_PER_ELEMENT
        NODE_SET(M)=NODE_SETS_F90(INDEX+M-1)+1
      ENDDO
      CALL CMISSMeshElements_NodesSet(Elements,N+1,NODE_SET,Err)
  ENDDO


  CALL CMISSMeshElements_CreateFinish(Elements,Err)

  CALL CMISSMesh_CreateFinish(Mesh,Err)

!Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,MeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,MeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,MeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(GeometricField,Err)

  DO N=0,NUMBER_OF_NODES-1
    DO M=0, NumberOfMeshDimensions-1
      CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,N+1 &
      & ,M+1,NODE_COORDINATES_F90(N*NumberOfMeshDimensions+M+1),Err)
    ENDDO
  ENDDO

  !Create the equations_set
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,CMISS_EQUATIONS_SET_BIOELECTRICS_CLASS, &
    & CMISS_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMISS_EQUATIONS_SET_NO_SUBTYPE,EquationsSetFieldUserNumber,EquationsSetField, &
    & EquationsSet,Err)
  !Set the equations set to be a standard Laplace problem

  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the equations set materials field variables
  CALL CMISSField_Initialise(MaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & 20.00_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
    & 20.00_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3, &
    & 20.00_CMISSDP, &
    & Err)

  !Create the equations set source field variables
  CALL CMISSField_Initialise(SourceField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSet_SourceCreateFinish(EquationsSet,Err)

  !Create the CellML environment
  CALL CMISSCellML_Initialise(CellML,Err)
  CALL CMISSCellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  ! and import LRd from a file
  !CALL CMISSCellML_ModelImport(CellML,"LRd.xml",LRdModelIndex,Err)
  CALL CMISSCellML_ModelImport(CellML,"FHN.xml",LRdModelIndex,Err)
 !CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,(/1,2,3,4,5/),"",(/"CELLML_CREATE_FIELD_TO_CELLML_MAP_C", &
   !& "CELLML_CREATE_CELLML_TO_FIELD_MAP_C"/),Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL CMISSCellML_VariableSetAsKnown(CellML,LRdModelIndex,"Main/I",Err)
  !   - to get from the CellML side
  !CALL CMISSCellML_VariableSetAsWanted(CellML,LRdModelIndex,"Main/alpha",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,LRdModelIndex,"membrane/stim_start",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,LRdModelIndex,"membrane/stim_end",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,LRdModelIndex,"membrane/i_Kp",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,LRdModelIndex,"membrane/i_b",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,LRdModelIndex,"membrane/i_si",Err)

  !   - and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
 !CALL CMISSDiagnosticsSetOff(Err)
  !Finish the CellML environment
  CALL CMISSCellML_CreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellML_FieldMapsCreateStart(CellML,Err)
  !Now we can set up the field variable compo!nent <--> CellML model variable mappings.
  !Map Vm
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & LRdModelIndex,"Main/v",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,LRdModelIndex,"Main/v",CMISS_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellML_FieldMapsCreateFinish(CellML,Err)

  !todo - get vm initialial value.
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & -92.50_CMISSDP,Err)

  !Start the creation of the CellML models field
  CALL CMISSField_Initialise(CellMLModelsField,Err)
  CALL CMISSCellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL CMISSCellML_ModelsFieldCreateFinish(CellML,Err)
  !Set up the models field
  !DO N=1,(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  !  IF(N < 5) THEN
  !    CELL_TYPE = 1
  !  ELSE
  !    CELL_TYPE = 2
  !  ENDIF
  !  CALL CMISSField_ParameterSetUpdateNode(CellMLModelsField, CMISS_FIELD_U_VARIABLE_TYPE, CMISS_FIELD_VALUES_SET_TYPE,1,N,1,CELL_TYPE,Err)
  !END DO
  !CALL CMISSField_ParameterSetUpdateStart(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !CALL CMISSField_ParameterSetUpdateFinish(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !CALL CMISSField_ComponentValuesInitialise(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2_CMISSIntg,Err)

  !Start the creation of the CellML state field
  CALL CMISSField_Initialise(CellMLStateField,Err)
  CALL CMISSCellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL CMISSCellML_StateFieldCreateFinish(CellML,Err)

  !!Start the creation of the CellML intermediate field
  !CALL CMISSField_Initialise(CellMLIntermediateField,Err)
  !CALL CMISSCellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  !!Finish the creation of the CellML intermediate field
  !CALL CMISSCellML_IntermediateFieldCreateFinish(CellML,Err)

  !Start the creation of CellML parameters field
  CALL CMISSField_Initialise(CellMLParametersField,Err)
  CALL CMISSCellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL CMISSCellML_ParametersFieldCreateFinish(CellML,Err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Start the creation of the equations set boundary conditions
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSet,BoundaryConditions,Err)
  !Set the first 3175 nodes to flux of 0.0
  DO N=1,3175
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,N,1, &
      & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
  ENDDO
  !Finish the creation of the equations set boundary conditions
  CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSet,Err)

  !Set the Stimulus at node 1
  CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4000,1, &
    & -10.0_CMISSDP,Err)

  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_BIOELECTRICS_CLASS,CMISS_PROBLEM_MONODOMAIN_EQUATION_TYPE, &
    & CMISS_PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  !Loop in time for 100 ms with the Stimulus applied.
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoop_TimesSet(ControlLoop,0.0_CMISSDP,0.10_CMISSDP,0.01_CMISSDP,Err)
  !Set the output
  !CALL CMISSControlLoop_OutputTypeSet(ControlLoop,CMISS_CONTROL_LOOP_TIMING_OUTPUT,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !Start the creation of the problem solvers
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  !Get the first (DAE) solver
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_DAETimeStepSet(Solver,0.00001_CMISSDP,Err)
  CALL CMISSSolver_DAESolverTypeSet(Solver,CMISS_SOLVER_DAE_EXTERNAL,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  !Get the second (Parabolic) solver
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver CellML equations
  CALL CMISSProblem_CellMLEquationsCreateStart(Problem,Err)
  !Get the first solver
  !Get the CellML equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSCellMLEquations_Initialise(CellMLEquations,Err)
  CALL CMISSSolver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL CMISSCellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
  !Finish the creation of the problem solver CellML equations
  CALL CMISSProblem_CellMLEquationsCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !Get the second solver
  !Get the solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_FULL_MATRICES,Err)
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Solve the problem for the first 100 ms
  CALL CMISSProblem_Solve(Problem,Err)

  !Now turn the stimulus off
  !Set the Stimulus at node 1
  CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4000,1, &
    & 0.0_CMISSDP,Err)

  !Set the time loop for another 900 ms
  CALL CMISSControlLoop_TimesSet(ControlLoop,0.1_CMISSDP,1.00_CMISSDP,0.01_CMISSDP,Err)

  !Solve the problem for the next 900 ms
  CALL CMISSProblem_Solve(Problem,Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"MonodomainCUDAExample","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"MonodomainCUDAExample","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF
!
!  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM MONODOMAINCUDAEXAMPLE
