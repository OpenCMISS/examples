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
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 3D
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,3,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)

  !Define basis function - tri-linear Lagrange
  CALL CMISSBasisTypeInitialise(Basis,Err)
  CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
  CALL CMISSBasisTypeSet(Basis,CMISSBasisSimplexType,Err)
  CALL CMISSBasisNumberOfXiSet(Basis,NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(Basis,(/CMISSBasisLinearSimplexInterpolation, &
    & CMISSBasisLinearSimplexInterpolation,CMISSBasisLinearSimplexInterpolation/),Err)
  CALL CMISSBasisCreateFinish(Basis,Err)


  CALL CMISSReadVTK("Atlas_Mesh.vtk", NODE_COORDINATES_C, NUMBER_OF_NODES, NODE_SETS_C, NUMBER_OF_ELEMENTS,Err)
  CALL C_F_POINTER(NODE_COORDINATES_C, NODE_COORDINATES_F90,[NUMBER_OF_NODES*NumberOfMeshDimensions])
  CALL C_F_POINTER(NODE_SETS_C, NODE_SETS_F90,[NUMBER_OF_ELEMENTS*NODES_PER_ELEMENT])

  !Create a mesh
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)
  CALL CMISSMeshNumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  CALL CMISSMeshNumberOfElementsSet(Mesh,NUMBER_OF_ELEMENTS,Err)

  !Define nodes for the mesh
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSNodesCreateStart(Region,NUMBER_OF_NODES,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)

  CALL CMISSMeshElementsTypeInitialise(Elements,Err)
  CALL CMISSMeshElementsCreateStart(Mesh,MeshComponentNumber,Basis,Elements,Err)
  DO N=0,NUMBER_OF_ELEMENTS-1
      INDEX=(NODES_PER_ELEMENT*N)+1
      DO M=1, NODES_PER_ELEMENT
        NODE_SET(M)=NODE_SETS_F90(INDEX+M-1)+1
      ENDDO
      CALL CMISSMeshElementsNodesSet(Elements,N+1,NODE_SET,Err)
  ENDDO


  CALL CMISSMeshElementsCreateFinish(Elements,Err)

  CALL CMISSMeshCreateFinish(Mesh,Err)

!Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,MeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,MeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,MeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  DO N=0,NUMBER_OF_NODES-1
    DO M=0, NumberOfMeshDimensions-1
      CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,N+1 &
      & ,M+1,NODE_COORDINATES_F90(N*NumberOfMeshDimensions+M+1),Err)
    ENDDO
  ENDDO

  !Create the equations_set
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetBioelectricsClass, &
    & CMISSEquationsSetMonodomainEquationType,CMISSEquationsSetNoSubtype,EquationsSetFieldUserNumber,EquationsSetField, &
    & EquationsSet,Err)
  !Set the equations set to be a standard Laplace problem

  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Create the equations set materials field variables
  CALL CMISSFieldTypeInitialise(MaterialsField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,20.00_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,20.00_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,20.00_CMISSDP,Err)

  !Create the equations set source field variables
  CALL CMISSFieldTypeInitialise(SourceField,Err)
  CALL CMISSEquationsSetSourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSetSourceCreateFinish(EquationsSet,Err)

  !Create the CellML environment
  CALL CMISSCellMLTypeInitialise(CellML,Err)
  CALL CMISSCellMLCreateStart(CellMLUserNumber,Region,CellML,Err)
  ! and import LRd from a file
  !CALL CMISSCellMLModelImport(CellML,"LRd.xml",LRdModelIndex,Err)
  CALL CMISSCellMLModelImport(CellML,"FHN.xml",LRdModelIndex,Err)
 !CALL CMISSDiagnosticsSetOn(CMISSInDiagType,(/1,2,3,4,5/),"",(/"CELLML_CREATE_FIELD_TO_CELLML_MAP_C", &
   !& "CELLML_CREATE_CELLML_TO_FIELD_MAP_C"/),Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL CMISSCellMLVariableSetAsKnown(CellML,LRdModelIndex,"Main/I",Err)
  !   - to get from the CellML side
  !CALL CMISSCellMLVariableSetAsWanted(CellML,LRdModelIndex,"Main/alpha",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,LRdModelIndex,"membrane/stim_start",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,LRdModelIndex,"membrane/stim_end",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,LRdModelIndex,"membrane/i_Kp",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,LRdModelIndex,"membrane/i_b",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,LRdModelIndex,"membrane/i_si",Err)

  !   - and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
 !CALL CMISSDiagnosticsSetOff(Err)
  !Finish the CellML environment
  CALL CMISSCellMLCreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellMLFieldMapsCreateStart(CellML,Err)
  !Now we can set up the field variable compo!nent <--> CellML model variable mappings.
  !Map Vm
  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType, &
    & LRdModelIndex,"Main/v",CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,LRdModelIndex,"Main/v",CMISSFieldValuesSetType, &
    & DependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType,Err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellMLFieldMapsCreateFinish(CellML,Err)

  !todo - get vm initialial value.
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,-92.50_CMISSDP,Err)

  !Start the creation of the CellML models field
  CALL CMISSFieldTypeInitialise(CellMLModelsField,Err)
  CALL CMISSCellMLModelsFieldCreateStart(CellMLModelsFieldUserNumber,CellML,CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL CMISSCellMLModelsFieldCreateFinish(CellML,Err)
  !Set up the models field
  !DO N=1,(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  !  IF(N < 5) THEN
  !    CELL_TYPE = 1
  !  ELSE
  !    CELL_TYPE = 2
  !  ENDIF
  !  CALL CMISSFieldParameterSetUpdateNode(CellMLModelsField, CMISSFieldUVariableType, CMISSFieldValuesSetType,1,N,1,CELL_TYPE,Err)
  !END DO
  !CALL CMISSFieldParameterSetUpdateStart(CellMLModelsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  !CALL CMISSFieldParameterSetUpdateFinish(CellMLModelsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  !CALL CMISSFieldComponentValuesInitialise(CellMLModelsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2_CMISSIntg,Err)

  !Start the creation of the CellML state field
  CALL CMISSFieldTypeInitialise(CellMLStateField,Err)
  CALL CMISSCellMLStateFieldCreateStart(CellMLStateFieldUserNumber,CellML,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL CMISSCellMLStateFieldCreateFinish(CellML,Err)

  !!Start the creation of the CellML intermediate field
  !CALL CMISSFieldTypeInitialise(CellMLIntermediateField,Err)
  !CALL CMISSCellMLIntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber,CellML,CellMLIntermediateField,Err)
  !!Finish the creation of the CellML intermediate field
  !CALL CMISSCellMLIntermediateFieldCreateFinish(CellML,Err)

  !Start the creation of CellML parameters field
  CALL CMISSFieldTypeInitialise(CellMLParametersField,Err)
  CALL CMISSCellMLParametersFieldCreateStart(CellMLParametersFieldUserNumber,CellML,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL CMISSCellMLParametersFieldCreateFinish(CellML,Err)

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)

  !Start the creation of the equations set boundary conditions
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSet,BoundaryConditions,Err)
  !Set the first 3175 nodes to flux of 0.0
  DO N=1,3175
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldDelUDelNVariableType,1,1,N,1, &
      & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  ENDDO
  !Finish the creation of the equations set boundary conditions
  CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSet,Err)

  !Set the Stimulus at node 1
  CALL CMISSFieldParameterSetUpdateNode(CellMLParametersField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4000,1, &
    & -10.0_CMISSDP,Err)

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemBioelectricsClass,CMISSProblemMonodomainEquationType, &
    & CMISSProblemMonodomainGudunovSplitSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  !Loop in time for 100 ms with the Stimulus applied.
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,0.0_CMISSDP,0.10_CMISSDP,0.01_CMISSDP,Err)
  !Set the output
  !CALL CMISSControlLoopOutputTypeSet(ControlLoop,CMISSControlLoopTimingOutput,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !Start the creation of the problem solvers
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  !Get the first (DAE) solver
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverDAETimeStepSet(Solver,0.00001_CMISSDP,Err)
  CALL CMISSSolverDAESolverTypeSet(Solver,CMISSSolverDAEExternal,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
  !Get the second (Parabolic) solver
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,2,Solver,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver CellML equations
  CALL CMISSProblemCellMLEquationsCreateStart(Problem,Err)
  !Get the first solver
  !Get the CellML equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSCellMLEquationsTypeInitialise(CellMLEquations,Err)
  CALL CMISSSolverCellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL CMISSCellMLEquationsCellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
  !Finish the creation of the problem solver CellML equations
  CALL CMISSProblemCellMLEquationsCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the second solver
  !Get the solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,2,Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Solve the problem for the first 100 ms
  CALL CMISSProblemSolve(Problem,Err)

  !Now turn the stimulus off
  !Set the Stimulus at node 1
  CALL CMISSFieldParameterSetUpdateNode(CellMLParametersField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4000,1, &
    & 0.0_CMISSDP,Err)

  !Set the time loop for another 900 ms
  CALL CMISSControlLoopTimesSet(ControlLoop,0.1_CMISSDP,1.00_CMISSDP,0.01_CMISSDP,Err)

  !Solve the problem for the next 900 ms
  CALL CMISSProblemSolve(Problem,Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"MonodomainCUDAExample","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"MonodomainCUDAExample","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
  ENDIF
!
!  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM MONODOMAINCUDAEXAMPLE
