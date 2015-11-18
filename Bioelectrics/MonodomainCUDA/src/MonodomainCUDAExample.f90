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

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  USE ISO_C_BINDING

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=18
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

  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField,MaterialsField,SourceField
  TYPE(cmfe_FieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_MeshElementsType) :: Elements
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations

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
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  CALL cmfe_OutputSetOn("MonodomainCUDA",Err)

  NUMBER_OF_DOMAINS=NumberOfComputationalNodes

  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 3D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Define basis function - tri-linear Lagrange
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(Basis,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(Basis,[CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
    & CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION,CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],Err)
  CALL cmfe_Basis_CreateFinish(Basis,Err)


  CALL cmfe_ReadVTK("Atlas_Mesh.vtk", NODE_COORDINATES_C, NUMBER_OF_NODES, NODE_SETS_C, NUMBER_OF_ELEMENTS,Err)
  CALL C_F_POINTER(NODE_COORDINATES_C, NODE_COORDINATES_F90,[NUMBER_OF_NODES*NumberOfMeshDimensions])
  CALL C_F_POINTER(NODE_SETS_C, NODE_SETS_F90,[NUMBER_OF_ELEMENTS*NODES_PER_ELEMENT])

  !Create a mesh
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NUMBER_OF_ELEMENTS,Err)

  !Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,NUMBER_OF_NODES,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)

  CALL cmfe_MeshElements_Initialise(Elements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,MeshComponentNumber,Basis,Elements,Err)
  DO N=0,NUMBER_OF_ELEMENTS-1
      INDEX=(NODES_PER_ELEMENT*N)+1
      DO M=1, NODES_PER_ELEMENT
        NODE_SET(M)=NODE_SETS_F90(INDEX+M-1)+1
      ENDDO
      CALL cmfe_MeshElements_NodesSet(Elements,N+1,NODE_SET,Err)
  ENDDO


  CALL cmfe_MeshElements_CreateFinish(Elements,Err)

  CALL cmfe_Mesh_CreateFinish(Mesh,Err)

!Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,MeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,MeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,MeshComponentNumber,Err)
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  DO N=0,NUMBER_OF_NODES-1
    DO M=0, NumberOfMeshDimensions-1
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,N+1 &
      & ,M+1,NODE_COORDINATES_F90(N*NumberOfMeshDimensions+M+1),Err)
    ENDDO
  ENDDO

  !Create the equations_set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
    & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_NO_SUBTYPE],EquationsSetFieldUserNumber,EquationsSetField, &
    & EquationsSet,Err)
  !Set the equations set to be a standard Laplace problem

  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the equations set materials field variables
  CALL cmfe_Field_Initialise(MaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & 20.00_CMISSRP, &
    & Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
    & 20.00_CMISSRP, &
    & Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
    & 20.00_CMISSRP, &
    & Err)

  !Create the equations set source field variables
  CALL cmfe_Field_Initialise(SourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(EquationsSet,Err)

  !Create the CellML environment
  CALL cmfe_CellML_Initialise(CellML,Err)
  CALL cmfe_CellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  ! and import LRd from a file
  !CALL cmfe_CellML_ModelImport(CellML,"LRd.xml",LRdModelIndex,Err)
  CALL cmfe_CellML_ModelImport(CellML,"FHN.xml",LRdModelIndex,Err)
 !CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[1,2,3,4,5],"",["CELLML_CREATE_FIELD_TO_CELLML_MAP_C", &
   !& "CELLML_CREATE_CELLML_TO_FIELD_MAP_C"],Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL cmfe_CellML_VariableSetAsKnown(CellML,LRdModelIndex,"Main/I",Err)
  !   - to get from the CellML side
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,LRdModelIndex,"Main/alpha",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,LRdModelIndex,"membrane/stim_start",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,LRdModelIndex,"membrane/stim_end",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,LRdModelIndex,"membrane/i_Kp",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,LRdModelIndex,"membrane/i_b",Err)
  !CALL cmfe_CellML_VariableSetAsWanted(CellML,LRdModelIndex,"membrane/i_si",Err)

  !   - and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
 !CALL cmfe_DiagnosticsSetOff(Err)
  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(CellML,Err)
  !Now we can set up the field variable compo!nent <--> CellML model variable mappings.
  !Map Vm
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & LRdModelIndex,"Main/v",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,LRdModelIndex,"Main/v",CMFE_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateFinish(CellML,Err)

  !todo - get vm initialial value.
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & -92.50_CMISSRP,Err)

  !Start the creation of the CellML models field
  CALL cmfe_Field_Initialise(CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL cmfe_CellML_ModelsFieldCreateFinish(CellML,Err)
  !Set up the models field
  !DO N=1,(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  !  IF(N < 5) THEN
  !    CELL_TYPE = 1
  !  ELSE
  !    CELL_TYPE = 2
  !  ENDIF
  !  CALL cmfe_Field_ParameterSetUpdateNode(CellMLModelsField, CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE,1,N,1,CELL_TYPE,Err)
  !END DO
  !CALL cmfe_Field_ParameterSetUpdateStart(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  !CALL cmfe_Field_ParameterSetUpdateFinish(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  !CALL cmfe_Field_ComponentValuesInitialise(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2_CMISSIntg,Err)

  !Start the creation of the CellML state field
  CALL cmfe_Field_Initialise(CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL cmfe_CellML_StateFieldCreateFinish(CellML,Err)

  !!Start the creation of the CellML intermediate field
  !CALL cmfe_Field_Initialise(CellMLIntermediateField,Err)
  !CALL cmfe_CellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  !!Finish the creation of the CellML intermediate field
  !CALL cmfe_CellML_IntermediateFieldCreateFinish(CellML,Err)

  !Start the creation of CellML parameters field
  CALL cmfe_Field_Initialise(CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL cmfe_CellML_ParametersFieldCreateFinish(CellML,Err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Set the Stimulus at node 1
  CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4000,1, &
    & -10.0_CMISSRP,Err)

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_BIOELECTRICS_CLASS,CMFE_PROBLEM_MONODOMAIN_EQUATION_TYPE, &
    & CMFE_PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  !Loop in time for 100 ms with the Stimulus applied.
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,0.0_CMISSRP,0.10_CMISSRP,0.01_CMISSRP,Err)
  !Set the output
  !CALL cmfe_ControlLoop_OutputTypeSet(ControlLoop,CMFE_CONTROL_LOOP_TIMING_OUTPUT,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !Start the creation of the problem solvers
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  !Get the first (DAE) solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_DAETimeStepSet(Solver,0.00001_CMISSRP,Err)
  CALL cmfe_Solver_DAESolverTypeSet(Solver,CMFE_SOLVER_DAE_EXTERNAL,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  !Get the second (Parabolic) solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)
  !Get the first solver
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
  !Finish the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the second solver
  !Get the solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

 !Start the creation of the boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the first 3175 nodes to flux of 0.0
  DO N=1,3175
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,N,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
  ENDDO
  !Finish the creation of the boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem for the first 100 ms
  CALL cmfe_Problem_Solve(Problem,Err)

  !Now turn the stimulus off
  !Set the Stimulus at node 1
  CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,4000,1, &
    & 0.0_CMISSRP,Err)

  !Set the time loop for another 900 ms
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,0.1_CMISSRP,1.00_CMISSRP,0.01_CMISSRP,Err)

  !Solve the problem for the next 900 ms
  CALL cmfe_Problem_Solve(Problem,Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"MonodomainCUDAExample","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"MonodomainCUDAExample","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF
!
!  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM MONODOMAINCUDAEXAMPLE
