!> \file
!> \author Sebastian Krittian
!> \brief This is an example program to solve a static Navier-Stokes equation using OpenCMISS calls.
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

!> \example FluidMechanics/NavierStokes/Static_FieldML/src/Static_FieldMLExample.f90
!! Example program to solve a static Navier-Stokes equation using OpenCMISS calls.
!! \htmlinclude FluidMechanics/NavierStokes/Static_FieldML/history.html
!!
!<

!> Main program

PROGRAM NAVIERSTOKESSTATICEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif
  USE FIELDML_API

#ifdef WIN32
  USE IFQWINCMISS
#endif

  !
  !================================================================================================================================
  !

  !PROGRAM VARIABLES AND TYPES

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberNavierStokes=7
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokes=8
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberNavierStokes=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberNavierStokes=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverNavierStokesUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesMu=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesRho=2

  INTEGER(CMISSIntg), PARAMETER :: basisNumberTrilinear=1
  INTEGER(CMISSIntg), PARAMETER :: basisNumberTriquadratic=2

  INTEGER(CMISSIntg), PARAMETER :: gaussQuadrature(3) = [3,3,3]

  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: inputFilename = "../../TOOLS/NavierStokesMeshes/HEX-M2-V2-P1_FE.xml"

  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputDirectory = "output"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputFilename = outputDirectory//"/HEX-M2-V2-P1_FE_out.xml"

  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: basename = "static_navier_stokes"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: dataFormat = "PLAIN_TEXT"
  !Program types

  !Program variables

  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_PRESSURE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_SPACE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_NAVIER_STOKES
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES

  INTEGER(CMISSIntg) :: EQUATIONS_NAVIER_STOKES_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: NODE_COUNTER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_NAVIER_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_NAVIER_STOKES

  REAL(CMISSRP) :: INITIAL_FIELD_NAVIER_STOKES(3)
  REAL(CMISSRP) :: BOUNDARY_CONDITIONS_NAVIER_STOKES(3)
  REAL(CMISSRP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSRP) :: RELATIVE_TOLERANCE
  REAL(CMISSRP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSRP) :: LINESEARCH_ALPHA
  REAL(CMISSRP) :: VALUE
  REAL(CMISSRP) :: MU_PARAM_NAVIER_STOKES
  REAL(CMISSRP) :: RHO_PARAM_NAVIER_STOKES

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG
  LOGICAL :: FIXED_WALL_NODES_NAVIER_STOKES_FLAG
  LOGICAL :: INLET_WALL_NODES_NAVIER_STOKES_FLAG

  !CMISS variables

  !Regions
  TYPE(cmfe_RegionType) :: Region
  TYPE(cmfe_RegionType) :: WorldRegion
  !Coordinate systems
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
  TYPE(cmfe_CoordinateSystemType) :: WorldCoordinateSystem
  !Nodes
  TYPE(cmfe_NodesType) :: Nodes
  !Meshes
  TYPE(cmfe_MeshType) :: Mesh
  !Decompositions
  TYPE(cmfe_DecompositionType) :: Decomposition
  !Fields
  TYPE(cmfe_FieldsType) :: Fields
  !Field types
  TYPE(cmfe_FieldType) :: GeometricField
  TYPE(cmfe_FieldType) :: EquationsSetField
  TYPE(cmfe_FieldType) :: DependentFieldNavierStokes
  TYPE(cmfe_FieldType) :: MaterialsFieldNavierStokes
  !Boundary conditions
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsNavierStokes
  !Equations sets
  TYPE(cmfe_EquationsSetType) :: EquationsSetNavierStokes
  !Equations
  TYPE(cmfe_EquationsType) :: EquationsNavierStokes
  !Problems
  TYPE(cmfe_ProblemType) :: Problem
  !Control loops
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  !Solvers
  TYPE(cmfe_SolverType) :: NonlinearSolverNavierStokes
  TYPE(cmfe_SolverType) :: LinearSolverNavierStokes
  !Solver equations
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsNavierStokes
  
  !FieldML parsing variables
  TYPE(cmfe_FieldMLIOType) :: fieldmlInfo, outputInfo
  
  INTEGER(CMISSIntg) :: meshComponentCount
  
  INTEGER(CMISSIntg) :: typeHandle
  INTEGER(CMISSIntg) :: coordinateCount
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables

  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: err

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

  !
  !================================================================================================================================
  !

  !PROBLEM CONTROL PANEL
  
  !Set initial values
  INITIAL_FIELD_NAVIER_STOKES(1)=0.0_CMISSRP
  INITIAL_FIELD_NAVIER_STOKES(2)=0.0_CMISSRP
  INITIAL_FIELD_NAVIER_STOKES(3)=0.0_CMISSRP
  !Set boundary conditions
  FIXED_WALL_NODES_NAVIER_STOKES_FLAG=.TRUE.
  INLET_WALL_NODES_NAVIER_STOKES_FLAG=.TRUE.
  IF(FIXED_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    NUMBER_OF_FIXED_WALL_NODES_NAVIER_STOKES=80
    ALLOCATE(FIXED_WALL_NODES_NAVIER_STOKES(NUMBER_OF_FIXED_WALL_NODES_NAVIER_STOKES))
    FIXED_WALL_NODES_NAVIER_STOKES=[1,2,3,4,5,7,9,10,11,12,13,14,17,20,24,28,29,30,31,32,33,34,35,37,39, & 
    & 41,44,46,47,48,50,51,52,53,54,57,60,64,65,66,67,68,70,72,74,76,77,78,79,80,83,86, & 
    & 89,90,91,92,93,94,95,97,99,101,102,103,104,105,106,107,108,111,114,115,116,117,118, & 
    & 120,122,123,124,125]
  ENDIF
  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES=9
    ALLOCATE(INLET_WALL_NODES_NAVIER_STOKES(NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES))
    INLET_WALL_NODES_NAVIER_STOKES=[6,15,16,23,36,42,81,82,96]
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_NAVIER_STOKES(1)=0.0_CMISSRP
    BOUNDARY_CONDITIONS_NAVIER_STOKES(2)=1.0_CMISSRP
    BOUNDARY_CONDITIONS_NAVIER_STOKES(3)=0.0_CMISSRP
  ENDIF
  !Set material parameters
  MU_PARAM_NAVIER_STOKES=1.0_CMISSRP
  RHO_PARAM_NAVIER_STOKES=1.0_CMISSRP
  !Set interpolation parameters
  BASIS_XI_GAUSS_SPACE=3
  BASIS_XI_GAUSS_VELOCITY=3
  BASIS_XI_GAUSS_PRESSURE=3
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMFE_SOLVER_NO_OUTPUT
  NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMFE_SOLVER_NO_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_NAVIER_STOKES_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT
  !Set solver parameters
  LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-05_CMISSRP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-10_CMISSRP
  DIVERGENCE_TOLERANCE=1.0E20 !default: 1.0E5
  MAXIMUM_ITERATIONS=100000 !default: 100000
  RESTART_VALUE=3000 !default: 30
  LINESEARCH_ALPHA=1.0


  !
  !================================================================================================================================
  !

  !INITIALISE OPENCMISS

  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion, err )

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  !
  !================================================================================================================================
  !

  CALL cmfe_FieldMLIO_Initialise( fieldmlInfo, err ) 
  CALL cmfe_FieldML_InputCreateFromFile( inputFilename, fieldmlInfo, err )

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise( CoordinateSystem, err )
  CALL cmfe_FieldML_InputCoordinateSystemCreateStart( fieldmlInfo, "test_mesh.coordinates", CoordinateSystem, &
    & CoordinateSystemUserNumber, err )
  CALL cmfe_CoordinateSystem_CreateFinish( CoordinateSystem, err )
  CALL cmfe_CoordinateSystem_DimensionGet( CoordinateSystem, coordinateCount, err )

  !
  !================================================================================================================================
  !

  !REGION

  !Start the creation of a new region
  CALL cmfe_Region_Initialise( Region, err )
  CALL cmfe_Region_CreateStart( RegionUserNumber, WorldRegion, Region, err )
  !Set the regions coordinate system as defined above
  CALL cmfe_Region_CoordinateSystemSet( Region, CoordinateSystem, err )
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish( Region, err )

  !
  !================================================================================================================================
  !

  !NODES
  CALL cmfe_FieldML_InputNodesCreateStart( fieldmlInfo, "test_mesh.nodes.argument", Region, nodes, err )
  CALL cmfe_Nodes_CreateFinish( Nodes, err )

  !
  !================================================================================================================================
  !

  !BASES
  CALL cmfe_FieldML_InputBasisCreateStart( fieldmlInfo, "test_mesh.trilinear_lagrange", basisNumberTrilinear, err )
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet( basisNumberTrilinear, gaussQuadrature, err )
  CALL cmfe_Basis_CreateFinish( basisNumberTrilinear, err )

  CALL cmfe_FieldML_InputBasisCreateStart( fieldmlInfo, "test_mesh.triquadratic_lagrange", basisNumberTriquadratic, err )
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet( basisNumberTriquadratic, gaussQuadrature, err )
  CALL cmfe_Basis_CreateFinish( basisNumberTriquadratic, err )
  
  !
  !================================================================================================================================
  !

  !MESH
  
  meshComponentCount = 2

  CALL cmfe_FieldML_InputMeshCreateStart( fieldmlInfo, "test_mesh.mesh.argument", Mesh, MeshUserNumber, Region, err )
  CALL cmfe_Mesh_NumberOfComponentsSet( Mesh, meshComponentCount, err )
  
  CALL cmfe_FieldML_InputCreateMeshComponent( fieldmlInfo, RegionUserNumber, MeshUserNumber, 1, &
    & "test_mesh.template.triquadratic", err )
  CALL cmfe_FieldML_InputCreateMeshComponent( fieldmlInfo, RegionUserNumber, MeshUserNumber, 2, &
    & "test_mesh.template.trilinear", err )
  
  MESH_COMPONENT_NUMBER_SPACE = 1
  MESH_COMPONENT_NUMBER_VELOCITY = 1
  MESH_COMPONENT_NUMBER_PRESSURE = 2

  !Finish the creation of the mesh
  CALL cmfe_Mesh_CreateFinish(Mesh, err )

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition, err )
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition, err )
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE, err )
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,DomainUserNumber, err )
  CALL cmfe_Decomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition, err )
  
  CALL cmfe_FieldML_InputFieldCreateStart( fieldmlInfo, Region, Decomposition, GeometricFieldUserNumber, GeometricField, &
    & CMFE_FIELD_U_VARIABLE_TYPE, "test_mesh.coordinates", err )
  CALL cmfe_Field_CreateFinish( RegionUserNumber, GeometricFieldUserNumber, err )

  CALL cmfe_FieldML_InputFieldParametersUpdate( fieldmlInfo, GeometricField, "test_mesh.node.coordinates", &
    & CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err )
  CALL cmfe_Field_ParameterSetUpdateStart( GeometricField, CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err )
  CALL cmfe_Field_ParameterSetUpdateFinish( GeometricField, CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err )

  !
  !================================================================================================================================
  !

  CALL cmfe_FieldMLIO_Finalise( fieldmlInfo, err )

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for static Navier-Stokes
  CALL cmfe_EquationsSet_Initialise(EquationsSetNavierStokes, err )
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField, &
    & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMFE_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE],EquationsSetFieldUserNumber,EquationsSetField,EquationsSetNavierStokes,err)
  !Set the equations set to be a static Navier-Stokes problem
  
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetNavierStokes, err )


  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for static Navier-Stokes
  CALL cmfe_Field_Initialise(DependentFieldNavierStokes, err )
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetNavierStokes,DependentFieldUserNumberNavierStokes, & 
    & DependentFieldNavierStokes, err )
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,coordinateCount
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY, err )
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY, err )
  ENDDO
  COMPONENT_NUMBER=coordinateCount+1
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE, err )
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE, err )
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetNavierStokes, err )

  !Initialise dependent field
  DO COMPONENT_NUMBER=1,coordinateCount
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & COMPONENT_NUMBER,INITIAL_FIELD_NAVIER_STOKES(COMPONENT_NUMBER), err )
  ENDDO


  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for static Navier-Stokes
  CALL cmfe_Field_Initialise(MaterialsFieldNavierStokes, err )
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetNavierStokes,MaterialsFieldUserNumberNavierStokes, & 
    & MaterialsFieldNavierStokes, err )
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetNavierStokes, err )
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberNavierStokesMu,MU_PARAM_NAVIER_STOKES, err )
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberNavierStokesRho,RHO_PARAM_NAVIER_STOKES, err )


  !
  !================================================================================================================================
  !

  !EQUATIONS


  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsNavierStokes, err )
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetNavierStokes,EquationsNavierStokes, err )
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(EquationsNavierStokes,CMFE_EQUATIONS_SPARSE_MATRICES, err )
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsNavierStokes,EQUATIONS_NAVIER_STOKES_OUTPUT, err )
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetNavierStokes, err )
  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem, err )
  CALL cmfe_ControlLoop_Initialise(ControlLoop, err )
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS,CMFE_PROBLEM_NAVIER_STOKES_EQUATION_TYPE, &
    & CMFE_PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE],Problem,err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem, err )
  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem, err )
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem, err )

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(NonlinearSolverNavierStokes, err )
  CALL cmfe_Solver_Initialise(LinearSolverNavierStokes, err )
  CALL cmfe_Problem_SolversCreateStart(Problem, err )
  !Get the nonlinear static solver
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,NonlinearSolverNavierStokes, err )
  !Set the nonlinear Jacobian type
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED, &
    &  err )
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(NonlinearSolverNavierStokes,NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE, err )
  !Set the solver settings
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(NonlinearSolverNavierStokes,ABSOLUTE_TOLERANCE, err )
  CALL cmfe_Solver_NewtonRelativeToleranceSet(NonlinearSolverNavierStokes,RELATIVE_TOLERANCE, err )
  !Get the nonlinear linear solver
  CALL cmfe_Solver_NewtonLinearSolverGet(NonlinearSolverNavierStokes,LinearSolverNavierStokes, err )
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(LinearSolverNavierStokes,LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE, err )


  !Set the solver settings
  IF(LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolverNavierStokes,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE, err )
    CALL cmfe_Solver_LibraryTypeSet(LinearSolverNavierStokes,CMFE_SOLVER_MUMPS_LIBRARY, err )
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolverNavierStokes,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE, err )
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolverNavierStokes,MAXIMUM_ITERATIONS, err )
    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolverNavierStokes,DIVERGENCE_TOLERANCE, err )
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolverNavierStokes,RELATIVE_TOLERANCE, err )
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolverNavierStokes,ABSOLUTE_TOLERANCE, err )
    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(LinearSolverNavierStokes,RESTART_VALUE, err )
  ENDIF
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem, err )

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(LinearSolverNavierStokes, err )
  CALL cmfe_SolverEquations_Initialise(SolverEquationsNavierStokes, err )
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem, err )
  !Get the linear solver equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,LinearSolverNavierStokes, err )
  CALL cmfe_Solver_SolverEquationsGet(LinearSolverNavierStokes,SolverEquationsNavierStokes, err )
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsNavierStokes,CMFE_SOLVER_SPARSE_MATRICES, err )
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsNavierStokes,EquationsSetNavierStokes,EquationsSetIndex, err )
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem, err )


  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  !Start the creation of the equations set boundary conditions for Navier-Stokes
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsNavierStokes, err )
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsNavierStokes,BoundaryConditionsNavierStokes, err )
  !Set fixed wall nodes
  IF(FIXED_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_FIXED_WALL_NODES_NAVIER_STOKES
      NODE_NUMBER=FIXED_WALL_NODES_NAVIER_STOKES(NODE_COUNTER)
      CONDITION=CMFE_BOUNDARY_CONDITION_FIXED_WALL
      DO COMPONENT_NUMBER=1,coordinateCount
        VALUE=0.0_CMISSRP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
          & CMFE_FIELD_U_VARIABLE_TYPE, &
          & 1, &
          & CMFE_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE, err )
      ENDDO
    ENDDO
  ENDIF
  !Set velocity boundary conditions
  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES
      NODE_NUMBER=INLET_WALL_NODES_NAVIER_STOKES(NODE_COUNTER)
      CONDITION=CMFE_BOUNDARY_CONDITION_FIXED_INLET
      DO COMPONENT_NUMBER=1,coordinateCount
        VALUE=BOUNDARY_CONDITIONS_NAVIER_STOKES(COMPONENT_NUMBER)
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
          & CMFE_FIELD_U_VARIABLE_TYPE, &
          & 1, &
          & CMFE_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE, err )
      ENDDO
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsNavierStokes, err )

  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn of PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL cmfe_Problem_Solve(Problem, err )
  WRITE(*,'(A)') "Problem solved!"

  !
  !================================================================================================================================
  !

  !OUTPUT
  CALL cmfe_FieldMLIO_Initialise( outputInfo, err )

  CALL cmfe_FieldML_OutputCreate( Mesh, outputDirectory, basename, dataFormat, outputInfo, err )
  
  CALL cmfe_FieldML_OutputAddImport( outputInfo, "coordinates.rc.3d", typeHandle, err )
  CALL cmfe_FieldML_OutputAddField( outputInfo, baseName//".geometric", dataFormat, GeometricField, &
    & CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err )
  
  CALL cmfe_FieldML_OutputAddFieldComponents( outputInfo, typeHandle, baseName//".velocity", dataFormat, &
    & DependentFieldNavierStokes, [1,2,3], CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err )
  
  CALL cmfe_FieldML_OutputAddImport( outputInfo, "real.1d", typeHandle, err )
  CALL cmfe_FieldML_OutputAddFieldComponents( outputInfo, typeHandle, baseName//".pressure", dataFormat, &
    & DependentFieldNavierStokes, [4], CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err )
  
  CALL cmfe_FieldML_OutputWrite( outputInfo, outputFilename, err )
  
  CALL cmfe_FieldMLIO_Finalise( outputInfo, err )
  
  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL cmfe_Fields_Initialise(Fields, err )
    CALL cmfe_Fields_Create(Region,Fields, err )
    CALL cmfe_Fields_NodesExport(Fields,"StaticNavierStokes","FORTRAN", err )
    CALL cmfe_Fields_ElementsExport(Fields,"StaticNavierStokes","FORTRAN", err )
    CALL cmfe_Fields_Finalise(Fields, err )
    WRITE(*,'(A)') "Field exported!"
  ENDIF

  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM NAVIERSTOKESSTATICEXAMPLE
