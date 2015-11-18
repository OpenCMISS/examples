!> \file
!> \author David Ladd
!> \brief This is an example program to numerically solve a 2D dynamic 
!>  nonlinear Taylor-Green vortex problem with the Navier-Stokes equations 
!>  and validate against the exact analytic form using OpenCMISS calls.
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

!> \example FluidMechanics/NavierStokes/TaylorGreen/src/TaylorGreenExample.f90
!! Example program to solve a dynamic Taylor Green vortex using OpenCMISS calls.
!<

!> Main program
PROGRAM TAYLORGREENEXAMPLE

  ! PROGRAM LIBRARIES
  USE OpenCMISS
  USE OpenCMISS_Iron
  USE FLUID_MECHANICS_IO_ROUTINES
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

!==============================================================


  REAL(CMISSRP), PARAMETER :: U_0=1.0_CMISSRP !Characteristic velocity (velocity boundary amplitude at t=0)
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP !Length scale (length of one side of the square)
  REAL(CMISSRP), PARAMETER :: VISCOSITY=0.01_CMISSRP
  REAL(CMISSRP), PARAMETER :: DENSITY=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: START_TIME=0.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: STOP_TIME=1.0001_CMISSRP
  REAL(CMISSRP), PARAMETER :: TIME_INCREMENT=0.01_CMISSRP
  REAL(CMISSRP), PARAMETER :: THETA=1.0_CMISSRP
  INTEGER(CMISSIntg), PARAMETER :: OUTPUT_FREQUENCY=1

!==============================================================


  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberNavierStokes=8
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberNavierStokes=9
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokes=10
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberNavierStokes=12
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=15

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: SolverNavierStokesUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesMu=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesRho=2

  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_DIMENSIONS,INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS
  INTEGER(CMISSIntg) :: component_idx
  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  REAL(CMISSRP) :: POSITION,TOTAL_LENGTH
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT

  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber,FirstNodeDomain,LastNodeDomain

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_SPACE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_PRESSURE
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_SPACE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_PRESSURE
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
!   INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NUMBER_OF_BOUNDARY_NODES_NAVIER_STOKES

  INTEGER(CMISSIntg) :: EQUATIONS_NAVIER_STOKES_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_COUNTER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER, ALLOCATABLE, DIMENSION(:):: BOUNDARY_NODES_NAVIER_STOKES

  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE

  REAL(CMISSRP) :: INITIAL_FIELD_NAVIER_STOKES(2)
  REAL(CMISSRP) :: BOUNDARY_CONDITIONS_NAVIER_STOKES(2)
  REAL(CMISSRP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSRP) :: RELATIVE_TOLERANCE
  REAL(CMISSRP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSRP) :: LINESEARCH_ALPHA
  REAL(CMISSRP) :: VALUE
  REAL(CMISSRP) :: MU_PARAM_NAVIER_STOKES
  REAL(CMISSRP) :: RHO_PARAM_NAVIER_STOKES

  REAL(CMISSRP) :: DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME
  REAL(CMISSRP) :: DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME
  REAL(CMISSRP) :: DYNAMIC_SOLVER_NAVIER_STOKES_THETA
  REAL(CMISSRP) :: DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG
  LOGICAL :: BOUNDARY_NODES_NAVIER_STOKES_FLAG

  !CMISS variables

  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsNavierStokes
  TYPE(cmfe_BasisType) :: BasisSpace,BasisVelocity,BasisPressure

  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: EquationsNavierStokes
  TYPE(cmfe_EquationsSetType) :: EquationsSetNavierStokes
  TYPE(cmfe_FieldType) :: GeometricField,DependentFieldNavierStokes,MaterialsFieldNavierStokes,AnalyticFieldNavierStokes
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_MeshElementsType) :: MeshElementsSpace,MeshElementsVelocity,MeshElementsPressure
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_SolverType) :: LinearSolverNavierStokes,NonlinearSolverNavierStokes,DynamicSolverNavierStokes
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsNavierStokes
  TYPE(cmfe_FieldType) :: EquationsSetField

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables

  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,BoundaryNodeDomain
  INTEGER(CMISSIntg) :: EquationsSetIndex
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

  !-----------------------------------------------------------------------------------------------------------
  !PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------
  !Import cmHeart mesh information
  CALL FLUID_MECHANICS_IO_READ_CMHEART(CM,Err)  
  BASIS_NUMBER_SPACE=CM%ID_M
  BASIS_NUMBER_VELOCITY=CM%ID_V
  BASIS_NUMBER_PRESSURE=CM%ID_P
  NUMBER_OF_DIMENSIONS=CM%D
  BASIS_TYPE=CM%IT_T
  BASIS_XI_INTERPOLATION_SPACE=CM%IT_M
  BASIS_XI_INTERPOLATION_VELOCITY=CM%IT_V
  BASIS_XI_INTERPOLATION_PRESSURE=CM%IT_P
  NUMBER_OF_NODES_SPACE=CM%N_M
  NUMBER_OF_NODES_VELOCITY=CM%N_V
  NUMBER_OF_NODES_PRESSURE=CM%N_P
  TOTAL_NUMBER_OF_NODES=CM%N_T
  TOTAL_NUMBER_OF_ELEMENTS=CM%E_T
  NUMBER_OF_ELEMENT_NODES_SPACE=CM%EN_M
  NUMBER_OF_ELEMENT_NODES_VELOCITY=CM%EN_V
  NUMBER_OF_ELEMENT_NODES_PRESSURE=CM%EN_P
  !Set initial values
  INITIAL_FIELD_NAVIER_STOKES(1)=0.0_CMISSRP
  INITIAL_FIELD_NAVIER_STOKES(2)=0.0_CMISSRP
  !Set boundary conditions
  BOUNDARY_NODES_NAVIER_STOKES_FLAG=.TRUE.
  IF(BOUNDARY_NODES_NAVIER_STOKES_FLAG) THEN

! ! 10x10 element unit square
!    NUMBER_OF_BOUNDARY_NODES_NAVIER_STOKES=80
!     ALLOCATE(BOUNDARY_NODES_NAVIER_STOKES(NUMBER_OF_BOUNDARY_NODES_NAVIER_STOKES))
!       BOUNDARY_NODES_NAVIER_STOKES=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &
!     & 17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32, &
!     & 33,34,35,36,37,38,39,40,122,123,124,125,126,127,128,129, &
!     & 130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145, &
!     & 146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161]

! 20x20 element unit square
    NUMBER_OF_BOUNDARY_NODES_NAVIER_STOKES=160
    ALLOCATE(BOUNDARY_NODES_NAVIER_STOKES(NUMBER_OF_BOUNDARY_NODES_NAVIER_STOKES))
      BOUNDARY_NODES_NAVIER_STOKES=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &
    & 17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32, &
    & 33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48, &
    & 49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64, &
    & 65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80, &
    & 442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457, &
    & 458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473, &
    & 474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489, &
    & 490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505, &
    & 506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521]

  ENDIF
  !Set material parameters
  MU_PARAM_NAVIER_STOKES=VISCOSITY
  RHO_PARAM_NAVIER_STOKES=DENSITY
  !Set interpolation parameters
  BASIS_XI_GAUSS_SPACE=2
  BASIS_XI_GAUSS_VELOCITY=2
  BASIS_XI_GAUSS_PRESSURE=2
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMFE_SOLVER_PROGRESS_OUTPUT
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMFE_SOLVER_TIMING_OUTPUT
  NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMFE_SOLVER_TIMING_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_NAVIER_STOKES_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT
  !Set time parameter
  DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME=START_TIME
  DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME=STOP_TIME 
  DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT=TIME_INCREMENT
  DYNAMIC_SOLVER_NAVIER_STOKES_THETA=THETA
  !Set result output parameter
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY=OUTPUT_FREQUENCY

  !Set solver parameters
  LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-5_CMISSRP !default: 1.0E-05_CMISSRP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-10_CMISSRP
  DIVERGENCE_TOLERANCE=1.0E5 !default: 1.0E5
  MAXIMUM_ITERATIONS=100000 !default: 100000
  RESTART_VALUE=3000 !default: 30
  LINESEARCH_ALPHA=1.0

  !-----------------------------------------------------------------------------------------------------------
  !COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)


  !-----------------------------------------------------------------------------------------------------------
  !REGION
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of a new region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system as defined above
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)


  !-----------------------------------------------------------------------------------------------------------
  !BASIS
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of new bases
  MESH_NUMBER_OF_COMPONENTS=1
  CALL cmfe_Basis_Initialise(BasisSpace,Err)
  CALL cmfe_Basis_CreateStart(BASIS_NUMBER_SPACE,BasisSpace,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL cmfe_Basis_TypeSet(BasisSpace,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL cmfe_Basis_NumberOfXiSet(BasisSpace,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  CALL cmfe_Basis_InterpolationXiSet(BasisSpace,[BASIS_XI_INTERPOLATION_SPACE,BASIS_XI_INTERPOLATION_SPACE],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisSpace,[BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE],Err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(BasisSpace,Err)
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisVelocity=BasisSpace
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new velocity basis
    CALL cmfe_Basis_Initialise(BasisVelocity,Err)
    !Start the creation of a basis
    CALL cmfe_Basis_CreateStart(BASIS_NUMBER_VELOCITY,BasisVelocity,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL cmfe_Basis_TypeSet(BasisVelocity,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL cmfe_Basis_NumberOfXiSet(BasisVelocity,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisVelocity,[BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisVelocity,[BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY],Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisVelocity,[BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY, & 
        & BASIS_XI_INTERPOLATION_VELOCITY],Err)                         
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisVelocity,[BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY, & 
        & BASIS_XI_GAUSS_VELOCITY],Err)
    ENDIF
    !Finish the creation of the basis
    CALL cmfe_Basis_CreateFinish(BasisVelocity,Err)
  ENDIF
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisPressure=BasisSpace
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_VELOCITY) THEN
    BasisPressure=BasisVelocity
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new pressure basis
    CALL cmfe_Basis_Initialise(BasisPressure,Err)
    !Start the creation of a basis
    CALL cmfe_Basis_CreateStart(BASIS_NUMBER_PRESSURE,BasisPressure,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL cmfe_Basis_TypeSet(BasisPressure,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL cmfe_Basis_NumberOfXiSet(BasisPressure,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisPressure,[BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisPressure,[BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE],Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisPressure,[BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE, & 
        & BASIS_XI_INTERPOLATION_PRESSURE],Err)                         
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisPressure,[BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE, & 
        & BASIS_XI_GAUSS_PRESSURE],Err)
    ENDIF
    !Finish the creation of the basis
    CALL cmfe_Basis_CreateFinish(BasisPressure,Err)
  ENDIF


  !-----------------------------------------------------------------------------------------------------------
  !MESH
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of mesh nodes
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Nodes_CreateStart(Region,TOTAL_NUMBER_OF_NODES,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)
  !Specify spatial mesh component
  CALL cmfe_MeshElements_Initialise(MeshElementsSpace,Err)
  CALL cmfe_MeshElements_Initialise(MeshElementsVelocity,Err)
  CALL cmfe_MeshElements_Initialise(MeshElementsPressure,Err)
  MESH_COMPONENT_NUMBER_SPACE=1
  MESH_COMPONENT_NUMBER_VELOCITY=1
  MESH_COMPONENT_NUMBER_PRESSURE=1
  CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_SPACE,BasisSpace,MeshElementsSpace,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL cmfe_MeshElements_NodesSet(MeshElementsSpace,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_SPACE),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(MeshElementsSpace,Err)
  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsVelocity=MeshElementsSpace
  ELSE
    MESH_COMPONENT_NUMBER_VELOCITY=MESH_COMPONENT_NUMBER_SPACE+1
    CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_VELOCITY,BasisVelocity,MeshElementsVelocity,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL cmfe_MeshElements_NodesSet(MeshElementsVelocity,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_VELOCITY),Err)
    ENDDO
    CALL cmfe_MeshElements_CreateFinish(MeshElementsVelocity,Err)
  ENDIF
  !Specify pressure mesh component
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsPressure=MeshElementsSpace
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_SPACE
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_VELOCITY) THEN
    MeshElementsPressure=MeshElementsVelocity
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_VELOCITY
  ELSE
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_VELOCITY+1
    CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_PRESSURE,BasisPressure,MeshElementsPressure,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL cmfe_MeshElements_NodesSet(MeshElementsPressure,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_PRESSURE),Err)
    ENDDO
    CALL cmfe_MeshElements_CreateFinish(MeshElementsPressure,Err)
  ENDIF
  !Finish the creation of the mesh
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)


  !-----------------------------------------------------------------------------------------------------------  
  !GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field type
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL cmfe_Field_ScalingTypeSet(GeometricField,CMFE_FIELD_NO_SCALING,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_SPACE
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
          & 1,CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
      ENDIF
    ENDDO
  ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)


  !-----------------------------------------------------------------------------------------------------------  
  !EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set for Navier-Stokes
  CALL cmfe_EquationsSet_Initialise(EquationsSetNavierStokes,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField, &
    & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMFE_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSetNavierStokes,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetNavierStokes,Err)


  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set dependent field variables for dynamic Navier-Stokes
  CALL cmfe_Field_Initialise(DependentFieldNavierStokes,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetNavierStokes,DependentFieldUserNumberNavierStokes, & 
    & DependentFieldNavierStokes,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
  ENDDO
  COMPONENT_NUMBER=NUMBER_OF_DIMENSIONS+1
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetNavierStokes,Err)

  !Initialise dependent field
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & COMPONENT_NUMBER,INITIAL_FIELD_NAVIER_STOKES(COMPONENT_NUMBER),Err)
  ENDDO

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIALS FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set materials field variables for dynamic Navier-Stokes
  CALL cmfe_Field_Initialise(MaterialsFieldNavierStokes,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetNavierStokes,MaterialsFieldUserNumberNavierStokes, & 
    & MaterialsFieldNavierStokes,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetNavierStokes,Err)
!  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
!    & MaterialsFieldUserNumberNavierStokesMu,MU_PARAM_NAVIER_STOKES,Err)
!  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
!    & MaterialsFieldUserNumberNavierStokesRho,RHO_PARAM_NAVIER_STOKES,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 1,VISCOSITY,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 2,DENSITY,Err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set analytic field variables
  CALL cmfe_Field_Initialise(AnalyticFieldNavierStokes,Err)
  CALL cmfe_EquationsSet_AnalyticCreateStart(EquationsSetNavierStokes, &
    & CMFE_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN, &
    & AnalyticFieldUserNumber,AnalyticFieldNavierStokes,Err)
  !Finish the equations set analytic field variables
  CALL cmfe_EquationsSet_AnalyticCreateFinish(EquationsSetNavierStokes,Err)
  !Set analytic field variables
  CALL cmfe_Field_ComponentValuesInitialise(AnalyticFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 1,U_0,Err)
  CALL cmfe_Field_ComponentValuesInitialise(AnalyticFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 2,LENGTH,Err)


  !-----------------------------------------------------------------------------------------------------------  
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsNavierStokes,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetNavierStokes,EquationsNavierStokes,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(EquationsNavierStokes,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsNavierStokes,EQUATIONS_NAVIER_STOKES_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetNavierStokes,Err)


  !-----------------------------------------------------------------------------------------------------------
  !PROBLEM
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  !Set the problem to be a transient Navier-Stokes problem
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS,CMFE_PROBLEM_NAVIER_STOKES_EQUATION_TYPE, &
    CMFE_PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME,DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME, & 
    & DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT,Err)
  !Set the output timing
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)


  !-----------------------------------------------------------------------------------------------------------
  !SOLVER
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(DynamicSolverNavierStokes,Err)
  CALL cmfe_Solver_Initialise(NonlinearSolverNavierStokes,Err)
  CALL cmfe_Solver_Initialise(LinearSolverNavierStokes,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)

  !Get the dynamic dymamic solver
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(DynamicSolverNavierStokes,DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set theta
  CALL cmfe_Solver_DynamicThetaSet(DynamicSolverNavierStokes,DYNAMIC_SOLVER_NAVIER_STOKES_THETA,Err)
  !Get the dynamic nonlinear solver
  CALL cmfe_Solver_DynamicNonlinearSolverGet(DynamicSolverNavierStokes,NonlinearSolverNavierStokes,Err)
  !Set the nonlinear Jacobian type
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes, &
    & CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(NonlinearSolverNavierStokes,NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(NonlinearSolverNavierStokes,ABSOLUTE_TOLERANCE,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(NonlinearSolverNavierStokes,RELATIVE_TOLERANCE,Err)
  !Get the dynamic nonlinear linear solver
  CALL cmfe_Solver_NewtonLinearSolverGet(NonlinearSolverNavierStokes,LinearSolverNavierStokes,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(LinearSolverNavierStokes,LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolverNavierStokes,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LibraryTypeSet(LinearSolverNavierStokes,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolverNavierStokes,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolverNavierStokes,MAXIMUM_ITERATIONS,Err)
    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolverNavierStokes,DIVERGENCE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolverNavierStokes,RELATIVE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolverNavierStokes,ABSOLUTE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(LinearSolverNavierStokes,RESTART_VALUE,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)


  !-----------------------------------------------------------------------------------------------------------
  !SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(DynamicSolverNavierStokes,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsNavierStokes,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the dynamic solver equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  CALL cmfe_Solver_SolverEquationsGet(DynamicSolverNavierStokes,SolverEquationsNavierStokes,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsNavierStokes,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsNavierStokes,EquationsSetNavierStokes,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)


  !-----------------------------------------------------------------------------------------------------------
  !BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of the equations set boundary conditions for Stokes
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsNavierStokes,Err)
!  CALL cmfe_EquationsSet_BoundaryConditionsCreateStart(EquationsSetNavierStokes,BoundaryConditionsNavierStokes,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsNavierStokes,BoundaryConditionsNavierStokes,Err)
  !Identify boundary nodes- BoundaryConditionFixed nodes will be set with analytic boundary conditions.
  IF(BOUNDARY_NODES_NAVIER_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_BOUNDARY_NODES_NAVIER_STOKES
      NODE_NUMBER=BOUNDARY_NODES_NAVIER_STOKES(NODE_COUNTER)
      CONDITION=CMFE_BOUNDARY_CONDITION_FIXED
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
          !Initialize boundary nodes to 0
          VALUE=0.0_CMISSRP
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
            & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsNavierStokes,Err)

  !-----------------------------------------------------------------------------------------------------------
  !SOLVE
  !-----------------------------------------------------------------------------------------------------------

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL cmfe_Problem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"


  !-----------------------------------------------------------------------------------------------------------
  !OUTPUT
  !-----------------------------------------------------------------------------------------------------------

  !Output Analytic analysis
  CALL cmfe_AnalyticAnalysis_Output(DependentFieldNavierStokes,"Navier_Stokes_Analytic_Taylor_Green_2D",Err)
  !Finialise CMISS
  CALL cmfe_Finalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

 END PROGRAM TAYLORGREENEXAMPLE
