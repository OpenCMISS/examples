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
!! Example program to solve a dynamic Taylor Green vortex using openCMISS calls.
!<

!> Main program
PROGRAM TAYLORGREENEXAMPLE

  ! PROGRAM LIBRARIES
  USE OPENCMISS
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


  REAL(CMISSDP), PARAMETER :: U_0=1.0_CMISSDP !Characteristic velocity (velocity boundary amplitude at t=0)
  REAL(CMISSDP), PARAMETER :: LENGTH=1.0_CMISSDP !Length scale (length of one side of the square)
  REAL(CMISSDP), PARAMETER :: VISCOSITY=0.01_CMISSDP
  REAL(CMISSDP), PARAMETER :: DENSITY=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: START_TIME=0.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: STOP_TIME=1.0001_CMISSDP
  REAL(CMISSDP), PARAMETER :: TIME_INCREMENT=0.01_CMISSDP
  REAL(CMISSDP), PARAMETER :: THETA=1.0_CMISSDP
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
  REAL(CMISSDP) :: POSITION,TOTAL_LENGTH
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

  REAL(CMISSDP) :: INITIAL_FIELD_NAVIER_STOKES(2)
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_NAVIER_STOKES(2)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE
  REAL(CMISSDP) :: MU_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: RHO_PARAM_NAVIER_STOKES

  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_THETA
  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG
  LOGICAL :: BOUNDARY_NODES_NAVIER_STOKES_FLAG

  !CMISS variables

  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsNavierStokes
  TYPE(CMISSBasisType) :: BasisSpace,BasisVelocity,BasisPressure

  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: EquationsNavierStokes
  TYPE(CMISSEquationsSetType) :: EquationsSetNavierStokes
  TYPE(CMISSFieldType) :: GeometricField,DependentFieldNavierStokes,MaterialsFieldNavierStokes,AnalyticFieldNavierStokes
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSMeshElementsType) :: MeshElementsSpace,MeshElementsVelocity,MeshElementsPressure
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSSolverType) :: LinearSolverNavierStokes,NonlinearSolverNavierStokes,DynamicSolverNavierStokes
  TYPE(CMISSSolverEquationsType) :: SolverEquationsNavierStokes
  TYPE(CMISSFieldType) :: EquationsSetField

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
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

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
  INITIAL_FIELD_NAVIER_STOKES(1)=0.0_CMISSDP
  INITIAL_FIELD_NAVIER_STOKES(2)=0.0_CMISSDP
  !Set boundary conditions
  BOUNDARY_NODES_NAVIER_STOKES_FLAG=.TRUE.
  IF(BOUNDARY_NODES_NAVIER_STOKES_FLAG) THEN

! ! 10x10 element unit square
!    NUMBER_OF_BOUNDARY_NODES_NAVIER_STOKES=80
!     ALLOCATE(BOUNDARY_NODES_NAVIER_STOKES(NUMBER_OF_BOUNDARY_NODES_NAVIER_STOKES))
!       BOUNDARY_NODES_NAVIER_STOKES=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &
!     & 17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32, &
!     & 33,34,35,36,37,38,39,40,122,123,124,125,126,127,128,129, &
!     & 130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145, &
!     & 146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161/)

! 20x20 element unit square
    NUMBER_OF_BOUNDARY_NODES_NAVIER_STOKES=160
    ALLOCATE(BOUNDARY_NODES_NAVIER_STOKES(NUMBER_OF_BOUNDARY_NODES_NAVIER_STOKES))
      BOUNDARY_NODES_NAVIER_STOKES=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &
    & 17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32, &
    & 33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48, &
    & 49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64, &
    & 65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80, &
    & 442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457, &
    & 458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473, &
    & 474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489, &
    & 490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505, &
    & 506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521/)

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
  LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISSSolverProgressOutput
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISSSolverTimingOutput
  NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISSSolverTimingOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_NAVIER_STOKES_OUTPUT=CMISSEquationsNoOutput
  !Set time parameter
  DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME=START_TIME
  DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME=STOP_TIME 
  DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT=TIME_INCREMENT
  DYNAMIC_SOLVER_NAVIER_STOKES_THETA=THETA
  !Set result output parameter
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY=OUTPUT_FREQUENCY

  !Set solver parameters
  LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-5_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E5 !default: 1.0E5
  MAXIMUM_ITERATIONS=100000 !default: 100000
  RESTART_VALUE=3000 !default: 30
  LINESEARCH_ALPHA=1.0

  !-----------------------------------------------------------------------------------------------------------
  !COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)


  !-----------------------------------------------------------------------------------------------------------
  !REGION
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of a new region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system as defined above
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)


  !-----------------------------------------------------------------------------------------------------------
  !BASIS
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of new bases
  MESH_NUMBER_OF_COMPONENTS=1
  CALL CMISSBasisTypeInitialise(BasisSpace,Err)
  CALL CMISSBasisCreateStart(BASIS_NUMBER_SPACE,BasisSpace,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasisTypeSet(BasisSpace,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL CMISSBasisNumberOfXiSet(BasisSpace,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  CALL CMISSBasisInterpolationXiSet(BasisSpace,(/BASIS_XI_INTERPOLATION_SPACE,BASIS_XI_INTERPOLATION_SPACE/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisSpace,(/BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE/),Err)
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(BasisSpace,Err)
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisVelocity=BasisSpace
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new velocity basis
    CALL CMISSBasisTypeInitialise(BasisVelocity,Err)
    !Start the creation of a basis
    CALL CMISSBasisCreateStart(BASIS_NUMBER_VELOCITY,BasisVelocity,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasisTypeSet(BasisVelocity,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL CMISSBasisNumberOfXiSet(BasisVelocity,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL CMISSBasisInterpolationXiSet(BasisVelocity,(/BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY/),Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisVelocity,(/BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL CMISSBasisInterpolationXiSet(BasisVelocity,(/BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY, & 
        & BASIS_XI_INTERPOLATION_VELOCITY/),Err)                         
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisVelocity,(/BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY, & 
        & BASIS_XI_GAUSS_VELOCITY/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(BasisVelocity,Err)
  ENDIF
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisPressure=BasisSpace
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_VELOCITY) THEN
    BasisPressure=BasisVelocity
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new pressure basis
    CALL CMISSBasisTypeInitialise(BasisPressure,Err)
    !Start the creation of a basis
    CALL CMISSBasisCreateStart(BASIS_NUMBER_PRESSURE,BasisPressure,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasisTypeSet(BasisPressure,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL CMISSBasisNumberOfXiSet(BasisPressure,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL CMISSBasisInterpolationXiSet(BasisPressure,(/BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE/),Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisPressure,(/BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL CMISSBasisInterpolationXiSet(BasisPressure,(/BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE, & 
        & BASIS_XI_INTERPOLATION_PRESSURE/),Err)                         
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisPressure,(/BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE, & 
        & BASIS_XI_GAUSS_PRESSURE/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(BasisPressure,Err)
  ENDIF


  !-----------------------------------------------------------------------------------------------------------
  !MESH
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of mesh nodes
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSNodesCreateStart(Region,TOTAL_NUMBER_OF_NODES,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMeshNumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMeshNumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)
  !Specify spatial mesh component
  CALL CMISSMeshElementsTypeInitialise(MeshElementsSpace,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsVelocity,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsPressure,Err)
  MESH_COMPONENT_NUMBER_SPACE=1
  MESH_COMPONENT_NUMBER_VELOCITY=1
  MESH_COMPONENT_NUMBER_PRESSURE=1
  CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_SPACE,BasisSpace,MeshElementsSpace,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL CMISSMeshElementsNodesSet(MeshElementsSpace,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_SPACE),Err)
  ENDDO
  CALL CMISSMeshElementsCreateFinish(MeshElementsSpace,Err)
  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsVelocity=MeshElementsSpace
  ELSE
    MESH_COMPONENT_NUMBER_VELOCITY=MESH_COMPONENT_NUMBER_SPACE+1
    CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_VELOCITY,BasisVelocity,MeshElementsVelocity,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElementsNodesSet(MeshElementsVelocity,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_VELOCITY),Err)
    ENDDO
    CALL CMISSMeshElementsCreateFinish(MeshElementsVelocity,Err)
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
    CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_PRESSURE,BasisPressure,MeshElementsPressure,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElementsNodesSet(MeshElementsPressure,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_PRESSURE),Err)
    ENDDO
    CALL CMISSMeshElementsCreateFinish(MeshElementsPressure,Err)
  ENDIF
  !Finish the creation of the mesh
  CALL CMISSMeshCreateFinish(Mesh,Err)


  !-----------------------------------------------------------------------------------------------------------  
  !GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field type
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL CMISSFieldScalingTypeSet(GeometricField,CMISSFieldNoScaling,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_SPACE
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
          & 1,CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
      ENDIF
    ENDDO
  ENDDO
  CALL CMISSFieldParameterSetUpdateStart(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)


  !-----------------------------------------------------------------------------------------------------------  
  !EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set for Navier-Stokes
  CALL CMISSEquationsSetTypeInitialise(EquationsSetNavierStokes,Err)
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField, &
    & CMISSEquationsSetFluidMechanicsClass,CMISSEquationsSetNavierStokesEquationType, &
    & CMISSEquationsSetTransientNavierStokesSubtype,EquationsSetFieldUserNumber,EquationsSetField,EquationsSetNavierStokes,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSetNavierStokes,Err)


  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set dependent field variables for dynamic Navier-Stokes
  CALL CMISSFieldTypeInitialise(DependentFieldNavierStokes,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetNavierStokes,DependentFieldUserNumberNavierStokes, & 
    & DependentFieldNavierStokes,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
  ENDDO
  COMPONENT_NUMBER=NUMBER_OF_DIMENSIONS+1
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetNavierStokes,Err)

  !Initialise dependent field
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSFieldComponentValuesInitialise(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
      & COMPONENT_NUMBER,INITIAL_FIELD_NAVIER_STOKES(COMPONENT_NUMBER),Err)
  ENDDO

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIALS FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set materials field variables for dynamic Navier-Stokes
  CALL CMISSFieldTypeInitialise(MaterialsFieldNavierStokes,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetNavierStokes,MaterialsFieldUserNumberNavierStokes, & 
    & MaterialsFieldNavierStokes,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetNavierStokes,Err)
!  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
!    & MaterialsFieldUserNumberNavierStokesMu,MU_PARAM_NAVIER_STOKES,Err)
!  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
!    & MaterialsFieldUserNumberNavierStokesRho,RHO_PARAM_NAVIER_STOKES,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,VISCOSITY,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 2,DENSITY,Err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set analytic field variables
  CALL CMISSFieldTypeInitialise(AnalyticFieldNavierStokes,Err)
  CALL CMISSEquationsSetAnalyticCreateStart(EquationsSetNavierStokes,CMISSEquationsSetNavierStokesTwoDimTaylorGreen, &
    & AnalyticFieldUserNumber,AnalyticFieldNavierStokes,Err)
  !Finish the equations set analytic field variables
  CALL CMISSEquationsSetAnalyticCreateFinish(EquationsSetNavierStokes,Err)
  !Set analytic field variables
  CALL CMISSFieldComponentValuesInitialise(AnalyticFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,U_0,Err)
  CALL CMISSFieldComponentValuesInitialise(AnalyticFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 2,LENGTH,Err)


  !-----------------------------------------------------------------------------------------------------------  
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------
  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsNavierStokes,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetNavierStokes,EquationsNavierStokes,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(EquationsNavierStokes,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(EquationsNavierStokes,EQUATIONS_NAVIER_STOKES_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetNavierStokes,Err)


  !-----------------------------------------------------------------------------------------------------------
  !PROBLEM
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a transient Navier-Stokes problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemFluidMechanicsClass,CMISSProblemNavierStokesEquationType, &
    & CMISSProblemTransientNavierStokesSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME,DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME, & 
    & DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT,Err)
  !Set the output timing
  CALL CMISSControlLoopTimeOutputSet(ControlLoop,DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)


  !-----------------------------------------------------------------------------------------------------------
  !SOLVER
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(DynamicSolverNavierStokes,Err)
  CALL CMISSSolverTypeInitialise(NonlinearSolverNavierStokes,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverNavierStokes,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)

  !Get the dynamic dymamic solver
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(DynamicSolverNavierStokes,DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set theta
  CALL CMISSSolverDynamicThetaSet(DynamicSolverNavierStokes,DYNAMIC_SOLVER_NAVIER_STOKES_THETA,Err)
!   CALL CMISSSolverDynamicDynamicSet(DynamicSolverNavierStokes,.TRUE.,Err)
  !Get the dynamic nonlinear solver
  CALL CMISSSolverDynamicNonlinearSolverGet(DynamicSolverNavierStokes,NonlinearSolverNavierStokes,Err)
  !Set the nonlinear Jacobian type
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes,CMISSSolverNewtonJacobianAnalyticCalculated,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(NonlinearSolverNavierStokes,NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  CALL CMISSSolverNewtonAbsoluteToleranceSet(NonlinearSolverNavierStokes,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolverNewtonRelativeToleranceSet(NonlinearSolverNavierStokes,RELATIVE_TOLERANCE,Err)
  !Get the dynamic nonlinear linear solver
  CALL CMISSSolverNewtonLinearSolverGet(NonlinearSolverNavierStokes,LinearSolverNavierStokes,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(LinearSolverNavierStokes,LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG) THEN
    CALL CMISSSolverLinearTypeSet(LinearSolverNavierStokes,CMISSSolverLinearDirectSolveType,Err)
    CALL CMISSSolverLibraryTypeSet(LinearSolverNavierStokes,CMISSSolverMUMPSLibrary,Err)
  ELSE
    CALL CMISSSolverLinearTypeSet(LinearSolverNavierStokes,CMISSSolverLinearIterativeSolveType,Err)
    CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolverNavierStokes,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolverLinearIterativeDivergenceToleranceSet(LinearSolverNavierStokes,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeRelativeToleranceSet(LinearSolverNavierStokes,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(LinearSolverNavierStokes,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeGMRESRestartSet(LinearSolverNavierStokes,RESTART_VALUE,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)


  !-----------------------------------------------------------------------------------------------------------
  !SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solver equations
  CALL CMISSSolverTypeInitialise(DynamicSolverNavierStokes,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsNavierStokes,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the dynamic solver equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  CALL CMISSSolverSolverEquationsGet(DynamicSolverNavierStokes,SolverEquationsNavierStokes,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsNavierStokes,CMISSSolverEquationsSparseMatrices,Err)
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsNavierStokes,EquationsSetNavierStokes,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)


  !-----------------------------------------------------------------------------------------------------------
  !BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------
  !Start the creation of the equations set boundary conditions for Stokes
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsNavierStokes,Err)
!  CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSetNavierStokes,BoundaryConditionsNavierStokes,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsNavierStokes,BoundaryConditionsNavierStokes,Err)
  !Identify boundary nodes- BoundaryConditionFixed nodes will be set with analytic boundary conditions.
  IF(BOUNDARY_NODES_NAVIER_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_BOUNDARY_NODES_NAVIER_STOKES
      NODE_NUMBER=BOUNDARY_NODES_NAVIER_STOKES(NODE_COUNTER)
      CONDITION=CMISSBoundaryConditionFixed
      CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
          !Initialize boundary nodes to 0
          VALUE=0.0_CMISSDP
          CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes,CMISSFieldUVariableType,1, &
            & CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsNavierStokes,Err)

  !-----------------------------------------------------------------------------------------------------------
  !SOLVE
  !-----------------------------------------------------------------------------------------------------------

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL CMISSProblemSolve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"


  !-----------------------------------------------------------------------------------------------------------
  !OUTPUT
  !-----------------------------------------------------------------------------------------------------------

  !Output Analytic analysis
  CALL CMISSAnalyticAnalysisOutput(DependentFieldNavierStokes,"Navier_Stokes_Analytic_Taylor_Green_2D",Err)
  !Finialise CMISS
  CALL CMISSFinalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

 END PROGRAM TAYLORGREENEXAMPLE
