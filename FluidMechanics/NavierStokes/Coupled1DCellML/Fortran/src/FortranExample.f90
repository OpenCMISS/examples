!> \file
!> \author David Ladd
!> \brief This is an example program to solve 1D Transient Navier-Stokes over a bifurcation
!>  with coupled 0D lumped models (resistance, RCR) defined in CellML.
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
!> \example FluidMechanics/NavierStokes/Coupled1DCellML/src/Coupled1DCellMLExample.f90
!!
!<

!> Main program
PROGRAM FortranExample

  ! Libraries

  USE OPENCMISS
  USE FLUID_MECHANICS_IO_ROUTINES
  USE MPI

#ifdef WIN32
  USE IFQWINCMISS
#endif


  IMPLICIT NONE

  !================================================================================================================================
  !Test program parameters
  !================================================================================================================================

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberNavierStokes=1337
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberCharacteristic=1338
  TYPE(CMISSFieldType) :: EquationsSetFieldNavierStokes
  TYPE(CMISSFieldType) :: EquationsSetFieldCharacteristic

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberNavierStokes=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberCharacteristic=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: Basis_Lagrange_Hermite_TP_Type=1
  INTEGER(CMISSIntg), PARAMETER :: Basis_Quadratic_Lagrange_Interpolation=2

  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberCellML=18

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverCharacteristicUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverNavierStokesUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMu=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberRho=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberK=3
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberBs=4
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberAs=5
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberRe=6
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberFr=7
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberSt=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberA0=9
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberBeta=10
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberE=11
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberH0=12


  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_SPACE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_AREA
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_AREA
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_AREA
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_SPACE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_AREA
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_AREA
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_AREA
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE

  INTEGER(CMISSIntg) :: EQUATIONS_NAVIER_STOKES_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: CONDITION,i

  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE

  INTEGER, ALLOCATABLE, DIMENSION(:):: OUTLET_WALL_NODES_NAVIER_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_NAVIER_STOKES

  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_NAVIER_STOKES(2)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE1
  REAL(CMISSDP) :: MU_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: RHO_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: E_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: H0_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: A0_PARAM(6),Beta(6),Q1,A1,Ts,As,Qs,Xs,St,Fr,Re,Bs,K
  REAL(CMISSDP) :: Q2,Q3,A2,A3
  REAL(CMISSDP) :: resistanceProximal,resistanceDistal,capacitance
  REAL(CMISSDP) :: pCellML,pPrevious,pVesselWall,pExternal,qPrevious
  REAL(CMISSDP) :: length, position
  REAL(CMISSDP) :: VALUE,X,Y,Z

  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_THETA
  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT

  INTEGER(CMISSIntg) :: pCellMLComponent,pPreviousComponent,pVesselWallComponent,pExternalComponent,qPreviousComponent
  INTEGER(CMISSIntg) :: EquationsSetSubtype
  INTEGER(CMISSIntg) :: EquationsSetCharacteristicSubtype
  INTEGER(CMISSIntg) :: ProblemSubtype
  INTEGER(CMISSIntg) :: coupledNodeNumber,coupledNodeNumber1,coupledNodeNumber2
  INTEGER(CMISSIntg) :: resistanceComponent,numberOfCoordinateDimensions
  INTEGER(CMISSIntg) :: nodeIdx,versionIdx,componentIdx,normalWave

  LOGICAL :: LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG
  LOGICAL :: OUTLET_WALL_NODES_NAVIER_STOKES_FLAG
  LOGICAL :: INLET_WALL_NODES_NAVIER_STOKES_FLAG
  LOGICAL :: cellmlFlag,versionsFlag,windkesselFlag

  !CMISS variables

  !Regions
  TYPE(CMISSRegionType) :: Region
  TYPE(CMISSRegionType) :: WorldRegion
  !Coordinate systems
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  !Basis
  TYPE(CMISSBasisType) :: BasisSpace
  TYPE(CMISSBasisType) :: BasisVelocity
  TYPE(CMISSBasisType) :: BasisArea
  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Elements
  TYPE(CMISSMeshElementsType) :: MeshElementsSpace
  TYPE(CMISSMeshElementsType) :: MeshElementsVelocity
  TYPE(CMISSMeshElementsType) :: MeshElementsArea
  !Meshes
  TYPE(CMISSMeshType) :: Mesh
  !Decompositions
  TYPE(CMISSDecompositionType) :: Decomposition
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: DependentFieldNavierStokes
  TYPE(CMISSFieldType) :: MaterialsFieldNavierStokes
  TYPE(CMISSFieldType) :: IndependentFieldNavierStokes
  ! CellML types
  TYPE(CMISSCellMLType) :: CellML
  TYPE(CMISSFieldType) :: MaterialsFieldCellML
  TYPE(CMISSCellMLEquationsType) :: CellMLEquations
  TYPE(CMISSFieldType) :: CellMLModelsField
  TYPE(CMISSFieldType) :: CellMLStateField
  TYPE(CMISSFieldType) :: CellMLIntermediateField
  TYPE(CMISSFieldType) :: CellMLParametersField
  TYPE(CMISSSolverType) :: CellMLSolver
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsNavierStokes
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsCharacteristic
  !Equations sets
  TYPE(CMISSEquationsSetType) :: EquationsSetNavierStokes
  TYPE(CMISSEquationsSetType) :: EquationsSetCharacteristic
  !Equations
  TYPE(CMISSEquationsType) :: EquationsNavierStokes
  TYPE(CMISSEquationsType) :: EquationsCharacteristic
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: DynamicSolverNavierStokes
  TYPE(CMISSSolverType) :: NonlinearSolverNavierStokes
  TYPE(CMISSSolverType) :: LinearSolverNavierStokes
  TYPE(CMISSSolverType) :: NonlinearSolverCharacteristic
  TYPE(CMISSSolverType) :: LinearSolverCharacteristic
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsNavierStokes
  TYPE(CMISSSolverEquationsType) :: SolverEquationsCharacteristic

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,BoundaryNodeDomain
  INTEGER(CMISSIntg) :: CellMLIndex
  INTEGER(CMISSIntg) :: ResistanceModelIndex
  INTEGER(CMISSIntg) :: WindkesselModelIndex
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
  

  !================================================================================================================================
  !Intialise OpenCMISS
  !================================================================================================================================

  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)
!  CALL CMISSDiagnosticsSetOn(CMISS_ALL_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",[""],Err)
!  CALL CMISSOutputSetOn("Testing",Err)
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !================================================================================================================================
  !PROBLEM CONTROL PANEL
  !================================================================================================================================

  cellmlFlag = .FALSE.
  windkesselFlag = .FALSE.
  versionsFlag = .FALSE.
  numberOfCoordinateDimensions=2
!  resistanceProximal=9.2119E+11_CMISSDP
!  resistanceProximal=9.125E+10_CMISSDP
!  resistanceProximal=1.3201E+11_CMISSDP
!   resistanceProximal=0.0_CMISSDP

  resistanceProximal=1.7025E+7_CMISSDP
!  resistanceProximal=1.0_CMISSDP
  resistanceDistal=0.0_CMISSDP
  capacitance=0.0_CMISSDP
  BASIS_NUMBER_SPACE=1
  BASIS_NUMBER_VELOCITY=2
  BASIS_NUMBER_AREA=3
  NUMBER_OF_DIMENSIONS=1
  BASIS_TYPE=1
  BASIS_XI_INTERPOLATION_SPACE=2
  BASIS_XI_INTERPOLATION_VELOCITY=2
  BASIS_XI_INTERPOLATION_AREA=2
  NUMBER_OF_NODES_SPACE=13
  NUMBER_OF_NODES_VELOCITY=NUMBER_OF_NODES_SPACE
  NUMBER_OF_NODES_AREA=NUMBER_OF_NODES_SPACE
  TOTAL_NUMBER_OF_NODES=NUMBER_OF_NODES_SPACE*3
  TOTAL_NUMBER_OF_ELEMENTS=6
  NUMBER_OF_ELEMENT_NODES_SPACE=TOTAL_NUMBER_OF_ELEMENTS
  NUMBER_OF_ELEMENT_NODES_VELOCITY=TOTAL_NUMBER_OF_ELEMENTS
  NUMBER_OF_ELEMENT_NODES_AREA=TOTAL_NUMBER_OF_ELEMENTS

  !Set material parameters
  MU_PARAM_NAVIER_STOKES=0.0035_CMISSDP !Mu(Pa.s)
!  MU_PARAM_NAVIER_STOKES=1.0_CMISSDP !Mu(Pa.s)
  RHO_PARAM_NAVIER_STOKES=1050.0_CMISSDP !Rho(kg/m3)
!  RHO_PARAM_NAVIER_STOKES=1.0_CMISSDP !Rho(kg/m3)

  E_PARAM_NAVIER_STOKES=0.8E+6_CMISSDP !Elasticity(Pa)
!  E_PARAM_NAVIER_STOKES=0.8E+10_CMISSDP !Elasticity(Pa)

  H0_PARAM_NAVIER_STOKES=0.5E-3_CMISSDP !Wall Thickness(m)

  ! DO i=1,TOTAL_NUMBER_OF_ELEMENTS
  !   A0_PARAM(i)=19.6e-6_CMISSDP  !Wall Area(m2)
  ! ENDDO

  DO i=1,2
    A0_PARAM(i)=19.6e-6_CMISSDP  !Wall Area(m2)
  ENDDO
  DO i=3,TOTAL_NUMBER_OF_ELEMENTS
    A0_PARAM(i)=12.8e-6_CMISSDP
  ENDDO

  !Reference Values
  Qs=10.0e-6_CMISSDP !(m3/s)
  As=19.6e-6_CMISSDP !(m2)
  Xs=0.1_CMISSDP     !(m)
  Ts=0.1_CMISSDP     !(s)
  K=4.0/3.0_CMISSDP  !Parabolic Flow Section
!  Bs=(4.0*1.7725*E_PARAM_NAVIER_STOKES*H0_PARAM_NAVIER_STOKES)/(3.0*As)
  Bs=((4.0_CMISSDP)*(3.1416_CMISSDP**(0.5_CMISSDP))*E_PARAM_NAVIER_STOKES*H0_PARAM_NAVIER_STOKES)/((3.0_CMISSDP)*As)
  St=(As*Xs)/(Ts*Qs)
!  Fr=((As**2.5)/(Qs**2.0))*(Bs/(2.0*RHO_PARAM_NAVIER_STOKES))
  Fr=((As**(2.5_CMISSDP))/(Qs**2))*(Bs/((2.0_CMISSDP)*RHO_PARAM_NAVIER_STOKES))
  Re=8.0*3.1416*(MU_PARAM_NAVIER_STOKES*Xs)/(Qs*RHO_PARAM_NAVIER_STOKES)
  DO i=1,TOTAL_NUMBER_OF_ELEMENTS
    Beta(i)=(4.0*1.7725*E_PARAM_NAVIER_STOKES*H0_PARAM_NAVIER_STOKES)/(3.0*A0_PARAM(i))
  ENDDO

  !Set initial condition
  Q1=7.0_CMISSDP
  Q2=3.50_CMISSDP
  Q3=3.50_CMISSDP

  ! Q1=7.0_CMISSDP
  ! Q2=3.0_CMISSDP
  ! Q3=4.0_CMISSDP

  ! Q1=4.0_CMISSDP
  ! Q2=2.0_CMISSDP
  ! Q3=2.0_CMISSDP

  ! A1=0.621_CMISSDP
  ! A2=0.385_CMISSDP
  ! A3=0.389_CMISSDP

  A1=1.0_CMISSDP
  A2=0.653_CMISSDP
  A3=0.653_CMISSDP

  ! A1=0.5_CMISSDP
  ! A2=0.3_CMISSDP
  ! A3=0.2_CMISSDP

  ! A1=1.0_CMISSDP
  ! A2=0.48115_CMISSDP
  ! A3=0.249883_CMISSDP

  OUTLET_WALL_NODES_NAVIER_STOKES_FLAG=.TRUE.
  INLET_WALL_NODES_NAVIER_STOKES_FLAG=.TRUE.

  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    INLET_WALL_NODES_NAVIER_STOKES=[1]
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_NAVIER_STOKES(1)=7.0_CMISSDP
  ENDIF

  IF(OUTLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    OUTLET_WALL_NODES_NAVIER_STOKES=[NUMBER_OF_NODES_SPACE]
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_NAVIER_STOKES(2)=A1
  ENDIF

  !Set last node as the terminal node, with either versions set so that method of characteristics may be used for 1D-0D coupling
  ! or just a regular terminal boundary condition.
  coupledNodeNumber1=11
  coupledNodeNumber2=13

  !Set interpolation parameters
  BASIS_XI_GAUSS_SPACE=3
  BASIS_XI_GAUSS_VELOCITY=3
  BASIS_XI_GAUSS_AREA=3
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISS_SOLVER_PROGRESS_OUTPUT
  NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISS_SOLVER_SOLVER_OUTPUT
  LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISS_SOLVER_PROGRESS_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_NAVIER_STOKES_OUTPUT=CMISS_EQUATIONS_NO_OUTPUT
  !Set time parameter
  DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME=0.0_CMISSDP
!  DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME=10.00001_CMISSDP
  DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME=0.5_CMISSDP
  DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT=0.01_CMISSDP
  DYNAMIC_SOLVER_NAVIER_STOKES_THETA=1.0_CMISSDP/2.0_CMISSDP
  !Set result output parameter
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E20 !default: 1.0E5
  MAXIMUM_ITERATIONS=100000 !default: 100000
  RESTART_VALUE=3000 !default: 30
  LINESEARCH_ALPHA=1.0

  IF(cellmlFlag) THEN
    ! New equations set type to store p values in the Equations Set Field
    EquationsSetSubtype=CMISS_EQUATIONS_SET_Coupled1D0D_NAVIER_STOKES_SUBTYPE 
    ! Characteristic (nodal/characteristic) solver remains the same
    EquationsSetCharacteristicSubtype=CMISS_EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE
    ! New problem type to execute the 1D-0D coupling subloop at each timestep
    ProblemSubtype=CMISS_PROBLEM_Coupled1D0D_NAVIER_STOKES_SUBTYPE
  ELSE
    EquationsSetSubtype=CMISS_EQUATIONS_SET_1DTRANSIENT_NAVIER_STOKES_SUBTYPE
    EquationsSetCharacteristicSubtype=CMISS_EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE
    ProblemSubtype=CMISS_PROBLEM_1DTRANSIENT_NAVIER_STOKES_SUBTYPE
  ENDIF

  !================================================================================================================================
  !Coordinate System
  !================================================================================================================================

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,numberOfCoordinateDimensions,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !================================================================================================================================
  !Region
  !================================================================================================================================

  !Start the creation of a new region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the region label
  CALL CMISSRegion_LabelSet(Region,'OpenCMISS',Err)
  !Set the regions coordinate system as defined above
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)


  !================================================================================================================================
  !Bases
  !================================================================================================================================

  !Start the creation of SPACE bases
  MESH_NUMBER_OF_COMPONENTS=1
  CALL CMISSBasis_Initialise(BasisSpace,Err)
  CALL CMISSBasis_CreateStart(BASIS_NUMBER_SPACE,BasisSpace,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasis_TypeSet(BasisSpace,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL CMISSBasis_NumberOfXiSet(BasisSpace,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  CALL CMISSBasis_InterpolationXiSet(BasisSpace,(/BASIS_XI_INTERPOLATION_SPACE/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisSpace,(/BASIS_XI_GAUSS_SPACE/),Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(BasisSpace,Err)

  !Start the creation of VELOCITY basis
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisVelocity=BasisSpace
  ENDIF

  !Start the creation of AREA basis
  IF(BASIS_XI_INTERPOLATION_AREA==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisArea=BasisSpace
  ENDIF


  !================================================================================================================================
  ! Mesh
  !================================================================================================================================

  !Start the creation of mesh nodes
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSNodes_CreateStart(Region,TOTAL_NUMBER_OF_NODES,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMesh_NumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)

  !Specify spatial mesh component
  CALL CMISSMeshElements_Initialise(MeshElementsSpace,Err)
  CALL CMISSMeshElements_Initialise(MeshElementsVelocity,Err)
  CALL CMISSMeshElements_Initialise(MeshElementsArea,Err)

  MESH_COMPONENT_NUMBER_SPACE=1
  MESH_COMPONENT_NUMBER_VELOCITY=1
  MESH_COMPONENT_NUMBER_AREA=1

  ! NODES
  !             7-10-11
  !            /
  !           6
  !          /
  ! 1-2-3-4-5
  !          \   
  !           8 
  !            \
  !             9-12-13
  !
  ! ELEMENTS
  !             --5--
  !            /
  !           3
  !          /
  ! --1---2--
  !          \   
  !           4 
  !            \
  !             --6--


  CALL CMISSMeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_SPACE,BasisSpace,MeshElementsSpace,Err)

  CALL CMISSMeshElements_NodesSet(MeshElementsSpace,1,[1,2,3],Err)
  CALL CMISSMeshElements_NodesSet(MeshElementsSpace,2,[3,4,5],Err)
  CALL CMISSMeshElements_NodesSet(MeshElementsSpace,3,[5,6,7],Err)
  CALL CMISSMeshElements_NodesSet(MeshElementsSpace,4,[5,8,9],Err)
  CALL CMISSMeshElements_NodesSet(MeshElementsSpace,5,[7,10,11],Err)
  CALL CMISSMeshElements_NodesSet(MeshElementsSpace,6,[9,12,13],Err)
  !Bifurcation node
  !(meshElements,globalElementNumber,versionNumber,derivativeNumber,localElementNodeNumber,err)
  CALL CMISSMeshElements_LocalElementNodeVersionSet(MeshElementsSpace,2,1,CMISS_NO_GLOBAL_DERIV,3,Err) 
  CALL CMISSMeshElements_LocalElementNodeVersionSet(MeshElementsSpace,3,2,CMISS_NO_GLOBAL_DERIV,1,Err) 
  CALL CMISSMeshElements_LocalElementNodeVersionSet(MeshElementsSpace,4,3,CMISS_NO_GLOBAL_DERIV,1,Err) 

  IF(versionsFlag) THEN
    !set versions at 2 terminal vessels element 5 and 6, nodes 11 and 13- will be on same element
    CALL CMISSMeshElements_LocalElementNodeVersionSet(MeshElementsSpace,5,1,CMISS_NO_GLOBAL_DERIV,3,Err) 
    CALL CMISSMeshElements_LocalElementNodeVersionSet(MeshElementsSpace,5,2,CMISS_NO_GLOBAL_DERIV,3,Err) 
    CALL CMISSMeshElements_LocalElementNodeVersionSet(MeshElementsSpace,6,1,CMISS_NO_GLOBAL_DERIV,3,Err) 
    CALL CMISSMeshElements_LocalElementNodeVersionSet(MeshElementsSpace,6,2,CMISS_NO_GLOBAL_DERIV,3,Err) 
  ENDIF

  CALL CMISSMeshElements_CreateFinish(MeshElementsSpace,Err)

  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsVelocity=MeshElementsSpace
  ENDIF

  !Specify area mesh component
  IF(BASIS_XI_INTERPOLATION_AREA==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsArea=MeshElementsSpace
  ENDIF

  !Finish the creation of the mesh
  CALL CMISSMesh_CreateFinish(Mesh,Err)


  !================================================================================================================================
  !Decomposition
  !================================================================================================================================


  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,DomainUserNumber,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)


  !================================================================================================================================
  !Geometric Field
  !================================================================================================================================

  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field label
  CALL CMISSField_VariableLabelSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,'Coordinates',Err)
  !Set the field type
  CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL CMISSField_ScalingTypeSet(GeometricField,CMISS_FIELD_NO_SCALING,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,1
    CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Default values to 0
  X=0.0_CMISSDP
  Y=0.0_CMISSDP
  Z=0.0_CMISSDP
  DO nodeIdx=1,NUMBER_OF_NODES_SPACE
    !(field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value,err)
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & 1,CMISS_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & 1,CMISS_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
  ENDDO

  !parent vessel
  X=-0.5_CMISSDP
  Y=0.0_CMISSDP
  Z=0.0_CMISSDP
  DO nodeIdx=1,5
    X=X+0.5_CMISSDP
     !(field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value,err)
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & 1,CMISS_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & 1,CMISS_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
  ENDDO
  ! node 5 versions 2,3 (same x/y values)
  nodeIdx=5
  DO versionIdx=2,3
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
  ENDDO

  DO nodeIdx=6,7
    X=X+0.5_CMISSDP
    Y=Y+0.5_CMISSDP
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & 1,CMISS_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & 1,CMISS_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
  ENDDO

  DO nodeIdx=10,11
    X=X+0.5_CMISSDP
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & 1,CMISS_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & 1,CMISS_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
  ENDDO
  IF(versionsFlag) THEN
  ! node 11 versions 2 (same x/y values)
    nodeIdx=11
    versionIdx=2
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
  ENDIF

  X=2.0_CMISSDP
  Y=0.0_CMISSDP
  DO nodeIdx=8,9
    X=X+0.5_CMISSDP
    Y=Y-0.5_CMISSDP
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & 1,CMISS_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & 1,CMISS_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
  ENDDO

  DO nodeIdx=12,13
    X=X+0.5_CMISSDP
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & 1,CMISS_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & 1,CMISS_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
  ENDDO
  IF(versionsFlag) THEN
  ! node 13 versions 2 (same x/y values)
    nodeIdx=13
    versionIdx=2
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
  ENDIF

    
  !================================================================================================================================
  !Equations Sets
  !================================================================================================================================

  !Create the equations set for 1D Navier-Stokes
  CALL CMISSEquationsSet_Initialise(EquationsSetNavierStokes,Err)
  CALL CMISSField_Initialise(EquationsSetFieldNavierStokes,Err)
  !Set the equations set to be a dynamic Navier-Stokes problem
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField, &
    & CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMISS_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE, &
    & EquationsSetSubtype,EquationsSetFieldUserNumberNavierStokes, &
    & EquationsSetFieldNavierStokes,EquationsSetNavierStokes,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSetNavierStokes,Err)

  !Create the equations set for Characteristic
  CALL CMISSEquationsSet_Initialise(EquationsSetCharacteristic,Err)
  CALL CMISSField_Initialise(EquationsSetFieldCharacteristic,Err)
  !Set the equations set to be a static Nonlinear problem
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumberCharacteristic,Region,GeometricField, &
    & CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMISS_EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE, &
    & EquationsSetCharacteristicSubtype,EquationsSetFieldUserNumberCharacteristic, &
    & EquationsSetFieldCharacteristic,EquationsSetCharacteristic,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSetCharacteristic,Err)


  !================================================================================================================================
  !Dependent Field
  !================================================================================================================================

  !Create the equations set dependent field variables for static nonlinear Characteristic solver
  CALL CMISSField_Initialise(DependentFieldNavierStokes,Err)

  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetCharacteristic,DependentFieldUserNumber, & 
    & DependentFieldNavierStokes,Err)
  !Set the field label
  CALL CMISSField_VariableLabelSet(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,'General',Err)
  CALL CMISSField_VariableLabelSet(DependentFieldNavierStokes,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,'Derivatives',Err)
  CALL CMISSField_VariableLabelSet(DependentFieldNavierStokes,CMISS_FIELD_V_VARIABLE_TYPE,'Characteristics',Err)
  CALL CMISSField_VariableLabelSet(DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE,'calculated pressure',Err)
!  CALL CMISSField_VariableLabelSet(DependentFieldNavierStokes,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,'Coupling Derivatives',Err)

  !Set the mesh component to be used by the field components.
  COMPONENT_NUMBER=1 ! Velocity
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
    & MESH_COMPONENT_NUMBER_VELOCITY,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldNavierStokes,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
    & MESH_COMPONENT_NUMBER_VELOCITY,Err)

  COMPONENT_NUMBER=2 !Area
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
    & MESH_COMPONENT_NUMBER_AREA,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldNavierStokes,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, &
    & MESH_COMPONENT_NUMBER_AREA,Err)

  DO COMPONENT_NUMBER=1,2 ! W (Characteristics)
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldNavierStokes,CMISS_FIELD_V_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO

  IF (cellmlFlag) THEN
    COMPONENT_NUMBER=1 ! pCellML
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDIF

  !Finish the equations set dependent field variables for Characteristic equations set
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetCharacteristic,Err)

  !Create the equations set dependent field variables for dynamic Navier-Stokes
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetNavierStokes,DependentFieldUserNumber, & 
    & DependentFieldNavierStokes,Err)

  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetNavierStokes,Err)

  ! IF ((versionsFlag .EQV. .FALSE.) .AND. cellmlFlag) THEN
  !   CALL CMISSField_ParameterSetCreate(DependentFieldNavierStokes,CMISS_FIELD_V_VARIABLE_TYPE, &
  !    & CMISS_FIELD_RETURNING_CHARACTERISTIC_TYPE,err)
  ! ENDIF

  versionIdx=1
  componentIdx=1
  ! Initialize Q values
  DO nodeIdx=1,13
    IF(nodeIdx<6) THEN
      VALUE=Q1
    ELSE IF ((nodeIdx<13 .AND. nodeIdx>9) .OR. (nodeIdx<8 .AND. nodeIdx>5)) THEN
      VALUE=Q2
    ELSE
      VALUE=Q3
    ENDIF
    !(field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value,err)
    CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
  ENDDO
  versionIdx=1
  componentIdx=2
  ! Initialize A values
  DO nodeIdx=1,13
    IF(nodeIdx<6) THEN
      VALUE=A1
    ELSE IF ((nodeIdx<13 .AND. nodeIdx>9) .OR. (nodeIdx<8 .AND. nodeIdx>5)) THEN
      VALUE=A2
    ELSE
      VALUE=A3
    ENDIF
    !(field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value,err)
    CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
  ENDDO

  ! Branch versions at node 5
  nodeIdx=5
  versionIdx=2
  componentIdx=1
  VALUE=Q2  
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
  componentIdx=2
  VALUE=A2  
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
  versionIdx=3
  componentIdx=1
  VALUE=Q3  
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
  componentIdx=2
  VALUE=A3  
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
  !Initialize dependent characteristic field W (V_type) to 0
  nodeIdx=5
  VALUE=0.0_CMISSDP
  DO componentIdx=1,2
    DO versionIdx=1,3
      CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_V_VARIABLE_TYPE, &
        & CMISS_FIELD_VALUES_SET_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    ENDDO
  ENDDO

  IF(versionsFlag) THEN

    ! Branch versions at node 11
    nodeIdx=11
    versionIdx=2
    componentIdx=1
    VALUE=Q2  
    CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    componentIdx=2
    VALUE=A2  
    CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    !Initialize dependent characteristic field W (V_type) to 0
    VALUE=0.0_CMISSDP
    componentIdx=1
    DO versionIdx=1,2
      CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_V_VARIABLE_TYPE, &
        & CMISS_FIELD_VALUES_SET_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    ENDDO

    ! Branch versions at node 13
    nodeIdx=13
    versionIdx=2
    componentIdx=1
    VALUE=Q3  
    CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    componentIdx=2
    VALUE=A3  
    CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    !Initialize dependent characteristic field W (V_type) to 0
    VALUE=0.0_CMISSDP
    componentIdx=1
    DO versionIdx=1,2
      CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_V_VARIABLE_TYPE, &
        & CMISS_FIELD_VALUES_SET_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    ENDDO
    
  ENDIF


  !================================================================================================================================
  !Materials Field
  !================================================================================================================================

  !Create the equations set materials field variables for dynamic Navier-Stokes
  CALL CMISSField_Initialise(MaterialsFieldNavierStokes,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetNavierStokes,MaterialsFieldUserNumber, & 
    & MaterialsFieldNavierStokes,Err)
  !Set the field label
  CALL CMISSField_VariableLabelSet(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,'Materials',Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetNavierStokes,Err)

  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetCharacteristic,MaterialsFieldUserNumber, & 
    & MaterialsFieldNavierStokes,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetCharacteristic,Err)

  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberMu,MU_PARAM_NAVIER_STOKES,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberRho,RHO_PARAM_NAVIER_STOKES,Err)

  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberK,K,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberBs,Bs,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberAs,As,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberRe,Re,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberFr,Fr,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberSt,St,Err)
  DO i=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & MaterialsFieldUserNumberE,E_PARAM_NAVIER_STOKES,Err)
    CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & MaterialsFieldUserNumberH0,H0_PARAM_NAVIER_STOKES,Err)
    CALL CMISSField_ParameterSetUpdateElement(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & i,MaterialsFieldUserNumberA0,A0_PARAM(i),err)
    CALL CMISSField_ParameterSetUpdateElement(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & i,MaterialsFieldUserNumberBeta,Beta(i),Err)
  ENDDO


  IF(cellmlFlag) THEN

    ! CellML Materials Field parameters setup
    !--------------------------------------------------------------------------------------------------------------------
    !  pVesselWall: a constant expressing physical features of the vessel wall
    !  pExternal: pressure external to the vessel
    !
    ! WARNING: Do not change these component index values- they are what OpenCMISS uses to identify the variables
      pVesselWallComponent=1
      pExternalComponent=2
    !-------------------------------------------------------------------------------------------------------------------  

    !User defined variable values (feel free to change these- already initialised to 0 otherwise)
    pVesselWall=0.0_CMISSDP
    pExternal=0.0_CMISSDP

    ! Set values at coupled node 1
    !(field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value,err)
    CALL CMISSField_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,&
     & CMISS_NO_GLOBAL_DERIV,coupledNodeNumber1,pVesselWallComponent,pVesselWall,err)
    CALL CMISSField_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,&
     & CMISS_NO_GLOBAL_DERIV,coupledNodeNumber1,pExternalComponent,pExternal,err)

    ! Set values at coupled node 2
    !(field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value,err)
    CALL CMISSField_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,&
     & CMISS_NO_GLOBAL_DERIV,coupledNodeNumber2,pVesselWallComponent,pVesselWall,err)
    CALL CMISSField_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,&
     & CMISS_NO_GLOBAL_DERIV,coupledNodeNumber2,pExternalComponent,pExternal,err)

  ENDIF 


  !================================================================================================================================x
  ! Independent Field - Characteristic wave normal direction
  !================================================================================================================================x

  !Create the equations set independent field variables for Characteristics Solver
  CALL CMISSField_Initialise(IndependentFieldNavierStokes,Err)

  CALL CMISSEquationsSet_IndependentCreateStart(EquationsSetCharacteristic,IndependentFieldUserNumber, & 
    & IndependentFieldNavierStokes,Err)
  CALL CMISSField_VariableLabelSet(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,'Normal wave direction',Err)
  !Set the mesh component to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,1, & 
    & MESH_COMPONENT_NUMBER_SPACE,Err)
  !Finish the equations set independent field variables
  CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSetCharacteristic,Err)

  !(field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value,err)

  ! normalWave node 5
  nodeIdx=5
  componentIdx=1 ! Incoming normal
  ! 1 inlet/parent from element 2
  versionIdx=1
  VALUE=1.0_CMISSDP
  CALL CMISSField_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
     & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)

  ! 2 outlet/daughters from elements 3 and 4
  componentIdx=2 ! outgoing
  VALUE=-1.0_CMISSDP
  DO versionIdx=2,3
    CALL CMISSField_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
  ENDDO

  IF (versionsFlag) THEN

    ! normalDirection node 11
    nodeIdx=11
    componentIdx=1
    ! 1 inlet/parent from element 5
    VALUE=1.0_CMISSDP
    versionIdx=1
    CALL CMISSField_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE, &
       & CMISS_FIELD_VALUES_SET_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    ! 1 outlet/daughter from coupled domain
    VALUE=-1.0_CMISSDP
    versionIdx=2
    CALL CMISSField_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE, &
      & CMISS_FIELD_VALUES_SET_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)

    ! normalDirection node 13
    nodeIdx=13
    ! 1 inlet/parent from element 6
    VALUE=1.0_CMISSDP
    versionIdx=1
    CALL CMISSField_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE, &
       & CMISS_FIELD_VALUES_SET_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    ! 1 outlet/daughter from coupled domain
    VALUE=-1.0_CMISSDP
    versionIdx=2
    CALL CMISSField_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE, &
      & CMISS_FIELD_VALUES_SET_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)

  ELSE IF ((versionsFlag .EQV. .FALSE.) .AND. cellmlFlag) THEN

    ! Incoming normals for 11,13
    nodeIdx=11
    componentIdx=1
    VALUE=1.0_CMISSDP
    versionIdx=1
    CALL CMISSField_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE, &
       & CMISS_FIELD_VALUES_SET_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    nodeIdx=13
    CALL CMISSField_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE, &
       & CMISS_FIELD_VALUES_SET_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)

    ! Outgoing normals for 11,13
    nodeIdx=11
    componentIdx=2
    VALUE=-1.0_CMISSDP
    versionIdx=1
    CALL CMISSField_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE, &
       & CMISS_FIELD_VALUES_SET_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    nodeIdx=13
    CALL CMISSField_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE, &
       & CMISS_FIELD_VALUES_SET_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)

  ENDIF

  !Create the equations set independent field variables for 1D Transient Solver
  CALL CMISSEquationsSet_IndependentCreateStart(EquationsSetNavierStokes,IndependentFieldUserNumber, & 
    & IndependentFieldNavierStokes,Err)
  !Finish the equations set independent field variables
  CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSetNavierStokes,Err)


  !================================================================================================================================
  !  C e l l M L    M o d e l    M a p s
  !=====================================ppp===========================================================================================

  IF (cellmlFlag) THEN

    !------------------------------------------------------------------------------------------------------------------------------
    ! Description
    !------------------------------------------------------------------------------------------------------------------------------
    ! A CellML OD model is used to provide the impedance from the downstream vascular bed beyond the termination
    ! point of the 1D model. This is iteratively coupled with the the 1D solver. In the case of a simple resistance
    ! model, P=RQ, which is analogous to Ohm's law: V=IR. A variable map copies the guess for the FlowRate, Q at 
    ! the boundary from the OpenCMISS Dependent Field to the CellML equation, which then returns presssure, P.
    ! The initial guess value for Q is taken from the previous time step or is 0 for t=0.  In OpenCMISS this P value is 
    ! then used to compute a new Area value based on the P-A relationship and the Riemann variable W_2, which gives a
    ! new value for Q until the values for Q and P converge within tolerance of the previous value.
    !------------------------------------------------------------------------------------------------------------------------------

     pCellMLComponent=1

     ! --- W i n d k e s s e l   M o d e l --- !
     IF (windkesselFlag) THEN

      !Create the CellML environment
      CALL CMISSCellML_Initialise(CellML,Err)
      CALL CMISSCellML_CreateStart(CellMLUserNumber,Region,CellML,Err)

      !Import an RCR windkessel model
      CALL CMISSCellML_ModelImport(CellML,"windkessel.xml",WindkesselModelIndex,Err)    

      ! - known (to OpenCMISS) variables 
      CALL CMISSCellML_VariableSetAsKnown(CellML,WindkesselModelIndex,"equations/Q",Err)
      CALL CMISSCellML_VariableSetAsKnown(CellML,WindkesselModelIndex,"equations/R_p",Err)
      CALL CMISSCellML_VariableSetAsKnown(CellML,WindkesselModelIndex,"equations/R_d",Err)
      CALL CMISSCellML_VariableSetAsKnown(CellML,WindkesselModelIndex,"equations/C",Err)
      CALL CMISSCellML_VariableSetAsKnown(CellML,WindkesselModelIndex,"equations/t",Err)
      ! - to get from the CellML side 
      CALL CMISSCellML_VariableSetAsWanted(CellML,WindkesselModelIndex,"equations/P",Err)

      CALL CMISSCellML_CreateFinish(CellML,Err)

      !Start the creation of CellML <--> OpenCMISS field maps
      CALL CMISSCellML_FieldMapsCreateStart(CellML,Err)
      !Now we can set up the field variable component <--> CellML model variable mappings.

      !Map the OpenCMISS boundary flow rate values --> CellML
      ! Q is component 1 of the DependentField
      CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_FIELD_VALUES_SET_TYPE,WindkesselModelIndex,"equations/Q",CMISS_FIELD_VALUES_SET_TYPE,Err)
      CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_FIELD_VALUES_SET_TYPE,WindkesselModelIndex,"equations/t",CMISS_FIELD_VALUES_SET_TYPE,Err)
      !Map the returned pressure values from CellML --> CMISS
      ! pCellML is component 1 of the Dependent field U1 variable
      CALL CMISSCellML_CreateCellMLToFieldMap(CellML,WindkesselModelIndex,"equations/Pressure",CMISS_FIELD_VALUES_SET_TYPE, &
        & DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE,pCellMLComponent,CMISS_FIELD_VALUES_SET_TYPE,Err)

      !Finish the creation of CellML <--> OpenCMISS field maps
      CALL CMISSCellML_FieldMapsCreateFinish(CellML,Err)

      !Create the CellML models field --- only 1 model here
      CALL CMISSField_Initialise(CellMLModelsField,Err)
      CALL CMISSCellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
      CALL CMISSCellML_ModelsFieldCreateFinish(CellML,Err)

      !Create the CellML parameters field --- will be the Resistance and Flow rate
      CALL CMISSField_Initialise(CellMLParametersField,Err)
      CALL CMISSCellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
      CALL CMISSCellML_ParametersFieldCreateFinish(CellML,Err)

      !Create the CellML intermediate field --- will be the pressure value returned from CellML to be used for 
      ! recalculation of the incoming Riemann variable W(2)
      CALL CMISSField_Initialise(CellMLIntermediateField,Err)
      CALL CMISSCellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
      CALL CMISSCellML_IntermediateFieldCreateFinish(CellML,Err)

      ! Initialise pCellML (and previous pCellML coupling iteration values) values to 0 at the outlet nodes
      pCellML=0.0_CMISSDP
      pPrevious=0.0_CMISSDP
      !(field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value,err)
      CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE, &
       & CMISS_FIELD_VALUES_SET_TYPE,1,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber1,pCellMLComponent,pCellML,err)
      CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE, &
       & CMISS_FIELD_VALUES_SET_TYPE,1,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber2,pCellMLComponent,pCellML,err)
      CALL CMISSField_ParameterSetCreate(DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE, &
       & CMISS_FIELD_PREVIOUS_VALUES_SET_TYPE,err)
      CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE, &
       & CMISS_FIELD_PREVIOUS_VALUES_SET_TYPE,1,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber1,pCellMLComponent,pPrevious,err)
      CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE, &
       & CMISS_FIELD_PREVIOUS_VALUES_SET_TYPE,1,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber2,pCellMLComponent,pPrevious,err)


     ! --- R e s i s t a n c e   M o d e l --- !
     ELSE 

      !Create the CellML environment
      CALL CMISSCellML_Initialise(CellML,Err)
      CALL CMISSCellML_CreateStart(CellMLUserNumber,Region,CellML,Err)

      !Import a simple resistance model P = RQ, analogous to V = IR
      CALL CMISSCellML_ModelImport(CellML,"resistance.xml",ResistanceModelIndex,Err)

      ! - known (to OpenCMISS) variables 
      CALL CMISSCellML_VariableSetAsKnown(CellML,ResistanceModelIndex,"equations/FlowRate",Err)
      CALL CMISSCellML_VariableSetAsKnown(CellML,ResistanceModelIndex,"equations/Resistance",Err)
      ! - to get from the CellML side 
      CALL CMISSCellML_VariableSetAsWanted(CellML,ResistanceModelIndex,"equations/Pressure",Err)

      CALL CMISSCellML_CreateFinish(CellML,Err)

      !Start the creation of CellML <--> OpenCMISS field maps
      CALL CMISSCellML_FieldMapsCreateStart(CellML,Err)
      !Now we can set up the field variable component <--> CellML model variable mappings.

      !Map the OpenCMISS boundary flow rate values --> CellML
      ! Q is component 1 of the DependentField
      CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,1, &
        & CMISS_FIELD_VALUES_SET_TYPE,ResistanceModelIndex,"equations/FlowRate",CMISS_FIELD_VALUES_SET_TYPE,Err)
      !Map the returned pressure values from CellML --> CMISS
      ! pCellML is component 2 of the Dependent field V variable
      CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ResistanceModelIndex,"equations/Pressure",CMISS_FIELD_VALUES_SET_TYPE, &
        & DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE,pCellMLComponent,CMISS_FIELD_VALUES_SET_TYPE,Err)

      !Finish the creation of CellML <--> OpenCMISS field maps
      CALL CMISSCellML_FieldMapsCreateFinish(CellML,Err)

      !Create the CellML models field --- only 1 model here
      CALL CMISSField_Initialise(CellMLModelsField,Err)
      CALL CMISSCellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
      CALL CMISSCellML_ModelsFieldCreateFinish(CellML,Err)

      !Create the CellML parameters field --- will be the Resistance and Flow rate
      CALL CMISSField_Initialise(CellMLParametersField,Err)
      CALL CMISSCellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
      CALL CMISSCellML_ParametersFieldCreateFinish(CellML,Err)

      !Create the CellML intermediate field --- will be the pressure value returned from CellML to be used for 
      ! recalculation of the incoming Riemann variable W(2)
      CALL CMISSField_Initialise(CellMLIntermediateField,Err)
      CALL CMISSCellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
      CALL CMISSCellML_IntermediateFieldCreateFinish(CellML,Err)

      ! Initialise pCellML (and previous pCellML coupling iteration values) values to 0 at the outlet nodes
      pCellML=0.0_CMISSDP
      pPrevious=0.0_CMISSDP
      !(field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value,err)
      CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE, &
       & CMISS_FIELD_VALUES_SET_TYPE,1,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber1,pCellMLComponent,pCellML,err)
      CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE, &
       & CMISS_FIELD_VALUES_SET_TYPE,1,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber2,pCellMLComponent,pCellML,err)
      CALL CMISSField_ParameterSetCreate(DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE, &
       & CMISS_FIELD_PREVIOUS_VALUES_SET_TYPE,err)
      CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE, &
       & CMISS_FIELD_PREVIOUS_VALUES_SET_TYPE,1,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber1,pCellMLComponent,pPrevious,err)
      CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U1_VARIABLE_TYPE, &
       & CMISS_FIELD_PREVIOUS_VALUES_SET_TYPE,1,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber2,pCellMLComponent,pPrevious,err)

    ENDIF ! windkessel/resistance
 
  ENDIF ! cellml flag

  !================================================================================================================================
  !Equations
  !================================================================================================================================

  !Create the equations set equations
  CALL CMISSEquations_Initialise(EquationsNavierStokes,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetNavierStokes,EquationsNavierStokes,Err)
  !Set the equations matrices sparsity type
  !CALL CMISSEquations_SparsityTypeSet(EquationsNavierStokes,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_SparsityTypeSet(EquationsNavierStokes,CMISS_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations lumping type
  CALL CMISSEquations_LumpingTypeSet(EquationsNavierStokes,CMISS_EQUATIONS_UNLUMPED_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(EquationsNavierStokes,EQUATIONS_NAVIER_STOKES_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetNavierStokes,Err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(EquationsCharacteristic,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetCharacteristic,EquationsCharacteristic,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(EquationsCharacteristic,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(EquationsCharacteristic,EQUATIONS_NAVIER_STOKES_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetCharacteristic,Err)


  !================================================================================================================================
  !Problems
  !================================================================================================================================

  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a dynamic Navier-Stokes problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_FLUID_MECHANICS_CLASS,CMISS_PROBLEM_NAVIER_STOKES_EQUATION_TYPE, &
    & ProblemSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoop_TimesSet(ControlLoop,DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME,DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME, & 
    & DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT,Err)
  !Set the output timing
  CALL CMISSControlLoop_TimeOutputSet(ControlLoop,DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)


  !================================================================================================================================
  !Solvers
  !================================================================================================================================

  !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(DynamicSolverNavierStokes,Err)
  CALL CMISSSolver_Initialise(NonlinearSolverNavierStokes,Err)
  CALL CMISSSolver_Initialise(LinearSolverNavierStokes,Err)

  CALL CMISSSolver_Initialise(NonlinearSolverCharacteristic,Err)
  CALL CMISSSolver_Initialise(LinearSolverCharacteristic,Err)

  CALL CMISSProblem_SolversCreateStart(Problem,Err)

  !Get the dynamic dynamic solver
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(DynamicSolverNavierStokes,DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set theta
  CALL CMISSSolver_DynamicThetaSet(DynamicSolverNavierStokes,DYNAMIC_SOLVER_NAVIER_STOKES_THETA,Err)
  !Get the dynamic nonlinear solver
  CALL CMISSSolver_DynamicNonlinearSolverGet(DynamicSolverNavierStokes,NonlinearSolverNavierStokes,Err)
  !Set the nonlinear Jacobian type
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes, &
    & CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
!  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(NonlinearSolverNavierStokes,NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(NonlinearSolverNavierStokes,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(NonlinearSolverNavierStokes,RELATIVE_TOLERANCE,Err)
  !Get the dynamic nonlinear linear solver
  CALL CMISSSolver_NewtonLinearSolverGet(NonlinearSolverNavierStokes,LinearSolverNavierStokes,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(LinearSolverNavierStokes,LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG) THEN
    CALL CMISSSolver_LinearTypeSet(LinearSolverNavierStokes,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL CMISSSolver_LibraryTypeSet(LinearSolverNavierStokes,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL CMISSSolver_LinearTypeSet(LinearSolverNavierStokes,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolverNavierStokes,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolver_LinearIterativeDivergenceToleranceSet(LinearSolverNavierStokes,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(LinearSolverNavierStokes,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(LinearSolverNavierStokes,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeGMRESRestartSet(LinearSolverNavierStokes,RESTART_VALUE,Err)
  ENDIF

  !Get the static Nonlinear solver
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverCharacteristicUserNumber,NonlinearSolverCharacteristic,Err)
  !Set the nonlinear Jacobian type
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolverCharacteristic, &
    & CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
!  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolverCharacteristic,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(NonlinearSolverCharacteristic,NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(NonlinearSolverCharacteristic,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(NonlinearSolverCharacteristic,RELATIVE_TOLERANCE,Err)
  !Get the nonlinear linear solver
  CALL CMISSSolver_NewtonLinearSolverGet(NonlinearSolverCharacteristic,LinearSolverCharacteristic,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(LinearSolverCharacteristic,LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG) THEN
    CALL CMISSSolver_LinearTypeSet(LinearSolverCharacteristic,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL CMISSSolver_LibraryTypeSet(LinearSolverCharacteristic,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL CMISSSolver_LinearTypeSet(LinearSolverCharacteristic,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolverCharacteristic,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolver_LinearIterativeDivergenceToleranceSet(LinearSolverCharacteristic,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(LinearSolverCharacteristic,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(LinearSolverCharacteristic,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeGMRESRestartSet(LinearSolverCharacteristic,RESTART_VALUE,Err)
  ENDIF

  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)


  ! C e l l M L   S o l v e r
  !-----------------------------
  IF (cellmlFlag) THEN

    !Create the problem solver CellML equations
    CALL CMISSSolver_Initialise(CellMLSolver,Err)
    CALL CMISSCellMLEquations_Initialise(CellMLEquations,Err)
    CALL CMISSProblem_CellMLEquationsCreateStart(Problem,Err)

    CALL CMISSSolver_NewtonCellMLSolverGet(DynamicSolverNavierStokes,CellMLSolver,Err) 
    CALL CMISSSolver_CellMLEquationsGet(CellMLSolver,CellMLEquations,Err)
    CALL CMISSCellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
    CALL CMISSProblem_CellMLEquationsCreateFinish(Problem,Err)

  ENDIF


  !================================================================================================================================
  !Solver Equations
  !================================================================================================================================

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(DynamicSolverNavierStokes,Err)
  CALL CMISSSolver_Initialise(NonlinearSolverCharacteristic,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsNavierStokes,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsCharacteristic,Err)

  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)

  !Get the dynamic Navier-Stokes solver equations
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  CALL CMISSSolver_SolverEquationsGet(DynamicSolverNavierStokes,SolverEquationsNavierStokes,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsNavierStokes,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsNavierStokes,EquationsSetNavierStokes,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations

  !Get the static Nonlinear solver equations
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverCharacteristicUserNumber,NonlinearSolverCharacteristic,Err)
  CALL CMISSSolver_SolverEquationsGet(NonlinearSolverCharacteristic,SolverEquationsCharacteristic,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsCharacteristic,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsCharacteristic,EquationsSetCharacteristic,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations

  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)


  !================================================================================================================================
  !Boundary Conditions
  !================================================================================================================================

  !Start the creation of the equations set boundary conditions for Navier-Stokes
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsNavierStokes,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsNavierStokes,BoundaryConditionsNavierStokes,Err)
  !Set velocity boundary conditions
  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    CONDITION=CMISS_BOUNDARY_CONDITION_FIXED_INLET
    VALUE1=BOUNDARY_CONDITIONS_NAVIER_STOKES(1)
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
      & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV,1,1,CONDITION,VALUE1,Err)
  ENDIF
  !Set area boundary conditions if not a CellML problem
  IF(cellmlFlag) THEN
    CONDITION=CMISS_BOUNDARY_CONDITION_FIXED_OUTLET
    VALUE1=A2
    versionIdx=1
    !(boundaryConditions,field,variableType,versionNumber,derivativeNumber,nodeUserNumber,componentNumber,condition,value,err)
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
      & CMISS_FIELD_U_VARIABLE_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber1,2,CONDITION,VALUE1,Err)
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
      & CMISS_FIELD_U_VARIABLE_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber2,2,CONDITION,VALUE1,Err)
  ELSE
    CONDITION=CMISS_BOUNDARY_CONDITION_FIXED
    VALUE1=A2    
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
      & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber1,2,CONDITION,VALUE1,Err)
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
      & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber2,2,CONDITION,VALUE1,Err)
    IF (versionsFlag) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
        & CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber1,2,CONDITION,VALUE1,Err)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
        & CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber2,2,CONDITION,VALUE1,Err)
    ENDIF
  ENDIF

  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsCharacteristic,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsCharacteristic,BoundaryConditionsCharacteristic,Err)
  IF(cellmlFlag) THEN
    CONDITION=CMISS_BOUNDARY_CONDITION_FIXED_OUTLET
    VALUE1=A2
    versionIdx=1
    !(boundaryConditions,field,variableType,versionNumber,derivativeNumber,nodeUserNumber,componentNumber,condition,value,err)
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsCharacteristic,DependentFieldNavierStokes, &
      & CMISS_FIELD_U_VARIABLE_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber1,2,CONDITION,VALUE1,Err)
    VALUE1=A3
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsCharacteristic,DependentFieldNavierStokes, &
      & CMISS_FIELD_U_VARIABLE_TYPE,versionIdx,CMISS_NO_GLOBAL_DERIV,coupledNodeNumber2,2,CONDITION,VALUE1,Err)
  ENDIF

  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsNavierStokes,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsCharacteristic,Err)

  !================================================================================================================================
  ! C e l l M L   P a r a m e t e r s
  !================================================================================================================================

  IF (cellmlFlag) THEN

    ! Set CellML model parameters (Resistance/Capacitance) at boundary nodes
    CALL CMISSCellML_FieldComponentGet(CellML,ResistanceModelIndex,CMISS_CELLML_PARAMETERS_FIELD,"equations/Resistance", &
      & resistanceComponent,Err)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,coupledNodeNumber1,1,BoundaryNodeDomain,Err)
    IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
      !Branch 1 R=R
      CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,&
        & coupledNodeNumber1,resistanceComponent,resistanceProximal,Err)
    ENDIF
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,coupledNodeNumber2,1,BoundaryNodeDomain,Err)
    IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
      !Branch 2 
      resistanceProximal= resistanceProximal*0.5_CMISSDP
      CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,&
        & coupledNodeNumber2,resistanceComponent,resistanceProximal,Err)
    ENDIF

  ENDIF

  !================================================================================================================================
  ! RUN SOLVERS
  !================================================================================================================================

  !Turn of PETSc error handling
!  CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL CMISSProblem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"

  !
  !================================================================================================================================
  !

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM FortranExample
