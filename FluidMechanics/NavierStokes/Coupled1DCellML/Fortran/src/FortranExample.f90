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
  USE OpenCMISS
  USE OpenCMISS_Iron
  USE FLUID_MECHANICS_IO_ROUTINES
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWINCMISS
#endif


  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !================================================================================================================================
  !Test program parameters
  !================================================================================================================================

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberNavierStokes=6
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberCharacteristic=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberNavierStokes=11
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberCharacteristic=12
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: Basis_Lagrange_Hermite_TP_Type=1
  INTEGER(CMISSIntg), PARAMETER :: Basis_Quadratic_Lagrange_Interpolation=2

  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberCellML=19

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMu=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberRho=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberK=3
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberAs=4
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberRe=5
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberFr=6
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberSt=7
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberA0=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberE=9
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberH0=10

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

  REAL(CMISSRP) :: BOUNDARY_CONDITIONS_NAVIER_STOKES(2)
  REAL(CMISSRP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSRP) :: RELATIVE_TOLERANCE
  REAL(CMISSRP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSRP) :: LINESEARCH_ALPHA
  REAL(CMISSRP) :: VALUE1
  REAL(CMISSRP) :: MU_PARAM_NAVIER_STOKES
  REAL(CMISSRP) :: RHO_PARAM_NAVIER_STOKES
  REAL(CMISSRP) :: E_PARAM_NAVIER_STOKES(13)
  REAL(CMISSRP) :: H0_PARAM_NAVIER_STOKES(13)
  REAL(CMISSRP) :: A0_PARAM(13),Q1,A1,Ts,As,Qs,Xs,St,Fr,Re,Bs,K
  REAL(CMISSRP) :: Q2,Q3,A2,A3
  REAL(CMISSRP) :: resistanceProximal,resistanceDistal,capacitance
  REAL(CMISSRP) :: pCellML,pPrevious,pVesselWall,pExternal,qPrevious
  REAL(CMISSRP) :: length, position
  REAL(CMISSRP) :: VALUE,X,Y,Z

  REAL(CMISSRP) :: DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME
  REAL(CMISSRP) :: DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME
  REAL(CMISSRP) :: DYNAMIC_SOLVER_NAVIER_STOKES_THETA
  REAL(CMISSRP) :: DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT

  INTEGER(CMISSIntg) :: pCellMLComponent,pPreviousComponent,pVesselWallComponent,pExternalComponent,qPreviousComponent
  INTEGER(CMISSIntg) :: modelIndex,numberOfCellmlModels
  INTEGER(CMISSIntg) :: EquationsSetSubtype
  INTEGER(CMISSIntg) :: EquationsSetCharacteristicSubtype
  INTEGER(CMISSIntg) :: ProblemSubtype
  INTEGER(CMISSIntg) :: coupledNodeNumber,coupledNodeNumber1,coupledNodeNumber2
  INTEGER(CMISSIntg) :: resistanceComponent,numberOfCoordinateDimensions
  INTEGER(CMISSIntg) :: nodeIdx,versionIdx,componentIdx,normalWave
  INTEGER(CMISSIntg) :: numberOfNodes

  LOGICAL :: LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG
  LOGICAL :: OUTLET_WALL_NODES_NAVIER_STOKES_FLAG
  LOGICAL :: INLET_WALL_NODES_NAVIER_STOKES_FLAG
  LOGICAL :: cellmlFlag,windkesselFlag,mixedFlag

  INTEGER(CMISSIntg) :: SolverDaeUserNumber
  INTEGER(CMISSIntg) :: SolverCharacteristicUserNumber
  INTEGER(CMISSIntg) :: SolverNavierStokesUserNumber

  !CMISS variables

  !Regions
  TYPE(cmfe_RegionType) :: Region
  TYPE(cmfe_RegionType) :: WorldRegion
  !Coordinate systems
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
  TYPE(cmfe_CoordinateSystemType) :: WorldCoordinateSystem
  !Basis
  TYPE(cmfe_BasisType) :: BasisSpace
  TYPE(cmfe_BasisType) :: BasisVelocity
  TYPE(cmfe_BasisType) :: BasisArea
  !Nodes
  TYPE(cmfe_NodesType) :: Nodes
  !Elements
  TYPE(cmfe_MeshElementsType) :: MeshElementsSpace
  TYPE(cmfe_MeshElementsType) :: MeshElementsVelocity
  TYPE(cmfe_MeshElementsType) :: MeshElementsArea
  !Meshes
  TYPE(cmfe_MeshType) :: Mesh
  !Decompositions
  TYPE(cmfe_DecompositionType) :: Decomposition
  !Field types
  TYPE(cmfe_FieldType) :: GeometricField
  TYPE(cmfe_FieldType) :: EquationsSetFieldNavierStokes
  TYPE(cmfe_FieldType) :: EquationsSetFieldCharacteristic
  TYPE(cmfe_FieldType) :: DependentFieldNavierStokes
  TYPE(cmfe_FieldType) :: MaterialsFieldNavierStokes
  TYPE(cmfe_FieldType) :: IndependentFieldNavierStokes
  ! CellML types
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_FieldType) :: MaterialsFieldCellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_FieldType) :: CellMLModelsField
  TYPE(cmfe_FieldType) :: CellMLStateField
  TYPE(cmfe_FieldType) :: CellMLIntermediateField
  TYPE(cmfe_FieldType) :: CellMLParametersField
  TYPE(cmfe_SolverType) :: CellMLSolver
  !Boundary conditions
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsNavierStokes
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsCharacteristic
  !Equations sets
  TYPE(cmfe_EquationsSetType) :: EquationsSetNavierStokes
  TYPE(cmfe_EquationsSetType) :: EquationsSetCharacteristic
  !Equations
  TYPE(cmfe_EquationsType) :: EquationsNavierStokes
  TYPE(cmfe_EquationsType) :: EquationsCharacteristic
  !Problems
  TYPE(cmfe_ProblemType) :: Problem
  !Control loops
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  !Solvers
  TYPE(cmfe_SolverType) :: DynamicSolverNavierStokes
  TYPE(cmfe_SolverType) :: NonlinearSolverNavierStokes
  TYPE(cmfe_SolverType) :: LinearSolverNavierStokes
  TYPE(cmfe_SolverType) :: NonlinearSolverCharacteristic
  TYPE(cmfe_SolverType) :: LinearSolverCharacteristic
  !Solver equations
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsNavierStokes
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsCharacteristic

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeDomain,BoundaryNodeDomain
  INTEGER(CMISSIntg) :: CellMLIndex
  INTEGER(CMISSIntg) :: CellMLModelIndex,CellMLModelIndex1,CellMLModelIndex2
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

  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
!  CALL cmfe_DiagnosticsSetOn(CMFE_ALL_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",[""],Err)
!  CALL cmfe_OutputSetOn("Testing",Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !================================================================================================================================
  !PROBLEM CONTROL PANEL
  !================================================================================================================================

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

  ! Set true to use CellML models
  !cellmlFlag = .TRUE.
  cellmlFlag = .FALSE.
  ! Set true to use a windkessel DAE solver
  windkesselFlag = .TRUE.
  ! Set true to set one boundary to resistance, one to windkessel
  mixedFlag = .FALSE.

  numberOfCoordinateDimensions=2
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
  MU_PARAM_NAVIER_STOKES=0.0035_CMISSRP !Mu(Pa.s)
  RHO_PARAM_NAVIER_STOKES=1050.0_CMISSRP !Rho(kg/m3)
  E_PARAM_NAVIER_STOKES(:)=0.8E+6_CMISSRP !Elasticity(Pa)
  H0_PARAM_NAVIER_STOKES(:)=0.5E-3_CMISSRP !Wall Thickness(m)
  DO i=1,5
    A0_PARAM(i)=19.6e-6_CMISSRP  !Wall Area(m2)
  ENDDO
  DO i=6,13
    A0_PARAM(i)=12.8e-6_CMISSRP
  ENDDO

  !Reference Values
  Qs=10.0e-6_CMISSRP !(m3/s)
  As=19.6e-6_CMISSRP !(m2)
!  As=100e-6_CMISSRP !(m2)
  Xs=0.1_CMISSRP     !(m)
  Ts=0.1_CMISSRP     !(s)
  K=4.0_CMISSRP/3.0_CMISSRP  !Parabolic Flow Section
!  Bs=(4.0*1.7725*E_PARAM_NAVIER_STOKES*H0_PARAM_NAVIER_STOKES)/(3.0*As)2
!  Bs=((4.0_CMISSRP)*(3.1416_CMISSRP**(0.5_CMISSRP))*E_PARAM_NAVIER_STOKES(1)*H0_PARAM_NAVIER_STOKES(1))/((3.0_CMISSRP)*As)
  St=(As*Xs)/(Ts*Qs)
!  Fr=((As**(2.5_CMISSRP))/(Qs**2))*(Bs/((2.0_CMISSRP)*RHO_PARAM_NAVIER_STOKES))
  Fr=((As**2.5_CMISSRP)/(Qs**2.0_CMISSRP))/(2.0_CMISSRP*RHO_PARAM_NAVIER_STOKES)
  Re=8.0_CMISSRP*3.1416_CMISSRP*MU_PARAM_NAVIER_STOKES*Xs/(Qs*RHO_PARAM_NAVIER_STOKES)

  ! Qs=1.0_CMISSRP
  ! As=1.0_CMISSRP
  ! Xs=1.0_CMISSRP
  ! Ts=1.0_CMISSRP
  ! K=4.0_CMISSRP/3.0_CMISSRP  !Parabolic Flow Section
  ! St=(As*Xs)/(Ts*Qs)
  ! Fr=1.0_CMISSRP
  ! Re=1.0_CMISSRP

  !Set initial condition
  Q1=7.0_CMISSRP
  Q2=3.5_CMISSRP
  Q3=3.5_CMISSRP

  A1=1.0_CMISSRP
  A2=0.653_CMISSRP
  A3=0.653_CMISSRP

  OUTLET_WALL_NODES_NAVIER_STOKES_FLAG=.TRUE.
  INLET_WALL_NODES_NAVIER_STOKES_FLAG=.TRUE.

  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_NAVIER_STOKES(1)=1.0_CMISSRP
  ENDIF

  IF(OUTLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
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
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMFE_SOLVER_PROGRESS_OUTPUT
  NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMFE_SOLVER_SOLVER_OUTPUT
  LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMFE_SOLVER_PROGRESS_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_NAVIER_STOKES_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT
  !Set time parameter
  DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME=0.0_CMISSRP
!  DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME=10.00001_CMISSRP
  DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME=1.00000001_CMISSRP
  DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT=0.001_CMISSRP
  DYNAMIC_SOLVER_NAVIER_STOKES_THETA=1.0_CMISSRP/2.0_CMISSRP
  !Set result output parameter
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-5_CMISSRP !default: 1.0E-05_CMISSRP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-10_CMISSRP
  DIVERGENCE_TOLERANCE=1.0E20 !default: 1.0E5
  MAXIMUM_ITERATIONS=10000 !default: 100000
  RESTART_VALUE=300 !default: 30
  LINESEARCH_ALPHA=1.0

  IF(cellmlFlag) THEN
    ! New equations set type to store p values in the Equations Set Field
    EquationsSetSubtype=CMFE_EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE
    ! Characteristic (nodal/characteristic) solver remains the same
    EquationsSetCharacteristicSubtype=CMFE_EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE
    IF(windkesselFlag) THEN
      SolverDAEUserNumber=1
      SolverCharacteristicUserNumber=2
      SolverNavierStokesUserNumber=3
      ProblemSubtype=CMFE_PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE
    ELSE
      SolverDAEUserNumber=0
      SolverCharacteristicUserNumber=1
      SolverNavierStokesUserNumber=2
      ! New problem type to execute the 1D-0D coupling subloop at each timestep
      ProblemSubtype=CMFE_PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE
    ENDIF
  ELSE
    SolverDAEUserNumber=0
    SolverCharacteristicUserNumber=1
    SolverNavierStokesUserNumber=2
    EquationsSetSubtype=CMFE_EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE
    EquationsSetCharacteristicSubtype=CMFE_EQUATIONS_SET_CHARACTERISTIC_SUBTYPE
    ProblemSubtype=CMFE_PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE
  ENDIF

  ! Default cellml model-specific flags to false if CellML is turned off
  IF (.NOT. cellmlFlag) THEN
    windkesselFlag = .FALSE.
    mixedFlag = .FALSE.
  ENDIF

  !================================================================================================================================
  !Coordinate System
  !================================================================================================================================

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,numberOfCoordinateDimensions,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !================================================================================================================================
  !Region
  !================================================================================================================================

  !Start the creation of a new region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the region label
  CALL cmfe_Region_LabelSet(Region,'OpenCMISS',Err)
  !Set the regions coordinate system as defined above
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)


  !================================================================================================================================
  !Bases
  !================================================================================================================================

  !Start the creation of SPACE bases
  MESH_NUMBER_OF_COMPONENTS=1
  CALL cmfe_Basis_Initialise(BasisSpace,Err)
  CALL cmfe_Basis_CreateStart(BASIS_NUMBER_SPACE,BasisSpace,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL cmfe_Basis_TypeSet(BasisSpace,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL cmfe_Basis_NumberOfXiSet(BasisSpace,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  CALL cmfe_Basis_InterpolationXiSet(BasisSpace,[BASIS_XI_INTERPOLATION_SPACE],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisSpace,[BASIS_XI_GAUSS_SPACE],Err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(BasisSpace,Err)

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
  CALL cmfe_MeshElements_Initialise(MeshElementsArea,Err)

  MESH_COMPONENT_NUMBER_SPACE=1
  MESH_COMPONENT_NUMBER_VELOCITY=1
  MESH_COMPONENT_NUMBER_AREA=1

  CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_SPACE,BasisSpace,MeshElementsSpace,Err)

  CALL cmfe_MeshElements_NodesSet(MeshElementsSpace,1,[1,2,3],Err)
  CALL cmfe_MeshElements_NodesSet(MeshElementsSpace,2,[3,4,5],Err)
  CALL cmfe_MeshElements_NodesSet(MeshElementsSpace,3,[5,6,7],Err)
  CALL cmfe_MeshElements_NodesSet(MeshElementsSpace,4,[5,8,9],Err)
  CALL cmfe_MeshElements_NodesSet(MeshElementsSpace,5,[7,10,11],Err)
  CALL cmfe_MeshElements_NodesSet(MeshElementsSpace,6,[9,12,13],Err)
  !Bifurcation node
  !(meshElements,globalElementNumber,versionNumber,derivativeNumber,localElementNodeNumber,err)
  CALL cmfe_MeshElements_LocalElementNodeVersionSet(MeshElementsSpace,2,1,CMFE_NO_GLOBAL_DERIV,3,Err) 
  CALL cmfe_MeshElements_LocalElementNodeVersionSet(MeshElementsSpace,3,2,CMFE_NO_GLOBAL_DERIV,1,Err) 
  CALL cmfe_MeshElements_LocalElementNodeVersionSet(MeshElementsSpace,4,3,CMFE_NO_GLOBAL_DERIV,1,Err) 

  CALL cmfe_MeshElements_CreateFinish(MeshElementsSpace,Err)

  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsVelocity=MeshElementsSpace
  ENDIF

  !Specify area mesh component
  IF(BASIS_XI_INTERPOLATION_AREA==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsArea=MeshElementsSpace
  ENDIF

  !Finish the creation of the mesh
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)


  !================================================================================================================================
  !Decomposition
  !================================================================================================================================

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)


  !================================================================================================================================
  !Geometric Field
  !================================================================================================================================

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field label
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,'Coordinates',Err)
  !Set the field type
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL cmfe_Field_ScalingTypeSet(GeometricField,CMFE_FIELD_NO_SCALING,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,1
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Default values to 0
  X=0.0_CMISSRP
  Y=0.0_CMISSRP
  Z=0.0_CMISSRP
  DO nodeIdx=1,NUMBER_OF_NODES_SPACE
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      !(field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value,err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
    ENDIF
  ENDDO

  !parent vessel
  X=-0.5_CMISSRP
  Y=0.0_CMISSRP
  Z=0.0_CMISSRP
  DO nodeIdx=1,5
    X=X+0.5_CMISSRP
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
    ENDIF
  ENDDO
  ! node 5 versions 2,3 (same x/y values)
  nodeIdx=5
  DO versionIdx=2,3
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
    ENDIF
  ENDDO

  DO nodeIdx=6,7
    X=X+0.5_CMISSRP
    Y=Y+0.5_CMISSRP
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
    ENDIF
  ENDDO

  DO nodeIdx=10,11
    X=X+0.5_CMISSRP
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
    ENDIF
  ENDDO

  X=2.0_CMISSRP
  Y=0.0_CMISSRP
  DO nodeIdx=8,9
    X=X+0.5_CMISSRP
    Y=Y-0.5_CMISSRP
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
    ENDIF
  ENDDO

  DO nodeIdx=12,13
    X=X+0.5_CMISSRP
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,1,Y,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,2,X,Err)
    ENDIF
  ENDDO
    
  !================================================================================================================================
  !Equations Sets
  !================================================================================================================================

  !Create the equations set for 1D Navier-Stokes
  CALL cmfe_EquationsSet_Initialise(EquationsSetNavierStokes,Err)
  CALL cmfe_Field_Initialise(EquationsSetFieldNavierStokes,Err)
  !Set the equations set to be a dynamic Navier-Stokes problem
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField, &
    & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMFE_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE,EquationsSetSubtype], &
    & EquationsSetFieldUserNumberNavierStokes,EquationsSetFieldNavierStokes,EquationsSetNavierStokes,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetNavierStokes,Err)

  !Create the equations set for Characteristic
  CALL cmfe_EquationsSet_Initialise(EquationsSetCharacteristic,Err)
  CALL cmfe_Field_Initialise(EquationsSetFieldCharacteristic,Err)
  !Set the equations set to be a static Nonlinear problem
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberCharacteristic,Region,GeometricField, &
    & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMFE_EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE, &
    & EquationsSetCharacteristicSubtype],EquationsSetFieldUserNumberCharacteristic,EquationsSetFieldCharacteristic, &
    & EquationsSetCharacteristic,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetCharacteristic,Err)


  !================================================================================================================================
  !Dependent Field
  !================================================================================================================================

  !Create the equations set dependent field variables for static nonlinear Characteristic solver
  CALL cmfe_Field_Initialise(DependentFieldNavierStokes,Err)

  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetCharacteristic,DependentFieldUserNumber, & 
    & DependentFieldNavierStokes,Err)
  !Set the field label
  CALL cmfe_Field_VariableLabelSet(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,'General',Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldNavierStokes,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,'Derivatives',Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldNavierStokes,CMFE_FIELD_V_VARIABLE_TYPE,'Characteristics',Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldNavierStokes,CMFE_FIELD_U1_VARIABLE_TYPE,'CellML Pressure',Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldNavierStokes,CMFE_FIELD_U2_VARIABLE_TYPE,'1D Solver Pressure',Err)

  !Set the mesh component to be used by the field components.
  COMPONENT_NUMBER=1 ! Velocity
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
    & MESH_COMPONENT_NUMBER_VELOCITY,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
    & MESH_COMPONENT_NUMBER_VELOCITY,Err)

  COMPONENT_NUMBER=2 !Area
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
    & MESH_COMPONENT_NUMBER_AREA,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, &
    & MESH_COMPONENT_NUMBER_AREA,Err)

  DO COMPONENT_NUMBER=1,2 ! W (Characteristics)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_V_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO

  IF (cellmlFlag) THEN
    COMPONENT_NUMBER=1 ! pCellML
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_U1_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDIF
  COMPONENT_NUMBER=1 ! calculated pressure
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldNavierStokes,CMFE_FIELD_U2_VARIABLE_TYPE,COMPONENT_NUMBER, & 
    & MESH_COMPONENT_NUMBER_SPACE,Err)

  !Finish the equations set dependent field variables for Characteristic equations set
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetCharacteristic,Err)

  !Create the equations set dependent field variables for dynamic Navier-Stokes
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetNavierStokes,DependentFieldUserNumber, & 
    & DependentFieldNavierStokes,Err)

  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetNavierStokes,Err)

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
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    ENDIF
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
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    ENDIF
  ENDDO

  ! Branch versions at node 5
  nodeIdx=5
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
  IF(NodeDomain==ComputationalNodeNumber) THEN
    versionIdx=2
    componentIdx=1
    ! VALUE=Q2  
    ! CALL cmfe_Field_ParameterSetUpdateNode(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    !   & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    componentIdx=2
    VALUE=A2  
    CALL cmfe_Field_ParameterSetUpdateNode(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    versionIdx=3
    componentIdx=1
    ! VALUE=Q3  
    ! CALL cmfe_Field_ParameterSetUpdateNode(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    !   & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    componentIdx=2
    VALUE=A3  
    CALL cmfe_Field_ParameterSetUpdateNode(DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    !Initialize dependent characteristic field W (V_type) to 0
    nodeIdx=5
    VALUE=0.0_CMISSRP
    DO componentIdx=1,2
      DO versionIdx=1,3
        CALL cmfe_Field_ParameterSetUpdateNode(DependentFieldNavierStokes,CMFE_FIELD_V_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)     
      ENDDO
    ENDDO
  ENDIF


  !================================================================================================================================
  !Materials Field
  !================================================================================================================================

  !Create the equations set materials field variables for dynamic Navier-Stokes
  CALL cmfe_Field_Initialise(MaterialsFieldNavierStokes,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetNavierStokes,MaterialsFieldUserNumber, & 
    & MaterialsFieldNavierStokes,Err)
  !Set the field label
  CALL cmfe_Field_VariableLabelSet(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,'Materials',Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetNavierStokes,Err)

  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetCharacteristic,MaterialsFieldUserNumber, & 
    & MaterialsFieldNavierStokes,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetCharacteristic,Err)

  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberMu,MU_PARAM_NAVIER_STOKES,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberRho,RHO_PARAM_NAVIER_STOKES,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberK,K,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberAs,As,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberRe,Re,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberFr,Fr,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberSt,St,Err)

  !node-dependent values
  versionIdx=1
  DO nodeIdx=1,NUMBER_OF_NODES_SPACE
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,MaterialsFieldUserNumberE,E_PARAM_NAVIER_STOKES(nodeIdx),Err)
      CALL cmfe_Field_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,MaterialsFieldUserNumberH0,H0_PARAM_NAVIER_STOKES(nodeIdx),Err)
      CALL cmfe_Field_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,MaterialsFieldUserNumberA0,A0_PARAM(nodeIdx),Err)
    ENDIF
  ENDDO
  nodeIdx=5
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
  IF(NodeDomain==ComputationalNodeNumber) THEN
    DO versionIdx=2,3
      CALL cmfe_Field_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,MaterialsFieldUserNumberE,E_PARAM_NAVIER_STOKES(nodeIdx+1),Err)
      CALL cmfe_Field_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,MaterialsFieldUserNumberH0,H0_PARAM_NAVIER_STOKES(nodeIdx+1),Err)
      CALL cmfe_Field_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,MaterialsFieldUserNumberA0,A0_PARAM(nodeIdx+1),Err)
    ENDDO
  ENDIF

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
    pVesselWall=0.0_CMISSRP
    pExternal=0.0_CMISSRP

    ! Set values at coupled node 1
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,coupledNodeNumber1,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,&
       & CMFE_NO_GLOBAL_DERIV,coupledNodeNumber1,pVesselWallComponent,pVesselWall,err)
      CALL cmfe_Field_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,&
       & CMFE_NO_GLOBAL_DERIV,coupledNodeNumber1,pExternalComponent,pExternal,err)
    ENDIF

    ! Set values at coupled node 2
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,coupledNodeNumber2,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,&
       & CMFE_NO_GLOBAL_DERIV,coupledNodeNumber2,pVesselWallComponent,pVesselWall,err)
      CALL cmfe_Field_ParameterSetUpdateNode(MaterialsFieldNavierStokes,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,&
       & CMFE_NO_GLOBAL_DERIV,coupledNodeNumber2,pExternalComponent,pExternal,err)
    ENDIF
  ENDIF 


  !================================================================================================================================x
  ! Independent Field - Characteristic wave normal direction
  !================================================================================================================================x

  !Create the equations set independent field variables for Characteristics Solver
  CALL cmfe_Field_Initialise(IndependentFieldNavierStokes,Err)

  CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetCharacteristic,IndependentFieldUserNumber, & 
    & IndependentFieldNavierStokes,Err)
  CALL cmfe_Field_VariableLabelSet(IndependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,'Normal wave direction',Err)
  !Set the mesh component to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(IndependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,1, & 
    & MESH_COMPONENT_NUMBER_SPACE,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(IndependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,2, & 
    & MESH_COMPONENT_NUMBER_SPACE,Err)
  !Finish the equations set independent field variables
  CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetCharacteristic,Err)

  !(field,variableType,fieldSetType,versionNumber,derivativeNumber,userNodeNumber,componentNumber,value,err)

  ! normalWave node 5
  nodeIdx=5
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
  IF(NodeDomain==ComputationalNodeNumber) THEN
    componentIdx=1 ! Incoming normal
    ! 1 inlet/parent from element 2
    versionIdx=1
    VALUE=1.0_CMISSRP
    CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
       & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)

    ! 2 outlet/daughters from elements 3 and 4
    componentIdx=2 ! outgoing
    VALUE=-1.0_CMISSRP
    DO versionIdx=2,3
      CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    ENDDO
  ENDIF

  IF (cellmlFlag) THEN

    ! Incoming normals for 11,13
    nodeIdx=11
    componentIdx=1
    VALUE=1.0_CMISSRP
    versionIdx=1
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN

      CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE, &
         & CMFE_FIELD_VALUES_SET_TYPE,versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    ENDIF
    nodeIdx=13
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    ENDIF

    ! Outgoing normals for 11,13
    nodeIdx=11
    componentIdx=2
    VALUE=-1.0_CMISSRP
    versionIdx=1
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE, &
         & CMFE_FIELD_VALUES_SET_TYPE,versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    ENDIF
    nodeIdx=13
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(IndependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE, &
         & CMFE_FIELD_VALUES_SET_TYPE,versionIdx,CMFE_NO_GLOBAL_DERIV,nodeIdx,componentIdx,VALUE,Err)
    ENDIF

  ENDIF

  !Create the equations set independent field variables for 1D Transient Solver
  CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetNavierStokes,IndependentFieldUserNumber, & 
    & IndependentFieldNavierStokes,Err)
  !Finish the equations set independent field variables
  CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetNavierStokes,Err)


  !================================================================================================================================
  !  C e l l M L    M o d e l    M a p s
  !================================================================================================================================

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

    !Create the CellML environment
    CALL cmfe_CellML_Initialise(CellML,Err)
    CALL cmfe_CellML_CreateStart(CellMLUserNumber,Region,CellML,Err)

    ! --- W i n d k e s s e l   M o d e l --- !
    IF(windkesselFlag) THEN
      IF(mixedFlag) THEN
        CALL cmfe_CellML_ModelImport(CellML,"CellMLModels/Windkessel/WindkesselMain.cellml",CellMLModelIndex1,Err)
        CALL cmfe_CellML_ModelImport(CellML,"CellMLModels/Resistance/resistance.xml",CellMLModelIndex2,Err)
!        CALL cmfe_CellML_ModelImport(CellML,"./CellMLModels/WResistance/WindkesselMain.cellml",CellMLModelIndex2,Err)
        numberOfCellmlModels = 2
      ELSE
        CALL cmfe_CellML_ModelImport(CellML,"./CellMLModels/Windkessel/WindkesselMain.cellml",CellMLModelIndex,Err)
        numberOfCellmlModels = 1
      ENDIF
    ! --- R e s i s t a n c e   M o d e l --- !
    ELSE
      CALL cmfe_CellML_ModelImport(CellML,"./CellMLModels/Resistance/resistance.xml",CellMLModelIndex,Err)
      numberOfCellmlModels = 1
    ENDIF

    DO modelIndex=1,numberOfCellMLModels
      IF(mixedFlag) THEN
        IF(modelIndex == 1) THEN
          CellMLModelIndex = CellMLModelIndex1
        ELSE IF(modelIndex == 2) THEN
          CellMLModelIndex = CellMLModelIndex2
        ENDIF
      ENDIF
      ! - known (to OpenCMISS) variables 
      CALL cmfe_CellML_VariableSetAsKnown(CellML,CellMLModelIndex,"interface/FlowRate",Err)
      ! - to get from the CellML side 
      CALL cmfe_CellML_VariableSetAsWanted(CellML,CellMLModelIndex,"interface/Pressure",Err)
    ENDDO ! modelIndex

    CALL cmfe_CellML_CreateFinish(CellML,Err)
    !Start the creation of CellML <--> OpenCMISS field maps
    CALL cmfe_CellML_FieldMapsCreateStart(CellML,Err)

    DO modelIndex=1,numberOfCellMLModels
      IF(mixedFlag) THEN
        IF(modelIndex == 1) THEN
          CellMLModelIndex = CellMLModelIndex1
        ELSE IF(modelIndex == 2) THEN
          CellMLModelIndex = CellMLModelIndex2
        ENDIF
      ENDIF
      !Now we can set up the field variable component <--> CellML model variable mappings.

      !Map the OpenCMISS boundary flow rate values --> CellML
      ! Q is component 1 of the DependentField
      CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentFieldNavierStokes,CMFE_FIELD_U_VARIABLE_TYPE,1, &
        & CMFE_FIELD_VALUES_SET_TYPE,CellMLModelIndex,"interface/FlowRate",CMFE_FIELD_VALUES_SET_TYPE,Err)
      !Map the returned pressure values from CellML --> CMISS
      ! pCellML is component 2 of the Dependent field U1 variable
      CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,CellMLModelIndex,"interface/Pressure",CMFE_FIELD_VALUES_SET_TYPE, &
        & DependentFieldNavierStokes,CMFE_FIELD_U1_VARIABLE_TYPE,pCellMLComponent,CMFE_FIELD_VALUES_SET_TYPE,Err)
    ENDDO ! modelIndex

    !Finish the creation of CellML <--> OpenCMISS field maps
    CALL cmfe_CellML_FieldMapsCreateFinish(CellML,Err)

    !Create the CellML models field
    CALL cmfe_Field_Initialise(CellMLModelsField,Err)
    CALL cmfe_CellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
    CALL cmfe_CellML_ModelsFieldCreateFinish(CellML,Err)
    ! If we have multiple CellML models (1 R, 1 RCR boundary) set models field at each DOF
    IF(mixedFlag) THEN
      DO nodeIdx=1,numberOfNodes
        IF(nodeIdx == coupledNodeNumber1) THEN
          CellMLModelIndex = CellMLModelIndex1
        ELSE IF(nodeIdx == coupledNodeNumber2) THEN        
          CellMLModelIndex = CellMLModelIndex2
        ELSE
          CellMLModelIndex = 0_CMISSIntg
        ENDIF
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeIdx,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL cmfe_Field_ParameterSetUpdateNode(CellMLModelsField, CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE,&
           & 1,CMFE_NO_GLOBAL_DERIV,nodeIdx,1,CellMLModelIndex,Err)
        ENDIF
      ENDDO
    ENDIF

    IF (windkesselFlag) THEN 
      !Start the creation of the CellML state field
      CALL cmfe_Field_Initialise(CellMLStateField,Err)
      CALL cmfe_CellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,Err)
      !Finish the creation of the CellML state field
      CALL cmfe_CellML_StateFieldCreateFinish(CellML,Err)
    ENDIF

    !Create the CellML parameters field
    CALL cmfe_Field_Initialise(CellMLParametersField,Err)
    CALL cmfe_CellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
    CALL cmfe_CellML_ParametersFieldCreateFinish(CellML,Err)

    !Create the CellML intermediate field --- will be the pressure value returned from CellML to be used for 
    ! recalculation of the incoming Riemann variable W(2)
    CALL cmfe_Field_Initialise(CellMLIntermediateField,Err)
    CALL cmfe_CellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
    CALL cmfe_CellML_IntermediateFieldCreateFinish(CellML,Err)

    ! Initialise pCellML (and previous pCellML coupling iteration values) values to 0 at the outlet nodes
    pCellML=0.0_CMISSRP
    pPrevious=0.0_CMISSRP

    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,coupledNodeNumber1,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(DependentFieldNavierStokes,CMFE_FIELD_U1_VARIABLE_TYPE, &
       & CMFE_FIELD_VALUES_SET_TYPE,1,CMFE_NO_GLOBAL_DERIV,coupledNodeNumber1,pCellMLComponent,pCellML,err)
    ENDIF
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,coupledNodeNumber2,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(DependentFieldNavierStokes,CMFE_FIELD_U1_VARIABLE_TYPE, &
       & CMFE_FIELD_VALUES_SET_TYPE,1,CMFE_NO_GLOBAL_DERIV,coupledNodeNumber2,pCellMLComponent,pCellML,err)
    ENDIF

    CALL cmfe_Field_ParameterSetCreate(DependentFieldNavierStokes,CMFE_FIELD_U1_VARIABLE_TYPE, &
     & CMFE_FIELD_PREVIOUS_VALUES_SET_TYPE,err)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,coupledNodeNumber1,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(DependentFieldNavierStokes,CMFE_FIELD_U1_VARIABLE_TYPE, &
       & CMFE_FIELD_PREVIOUS_VALUES_SET_TYPE,1,CMFE_NO_GLOBAL_DERIV,coupledNodeNumber1,pCellMLComponent,pPrevious,err)
    ENDIF
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,coupledNodeNumber2,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(DependentFieldNavierStokes,CMFE_FIELD_U1_VARIABLE_TYPE, &
       & CMFE_FIELD_PREVIOUS_VALUES_SET_TYPE,1,CMFE_NO_GLOBAL_DERIV,coupledNodeNumber2,pCellMLComponent,pPrevious,err)
    ENDIF
 
  ENDIF ! cellml flag

  !================================================================================================================================
  !Equations
  !================================================================================================================================

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsNavierStokes,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetNavierStokes,EquationsNavierStokes,Err)
  !Set the equations matrices sparsity type
  !CALL cmfe_Equations_SparsityTypeSet(EquationsNavierStokes,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_SparsityTypeSet(EquationsNavierStokes,CMFE_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations lumping type
  CALL cmfe_Equations_LumpingTypeSet(EquationsNavierStokes,CMFE_EQUATIONS_UNLUMPED_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsNavierStokes,EQUATIONS_NAVIER_STOKES_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetNavierStokes,Err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsCharacteristic,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetCharacteristic,EquationsCharacteristic,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(EquationsCharacteristic,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsCharacteristic,EQUATIONS_NAVIER_STOKES_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetCharacteristic,Err)


  !================================================================================================================================
  !Problems
  !================================================================================================================================

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS,CMFE_PROBLEM_NAVIER_STOKES_EQUATION_TYPE, &
    & ProblemSubtype],Problem,Err)
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


  !================================================================================================================================
  !Solvers
  !================================================================================================================================

  CALL cmfe_Problem_SolversCreateStart(Problem,Err)

  ! CellML DAE solver
  IF (windkesselFlag) THEN
    CALL cmfe_Solver_Initialise(CellMLSolver,Err)
    CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverDAEUserNumber,CellMLSolver,Err)
    CALL cmfe_Solver_DAETimeStepSet(CellMLSolver,DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT,Err)
    CALL cmfe_Solver_OutputTypeSet(CellMLSolver,CMFE_SOLVER_NO_OUTPUT,Err)
  ENDIF

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(DynamicSolverNavierStokes,Err)
  CALL cmfe_Solver_Initialise(NonlinearSolverNavierStokes,Err)
  CALL cmfe_Solver_Initialise(LinearSolverNavierStokes,Err)

  CALL cmfe_Solver_Initialise(NonlinearSolverCharacteristic,Err)
  CALL cmfe_Solver_Initialise(LinearSolverCharacteristic,Err)

  !Get the dynamic dynamic solver
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
!  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
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

  !Get the static Nonlinear solver
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverCharacteristicUserNumber,NonlinearSolverCharacteristic,Err)
  !Set the nonlinear Jacobian type
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(NonlinearSolverCharacteristic, &
    & CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
!  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(NonlinearSolverCharacteristic,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(NonlinearSolverCharacteristic,NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(NonlinearSolverCharacteristic,ABSOLUTE_TOLERANCE,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(NonlinearSolverCharacteristic,RELATIVE_TOLERANCE,Err)
  !Get the nonlinear linear solver
  CALL cmfe_Solver_NewtonLinearSolverGet(NonlinearSolverCharacteristic,LinearSolverCharacteristic,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(LinearSolverCharacteristic,LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolverCharacteristic,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LibraryTypeSet(LinearSolverCharacteristic,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolverCharacteristic,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolverCharacteristic,MAXIMUM_ITERATIONS,Err)
    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolverCharacteristic,DIVERGENCE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolverCharacteristic,RELATIVE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolverCharacteristic,ABSOLUTE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(LinearSolverCharacteristic,RESTART_VALUE,Err)
  ENDIF

  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)


  ! C e l l M L   S o l v e r
  !-----------------------------
  IF (cellmlFlag) THEN

    !Create the problem solver CellML equations
    CALL cmfe_Solver_Initialise(CellMLSolver,Err)
    CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
    CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)

    IF(windkesselFlag) THEN
      !Get the DAE solver  
      CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverDAEUserNumber,CellMLSolver,Err)
    ELSE
      CALL cmfe_Solver_NewtonCellMLSolverGet(DynamicSolverNavierStokes,CellMLSolver,Err) 
    ENDIF
    CALL cmfe_Solver_CellMLEquationsGet(CellMLSolver,CellMLEquations,Err)
    CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
    CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)

  ENDIF


  !================================================================================================================================
  !Solver Equations
  !================================================================================================================================

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(DynamicSolverNavierStokes,Err)
  CALL cmfe_Solver_Initialise(NonlinearSolverCharacteristic,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsNavierStokes,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsCharacteristic,Err)

  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)

  !Get the dynamic Navier-Stokes solver equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  CALL cmfe_Solver_SolverEquationsGet(DynamicSolverNavierStokes,SolverEquationsNavierStokes,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsNavierStokes,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsNavierStokes,EquationsSetNavierStokes,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations

  !Get the static Nonlinear solver equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverCharacteristicUserNumber,NonlinearSolverCharacteristic,Err)
  CALL cmfe_Solver_SolverEquationsGet(NonlinearSolverCharacteristic,SolverEquationsCharacteristic,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsCharacteristic,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsCharacteristic,EquationsSetCharacteristic,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations

  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)


  !================================================================================================================================
  !Boundary Conditions
  !================================================================================================================================

  !Start the creation of the equations set boundary conditions for Navier-Stokes
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsNavierStokes,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsNavierStokes,BoundaryConditionsNavierStokes,Err)
  !Set velocity boundary conditions
  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    CONDITION=CMFE_BOUNDARY_CONDITION_FIXED_INLET
    VALUE1=BOUNDARY_CONDITIONS_NAVIER_STOKES(1)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,1,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV,1,1,CONDITION,VALUE1,Err)
    ENDIF
  ENDIF
  !Set area boundary conditions if not a CellML problem
  IF(cellmlFlag) THEN
    CONDITION=CMFE_BOUNDARY_CONDITION_FIXED_OUTLET
    VALUE1=A2
    versionIdx=1
    !(boundaryConditions,field,variableType,versionNumber,derivativeNumber,nodeUserNumber,componentNumber,condition,value,err)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,coupledNodeNumber1,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
        & CMFE_FIELD_U_VARIABLE_TYPE,versionIdx,CMFE_NO_GLOBAL_DERIV,coupledNodeNumber1,2,CONDITION,VALUE1,Err)
    ENDIF
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,coupledNodeNumber2,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
        & CMFE_FIELD_U_VARIABLE_TYPE,versionIdx,CMFE_NO_GLOBAL_DERIV,coupledNodeNumber2,2,CONDITION,VALUE1,Err)
    ENDIF
  ELSE
    CONDITION=CMFE_BOUNDARY_CONDITION_FIXED
    VALUE1=A2    
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,coupledNodeNumber1,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV,coupledNodeNumber1,2,CONDITION,VALUE1,Err)
    ENDIF
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,coupledNodeNumber1,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV,coupledNodeNumber2,2,CONDITION,VALUE1,Err)
    ENDIF
  ENDIF

  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsCharacteristic,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsCharacteristic,BoundaryConditionsCharacteristic,Err)
  IF(cellmlFlag) THEN
    CONDITION=CMFE_BOUNDARY_CONDITION_FIXED_OUTLET
    VALUE1=A2
    versionIdx=1
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,coupledNodeNumber1,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsCharacteristic,DependentFieldNavierStokes, &
        & CMFE_FIELD_U_VARIABLE_TYPE,versionIdx,CMFE_NO_GLOBAL_DERIV,coupledNodeNumber1,2,CONDITION,VALUE1,Err)
    ENDIF
    VALUE1=A3
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,coupledNodeNumber2,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsCharacteristic,DependentFieldNavierStokes, &
        & CMFE_FIELD_U_VARIABLE_TYPE,versionIdx,CMFE_NO_GLOBAL_DERIV,coupledNodeNumber2,2,CONDITION,VALUE1,Err)
    ENDIF
  ENDIF

  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsNavierStokes,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsCharacteristic,Err)


  ! Note: CellML Parameters (e.g. resistance, capacitance) should be set within each CellML model file
  !================================================================================================================================
  ! RUN SOLVERS
  !================================================================================================================================

  !Turn of PETSc error handling
!  CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL cmfe_Problem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"

  !
  !================================================================================================================================
  !

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM FortranExample
