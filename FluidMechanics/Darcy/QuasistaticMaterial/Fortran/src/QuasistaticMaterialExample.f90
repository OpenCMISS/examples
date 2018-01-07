!> \file
!> \author Christian Michler
!> \brief This is an example program to solve a Darcy equation using OpenCMISS calls.
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

!> \example FluidMechanics/Darcy/QuasistaticMaterial/src/QuasistaticMaterialExample.f90
!! Example program to solve a QuasistaticMaterial ALE Darcy equation using OpenCMISS calls.
!!
!! !! \htmlinclude FluidMechanics/Darcy/QuasistaticMaterial/history.html
!<

! ! 
! !  This example considers a moving mesh, time-dependent boundary conditions and
! !  mesh-deformation dependent material properties and, thus, we solve it as a quasi-static ALE Darcy problem.
! ! 

!> Main program

PROGRAM DARCYQUASISTATICMATERIALEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE FLUID_MECHANICS_IO_ROUTINES
#ifndef NOMPIMOD
  USE MPI
#endif

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
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberDarcy=6
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberMatProperties=42
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcy=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMatProperties=9
!   INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberDarcy=10
!   INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberMatProperties=11
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberDarcy=12
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberMatProperties=13
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberDarcy=21
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberMatProperties=22

  INTEGER(CMISSIntg), PARAMETER :: SolverMatPropertiesUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverDarcyUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPorosity=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPermOverVis=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMatPropertiesPorosity=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMatPropertiesPermOverVis=2

  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_PRESSURE
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_GEOMETRY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_GEOMETRY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_GEOMETRY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_PRESSURE
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
!   INTEGER(CMISSIntg) :: MPI_IERROR

  INTEGER(CMISSIntg) :: EQUATIONS_DARCY_OUTPUT
  INTEGER(CMISSIntg) :: EQUATIONS_MAT_PROPERTIES_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: LINEAR_SOLVER_DARCY_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_DARCY_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_MAT_PROPERTIES_OUTPUT_TYPE

  REAL(CMISSRP) :: COORD_X, COORD_Y, COORD_Z
  REAL(CMISSRP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2
  REAL(CMISSRP) :: GEOMETRY_TOLERANCE

  REAL(CMISSRP) :: INITIAL_FIELD_DARCY(3)
  REAL(CMISSRP) :: INITIAL_FIELD_MAT_PROPERTIES(3)
  REAL(CMISSRP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSRP) :: RELATIVE_TOLERANCE
  REAL(CMISSRP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSRP) :: LINESEARCH_ALPHA
  REAL(CMISSRP) :: VALUE
  REAL(CMISSRP) :: POROSITY_PARAM_MAT_PROPERTIES, PERM_OVER_VIS_PARAM_MAT_PROPERTIES
  REAL(CMISSRP) :: POROSITY_PARAM_DARCY, PERM_OVER_VIS_PARAM_DARCY

  REAL(CMISSRP) :: LINEAR_SOLVER_DARCY_START_TIME
  REAL(CMISSRP) :: LINEAR_SOLVER_DARCY_STOP_TIME
  REAL(CMISSRP) :: LINEAR_SOLVER_DARCY_TIME_INCREMENT

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_DARCY_DIRECT_FLAG
  LOGICAL :: LINEAR_SOLVER_MAT_PROPERTIES_DIRECT_FLAG

  !CMISS variables

  !Regions
  TYPE(cmfe_RegionType) :: Region
  TYPE(cmfe_RegionType) :: WorldRegion
  !Coordinate systems
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
  TYPE(cmfe_CoordinateSystemType) :: WorldCoordinateSystem
  !Basis
  TYPE(cmfe_BasisType) :: BasisGeometry
  TYPE(cmfe_BasisType) :: BasisVelocity
  TYPE(cmfe_BasisType) :: BasisPressure
  !Nodes
  TYPE(cmfe_NodesType) :: Nodes
  !Elements
  TYPE(cmfe_MeshElementsType) :: MeshElementsGeometry
  TYPE(cmfe_MeshElementsType) :: MeshElementsVelocity
  TYPE(cmfe_MeshElementsType) :: MeshElementsPressure
  !Meshes
  TYPE(cmfe_MeshType) :: Mesh
  !Decompositions
  TYPE(cmfe_DecompositionType) :: Decomposition
  !Fields
  TYPE(cmfe_FieldsType) :: Fields
  !Field types
  TYPE(cmfe_FieldType) :: GeometricField
  TYPE(cmfe_FieldType) :: DependentFieldDarcy
  TYPE(cmfe_FieldType) :: DependentFieldMatProperties
  TYPE(cmfe_FieldType) :: MaterialsFieldDarcy
  TYPE(cmfe_FieldType) :: MaterialsFieldMatProperties
  TYPE(cmfe_FieldType) :: EquationsSetFieldDarcy
  TYPE(cmfe_FieldType) :: EquationsSetFieldMatProperties  
!   TYPE(cmfe_FieldType) :: IndependentFieldDarcy
!   TYPE(cmfe_FieldType) :: IndependentFieldMatProperties
  !Boundary conditions
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsDarcy
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsMatProperties 
  !Equations sets
  TYPE(cmfe_EquationsSetType) :: EquationsSetDarcy
  TYPE(cmfe_EquationsSetType) :: EquationsSetMatProperties
  !Equations
  TYPE(cmfe_EquationsType) :: EquationsDarcy
  TYPE(cmfe_EquationsType) :: EquationsMatProperties
  !Problems
  TYPE(cmfe_ProblemType) :: Problem
  !Control loops
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  !Solvers
  TYPE(cmfe_SolverType) :: LinearSolverDarcy
  TYPE(cmfe_SolverType) :: LinearSolverMatProperties
  !Solver equations
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsDarcy
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsMatProperties
  !Computational information
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,NodeDomain

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: Err


  INTEGER(CMISSIntg) :: DIAG_LEVEL_LIST(5)
  CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(6) !,TIMING_ROUTINE_LIST(1)

  
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

  !INITIALISE OPENCMISS

  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !
  !================================================================================================================================
  !

  !PROBLEM CONTROL PANEL

  !Import cmHeart mesh information
  CALL FLUID_MECHANICS_IO_READ_CMHEART(CM,Err)  
  BASIS_NUMBER_GEOMETRY=CM%ID_M
  BASIS_NUMBER_VELOCITY=CM%ID_V
  BASIS_NUMBER_PRESSURE=CM%ID_P
  NUMBER_OF_DIMENSIONS=CM%D
  BASIS_TYPE=CM%IT_T
  BASIS_XI_INTERPOLATION_GEOMETRY=CM%IT_M
  BASIS_XI_INTERPOLATION_VELOCITY=CM%IT_V
  BASIS_XI_INTERPOLATION_PRESSURE=CM%IT_P
  NUMBER_OF_NODES_GEOMETRY=CM%N_M
  NUMBER_OF_NODES_VELOCITY=CM%N_V
  NUMBER_OF_NODES_PRESSURE=CM%N_P
  TOTAL_NUMBER_OF_NODES=CM%N_T
  TOTAL_NUMBER_OF_ELEMENTS=CM%E_T
  NUMBER_OF_ELEMENT_NODES_GEOMETRY=CM%EN_M
  NUMBER_OF_ELEMENT_NODES_VELOCITY=CM%EN_V
  NUMBER_OF_ELEMENT_NODES_PRESSURE=CM%EN_P
  !Set domain dimensions
  DOMAIN_X1 = -5.0_CMISSRP
  DOMAIN_X2 =  5.0_CMISSRP
  DOMAIN_Y1 = -5.0_CMISSRP
  DOMAIN_Y2 =  5.0_CMISSRP
  DOMAIN_Z1 = -5.0_CMISSRP
  DOMAIN_Z2 =  5.0_CMISSRP
  !Set geometric tolerance
  GEOMETRY_TOLERANCE = 1.0E-12_CMISSRP
  !Set initial values
  INITIAL_FIELD_DARCY(1)=0.0_CMISSRP
  INITIAL_FIELD_DARCY(2)=0.0_CMISSRP
  INITIAL_FIELD_DARCY(3)=0.0_CMISSRP
  INITIAL_FIELD_MAT_PROPERTIES(1)=0.0_CMISSRP
  INITIAL_FIELD_MAT_PROPERTIES(2)=0.0_CMISSRP
  INITIAL_FIELD_MAT_PROPERTIES(3)=0.0_CMISSRP
  !Set material parameters
  POROSITY_PARAM_DARCY=0.3_CMISSRP
  PERM_OVER_VIS_PARAM_DARCY=1.0_CMISSRP
  POROSITY_PARAM_MAT_PROPERTIES=POROSITY_PARAM_DARCY
  PERM_OVER_VIS_PARAM_MAT_PROPERTIES=PERM_OVER_VIS_PARAM_DARCY
  !Set number of Gauss points (Mind that also material field may be interpolated)
  BASIS_XI_GAUSS_GEOMETRY=3 !4
  BASIS_XI_GAUSS_VELOCITY=3 !4
  BASIS_XI_GAUSS_PRESSURE=3 !4
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_MAT_PROPERTIES_OUTPUT_TYPE=CMFE_SOLVER_NO_OUTPUT
  LINEAR_SOLVER_DARCY_OUTPUT_TYPE=CMFE_SOLVER_PROGRESS_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DARCY_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT
  EQUATIONS_MAT_PROPERTIES_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT
  !Set time parameter
  LINEAR_SOLVER_DARCY_START_TIME=0.0_CMISSRP
  LINEAR_SOLVER_DARCY_STOP_TIME=1.125_CMISSRP
  LINEAR_SOLVER_DARCY_TIME_INCREMENT=0.125_CMISSRP
  !Set result output parameter
  LINEAR_SOLVER_DARCY_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_MAT_PROPERTIES_DIRECT_FLAG=.FALSE.
  LINEAR_SOLVER_DARCY_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-05_CMISSRP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-10_CMISSRP
  DIVERGENCE_TOLERANCE=1.0E5_CMISSRP !default: 1.0E5
  MAXIMUM_ITERATIONS=10000_CMISSIntg !default: 100000
  RESTART_VALUE=3000_CMISSIntg !default: 30
  LINESEARCH_ALPHA=1.0_CMISSRP


  !
  !================================================================================================================================
  !

  !Set diagnostics

  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5

  DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_FINITE_ELEMENT_CALCULATE"
  DIAG_ROUTINE_LIST(2)="DARCY_EQUATION_PRE_SOLVE_STORE_REFERENCE_DATA"
  DIAG_ROUTINE_LIST(3)="DARCY_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH"
  DIAG_ROUTINE_LIST(4)="DARCY_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS"
  DIAG_ROUTINE_LIST(5)="DARCY_EQUATION_PRE_SOLVE_MAT_PROPERTIES"
  DIAG_ROUTINE_LIST(6)="FITTING_FINITE_ELEMENT_CALCULATE"

  !CMFE_ALL_DIAG_TYPE/CMFE_IN_DIAG_TYPE/CMFE_FROM_DIAG_TYPE
!   CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,DIAG_LEVEL_LIST,"DarcyDiagnostics",DIAG_ROUTINE_LIST,Err)

  !CMFE_ALL_TIMING_TYPE/CMFE_IN_TIMING_TYPE/CMFE_FROM_TIMING_TYPE
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !REGION

  !Start the creation of a new region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system as defined above
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !BASES

  !Start the creation of new bases
  MESH_NUMBER_OF_COMPONENTS=1
  CALL cmfe_Basis_Initialise(BasisGeometry,Err)
  CALL cmfe_Basis_CreateStart(BASIS_NUMBER_GEOMETRY,BasisGeometry,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL cmfe_Basis_TypeSet(BasisGeometry,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL cmfe_Basis_NumberOfXiSet(BasisGeometry,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  IF(NUMBER_OF_DIMENSIONS==2) THEN
    CALL cmfe_Basis_InterpolationXiSet(BasisGeometry,[BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisGeometry,[BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY],Err)
  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
    CALL cmfe_Basis_InterpolationXiSet(BasisGeometry,[BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY, & 
      & BASIS_XI_INTERPOLATION_GEOMETRY],Err)                         
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisGeometry,[BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY, &
      & BASIS_XI_GAUSS_GEOMETRY],Err)
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(BasisGeometry,Err)
  !
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisVelocity=BasisGeometry
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
  !
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisPressure=BasisGeometry
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

  !
  !================================================================================================================================
  !

  !MESH

  !Start the creation of mesh nodes
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,TOTAL_NUMBER_OF_NODES,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)
  !Specify spatial mesh component
  CALL cmfe_MeshElements_Initialise(MeshElementsGeometry,Err)
  CALL cmfe_MeshElements_Initialise(MeshElementsVelocity,Err)
  CALL cmfe_MeshElements_Initialise(MeshElementsPressure,Err)
  MESH_COMPONENT_NUMBER_GEOMETRY=1
  MESH_COMPONENT_NUMBER_VELOCITY=1
  MESH_COMPONENT_NUMBER_PRESSURE=1
  CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_GEOMETRY,BasisGeometry,MeshElementsGeometry,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL cmfe_MeshElements_NodesSet(MeshElementsGeometry,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(MeshElementsGeometry,Err)
  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    MeshElementsVelocity=MeshElementsGeometry
  ELSE
    MESH_COMPONENT_NUMBER_VELOCITY=MESH_COMPONENT_NUMBER_GEOMETRY+1
    CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_VELOCITY,BasisVelocity,MeshElementsVelocity,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL cmfe_MeshElements_NodesSet(MeshElementsVelocity,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_VELOCITY),Err)
    ENDDO
    CALL cmfe_MeshElements_CreateFinish(MeshElementsVelocity,Err)
  ENDIF
  !Specify pressure mesh component
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    MeshElementsPressure=MeshElementsGeometry
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_GEOMETRY
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

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

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
      & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
  ENDDO

  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
        & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for ALE Darcy
  CALL cmfe_Field_Initialise(EquationsSetFieldDarcy,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetDarcy,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberDarcy,Region,GeometricField,[CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & CMFE_EQUATIONS_SET_DARCY_EQUATION_TYPE,CMFE_EQUATIONS_SET_ALE_DARCY_SUBTYPE],EquationsSetFieldUserNumberDarcy, &
    & EquationsSetFieldDarcy,EquationsSetDarcy,Err)
  !Set the equations set to be a ALE Darcy problem
!   CALL cmfe_EquationsSet_SpecificationSet(EquationsSetDarcy,CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
!     & CMFE_EQUATIONS_SET_DARCY_EQUATION_TYPE,CMFE_EQUATIONS_SET_ALE_DARCY_SUBTYPE,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetDarcy,Err)

  !Create the equations set for deformation-dependent material properties
  CALL cmfe_Field_Initialise(EquationsSetFieldMatProperties,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetMatProperties,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberMatProperties,Region,GeometricField,[CMFE_EQUATIONS_SET_FITTING_CLASS, &
    & CMFE_EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE,CMFE_EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE], &
    & EquationsSetFieldUserNumberMatProperties,EquationsSetFieldMatProperties,EquationsSetMatProperties,Err)
  !Set the equations set to be a deformation-dependent material properties problem
!   CALL cmfe_EquationsSet_SpecificationSet(EquationsSetMatProperties,CMFE_EQUATIONS_SET_FITTING_CLASS, &
!     & CMFE_EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE,CMFE_EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetMatProperties,Err)

  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for ALE Darcy
  CALL cmfe_Field_Initialise(DependentFieldDarcy,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetDarcy,DependentFieldUserNumberDarcy, & 
    & DependentFieldDarcy,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDarcy,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
  ENDDO
  COMPONENT_NUMBER=NUMBER_OF_DIMENSIONS+1
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDarcy,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetDarcy,Err)

  !Initialise dependent field (velocity components)
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
  ENDDO

  !Create the equations set dependent field variables for deformation-dependent material properties
  CALL cmfe_Field_Initialise(DependentFieldMatProperties,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetMatProperties,DependentFieldUserNumberMatProperties, & 
    & DependentFieldMatProperties,Err)
  !Set the mesh component to be used by the field components.
  NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES = 2
  DO COMPONENT_NUMBER=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldMatProperties,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldMatProperties,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
  ENDDO
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetMatProperties,Err)

  !Initialise dependent field
  DO COMPONENT_NUMBER=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldMatProperties,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & COMPONENT_NUMBER,INITIAL_FIELD_MAT_PROPERTIES(COMPONENT_NUMBER),Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for ALE Darcy
  CALL cmfe_Field_Initialise(MaterialsFieldDarcy,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetDarcy,MaterialsFieldUserNumberDarcy, & 
    & MaterialsFieldDarcy,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetDarcy,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)
  !Create the equations set materials field variables for deformation-dependent material properties
  CALL cmfe_Field_Initialise(MaterialsFieldMatProperties,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetMatProperties,MaterialsFieldUserNumberMatProperties, & 
    & MaterialsFieldMatProperties,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetMatProperties,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldMatProperties,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberMatPropertiesPorosity,POROSITY_PARAM_MAT_PROPERTIES,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldMatProperties,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberMatPropertiesPermOverVis,PERM_OVER_VIS_PARAM_MAT_PROPERTIES,Err)

  !
  !================================================================================================================================
  !

!   !INDEPENDENT FIELDS
! 
!   !Create the equations set independent field variables for ALE Darcy
!   CALL cmfe_Field_Initialise(IndependentFieldDarcy,Err)
!   CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetDarcy,IndependentFieldUserNumberDarcy, & 
!     & IndependentFieldDarcy,Err)
!   !Set the mesh component to be used by the field components.
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS+1
!     CALL cmfe_Field_ComponentMeshComponentSet(IndependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
!   ENDDO
!   !Finish the equations set independent field variables
!   CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetDarcy,Err)
!   !Create the equations set independent field variables for deformation-dependent material properties
!   CALL cmfe_Field_Initialise(IndependentFieldMatProperties,Err)
!   CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetMatProperties,IndependentFieldUserNumberMatProperties, & 
!     & IndependentFieldMatProperties,Err)
!   !Set the mesh component to be used by the field components.
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!     CALL cmfe_Field_ComponentMeshComponentSet(IndependentFieldMatProperties,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
!   ENDDO
!   !Finish the equations set independent field variables
!   CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetMatProperties,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsDarcy,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetDarcy,EquationsDarcy,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(EquationsDarcy,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
!   !Set the equations lumping type
!   CALL cmfe_Equations_LumpingTypeSet(EquationsDarcy,CMFE_EQUATIONS_UNLUMPED_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsDarcy,EQUATIONS_DARCY_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetDarcy,Err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsMatProperties,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetMatProperties,EquationsMatProperties,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(EquationsMatProperties,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsMatProperties,EQUATIONS_MAT_PROPERTIES_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetMatProperties,Err)

  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS,CMFE_PROBLEM_DARCY_EQUATION_TYPE, &
    & CMFE_PROBLEM_ALE_DARCY_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,LINEAR_SOLVER_DARCY_START_TIME,LINEAR_SOLVER_DARCY_STOP_TIME, & 
    & LINEAR_SOLVER_DARCY_TIME_INCREMENT,Err)
  !Set the output timing
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoop,CMFE_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,LINEAR_SOLVER_DARCY_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(LinearSolverMatProperties,Err)
  CALL cmfe_Solver_Initialise(LinearSolverDarcy,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  !Get the deformation-dependent material properties solver
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverMatPropertiesUserNumber,LinearSolverMatProperties,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(LinearSolverMatProperties,LINEAR_SOLVER_MAT_PROPERTIES_OUTPUT_TYPE,Err)
!   CALL cmfe_Solver_OutputTypeSet(LinearSolverMatProperties,4,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_MAT_PROPERTIES_DIRECT_FLAG) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolverMatProperties,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LibraryTypeSet(LinearSolverMatProperties,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolverMatProperties,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolverMatProperties,MAXIMUM_ITERATIONS,Err)
    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolverMatProperties,DIVERGENCE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolverMatProperties,RELATIVE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolverMatProperties,ABSOLUTE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(LinearSolverMatProperties,RESTART_VALUE,Err)
  ENDIF
  !Get the Darcy solver
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverDarcyUserNumber,LinearSolverDarcy,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(LinearSolverDarcy,LINEAR_SOLVER_DARCY_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_DARCY_DIRECT_FLAG) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolverDarcy,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LibraryTypeSet(LinearSolverDarcy,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolverDarcy,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolverDarcy,MAXIMUM_ITERATIONS,Err)
    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolverDarcy,DIVERGENCE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolverDarcy,RELATIVE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolverDarcy,ABSOLUTE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(LinearSolverDarcy,RESTART_VALUE,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(LinearSolverMatProperties,Err)
  CALL cmfe_Solver_Initialise(LinearSolverDarcy,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsMatProperties,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsDarcy,Err)

  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the deformation-dependent material properties solver equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverMatPropertiesUserNumber,LinearSolverMatProperties,Err)
  CALL cmfe_Solver_SolverEquationsGet(LinearSolverMatProperties,SolverEquationsMatProperties,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsMatProperties,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsMatProperties,EquationsSetMatProperties,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  !Get the Darcy solver equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverDarcyUserNumber,LinearSolverDarcy,Err)
  CALL cmfe_Solver_SolverEquationsGet(LinearSolverDarcy,SolverEquationsDarcy,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsDarcy,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsDarcy,EquationsSetDarcy,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS
  !Start the creation of the equations set boundary conditions for Darcy
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsDarcy,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)

  !--- BCs on normal velocity only
  CONDITION = CMFE_BOUNDARY_CONDITION_MOVED_WALL

  IF( CM%D==2_CMISSIntg ) THEN
    DO NODE_NUMBER=1_CMISSIntg,NUMBER_OF_NODES_GEOMETRY
      COORD_X = CM%N(NODE_NUMBER,1_CMISSIntg)
      COORD_Y = CM%N(NODE_NUMBER,2_CMISSIntg)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,MESH_COMPONENT_NUMBER_VELOCITY,NodeDomain,Err)

      IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) .AND. ComputationalNodeNumber==NodeDomain ) THEN
        !x-velocity
        VALUE = 1.0_CMISSRP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_NO_GLOBAL_DERIV,1,NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
      END IF
      !
      IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) .AND. ComputationalNodeNumber==NodeDomain ) THEN
        !x-velocity
        VALUE = 1.0_CMISSRP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_NO_GLOBAL_DERIV,1,NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
      END IF
      !
      IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) .AND. ComputationalNodeNumber==NodeDomain ) THEN
        !y-velocity
        VALUE = 2.0_CMISSRP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_NO_GLOBAL_DERIV,1,NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
      END IF
      !
      IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) .AND. ComputationalNodeNumber==NodeDomain ) THEN
        !y-velocity
        VALUE = 2.0_CMISSRP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_NO_GLOBAL_DERIV,1,NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
      END IF
    END DO
  ELSE IF( CM%D==3_CMISSIntg ) THEN
    DO NODE_NUMBER=1_CMISSIntg,NUMBER_OF_NODES_GEOMETRY
      COORD_X = CM%N(NODE_NUMBER,1_CMISSIntg)
      COORD_Y = CM%N(NODE_NUMBER,2_CMISSIntg)
      COORD_Z = CM%N(NODE_NUMBER,3_CMISSIntg)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,MESH_COMPONENT_NUMBER_VELOCITY,NodeDomain,Err)

      IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) .AND. ComputationalNodeNumber==NodeDomain ) THEN
        !x-velocity
        VALUE = 1.0_CMISSRP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_NO_GLOBAL_DERIV,1,NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
      END IF
      !
      IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) .AND. ComputationalNodeNumber==NodeDomain ) THEN
        !x-velocity
        VALUE = 1.0_CMISSRP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_NO_GLOBAL_DERIV,1,NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
      END IF
      !
      IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) .AND. ComputationalNodeNumber==NodeDomain ) THEN
        !y-velocity
        VALUE = 1.0_CMISSRP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_NO_GLOBAL_DERIV,1,NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
      END IF
      !
      IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) .AND. ComputationalNodeNumber==NodeDomain ) THEN
        !y-velocity
        VALUE = 1.0_CMISSRP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_NO_GLOBAL_DERIV,1,NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
      END IF
      !
      IF( (ABS(COORD_Z-DOMAIN_Z1) < GEOMETRY_TOLERANCE) .AND. ComputationalNodeNumber==NodeDomain ) THEN
        !z-velocity
        VALUE = 1.0_CMISSRP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_NO_GLOBAL_DERIV,1,NODE_NUMBER,3_CMISSIntg,CONDITION,VALUE,Err)
      END IF
      !
      IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) .AND. ComputationalNodeNumber==NodeDomain ) THEN
        !z-velocity
        VALUE = 1.0_CMISSRP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_NO_GLOBAL_DERIV,1,NODE_NUMBER,3_CMISSIntg,CONDITION,VALUE,Err)
      END IF
    END DO
  END IF

  !Finish the creation of the equations set boundary conditions for Darcy
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsDarcy,Err)
  !
  !Start the creation of the equations set boundary conditions for deformation-dependent material properties
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsMatProperties,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsMatProperties,BoundaryConditionsMatProperties,Err)
  !(No boundary conditions requrired for deformation-dependent material properties)
  !Finish the creation of the equations set boundary conditions for deformation-dependent material properties
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsMatProperties,Err)

  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn of PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL cmfe_Problem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"

  !
  !================================================================================================================================
  !

  !OUTPUT

  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"DarcyQuasistaticMaterial","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"DarcyQuasistaticMaterial","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF

  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM DARCYQUASISTATICMATERIALEXAMPLE
