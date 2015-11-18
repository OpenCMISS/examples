!> \file
!> \author Christian Michler
!> \brief This is an example program to solve a static Darcy equation using OpenCMISS calls.
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

!> \example FluidMechanics/Darcy/Static/src/StaticExample.f90
!! Example program to solve a static Darcy equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FluidMechanics/Darcy/Static/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FluidMechanics/Darcy/Static/build-intel'>Linux GNU Build</a>
!!
!<

! ! 
! !  This example considers a static Darcy problem.
! ! 

!> Main program

PROGRAM DARCYSTATICEXAMPLE

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
  INTEGER(CMISSIntg) :: MaterialsFieldUserNumberDarcy
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=9

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverDarcyUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPorosity=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPermOverVis=2


  REAL(CMISSRP), PARAMETER :: PI=3.141592653589793238462643383279502884197_CMISSRP


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
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: LINEAR_SOLVER_DARCY_OUTPUT_TYPE

  REAL(CMISSRP) :: COORD_X, COORD_Y, COORD_Z
  REAL(CMISSRP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2
  REAL(CMISSRP) :: GEOMETRY_TOLERANCE

  REAL(CMISSRP) :: INITIAL_FIELD_DARCY(3)
  REAL(CMISSRP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSRP) :: RELATIVE_TOLERANCE
  REAL(CMISSRP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSRP) :: VALUE
  REAL(CMISSRP) :: POROSITY_PARAM_DARCY, PERM_OVER_VIS_PARAM_DARCY

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_DARCY_DIRECT_FLAG

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
  TYPE(cmfe_FieldType), ALLOCATABLE :: MaterialsFieldDarcy(:)
  !Boundary conditions
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsDarcy
  !Equations sets
!   TYPE(cmfe_EquationsSetType) :: EquationsSetDarcy
  TYPE(cmfe_EquationsSetType), ALLOCATABLE :: EquationsSetDarcy(:)
  TYPE(cmfe_EquationsSetType) :: EquationsSetMatProperties
  !Equations
  TYPE(cmfe_EquationsType), ALLOCATABLE :: EquationsDarcy(:)
  !Problems
  TYPE(cmfe_ProblemType) :: Problem
  !Control loops
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  !Solvers
  TYPE(cmfe_SolverType) :: LinearSolverDarcy
  !Solver equations
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsDarcy

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg), ALLOCATABLE :: EquationsSetIndices(:)
  INTEGER(CMISSIntg) :: Err
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:) :: VariableTypes


  INTEGER(CMISSIntg) :: DIAG_LEVEL_LIST(5)
  CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(1) !,TIMING_ROUTINE_LIST(1)


  TYPE(cmfe_FieldType), ALLOCATABLE :: EquationsSetFieldDarcy(:)
  INTEGER(CMISSIntg) :: EquationsSetFieldDarcyUserNumber
  INTEGER(CMISSIntg) :: EquationsSetUserNumberDarcy
  INTEGER(CMISSIntg) :: icomp,Ncompartments,num_var,icompartment

  
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
  !Set material parameters
  POROSITY_PARAM_DARCY=0.3_CMISSRP
  PERM_OVER_VIS_PARAM_DARCY=0.8_CMISSRP 
  !Set number of Gauss points (Mind that also material field may be interpolated)
  BASIS_XI_GAUSS_GEOMETRY=3 !4
  BASIS_XI_GAUSS_VELOCITY=3 !4
  BASIS_XI_GAUSS_PRESSURE=3 !4
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_DARCY_OUTPUT_TYPE=CMFE_SOLVER_SOLVER_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DARCY_OUTPUT=CMFE_EQUATIONS_MATRIX_OUTPUT
  !Set solver parameters
  LINEAR_SOLVER_DARCY_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-05_CMISSRP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-10_CMISSRP
  DIVERGENCE_TOLERANCE=1.0E5_CMISSRP !default: 1.0E5
  MAXIMUM_ITERATIONS=10000_CMISSIntg !default: 100000
  RESTART_VALUE=3000_CMISSIntg !default: 30

  Ncompartments = 4_CMISSIntg

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

  !Start the creation of new bases: Geometry
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
  !Start the creation of another basis: Velocity
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
  !Start the creation of another basis: Pressure
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
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,DomainUserNumber,Err)
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

  IF(.NOT.ALLOCATED(EquationsSetFieldDarcy)) ALLOCATE(EquationsSetFieldDarcy(Ncompartments))
  IF(.NOT.ALLOCATED(EquationsSetDarcy)) ALLOCATE(EquationsSetDarcy(Ncompartments))
  IF(.NOT.ALLOCATED(MaterialsFieldDarcy)) ALLOCATE(MaterialsFieldDarcy(Ncompartments))
  IF(.NOT.ALLOCATED(EquationsSetIndices)) ALLOCATE(EquationsSetIndices(Ncompartments))
  IF(.NOT.ALLOCATED(EquationsDarcy)) ALLOCATE(EquationsDarcy(Ncompartments))


  DO icomp = 1,Ncompartments
    EquationsSetFieldDarcyUserNumber = 100_CMISSIntg + icomp
    EquationsSetUserNumberDarcy = 200_CMISSIntg + icomp

    !Create the equations set field for Darcy - *** MULTIPLE COMPARTMENTS ***
    CALL cmfe_Field_Initialise(EquationsSetFieldDarcy(icomp),Err)

!   CALL cmfe_Field_CreateStart(FieldFibreSolidUserNumber,Region,EquationsSetFieldDarcy,Err)
!   CALL cmfe_Field_TypeSet(EquationsSetFieldDarcy,CMFE_FIELD_FIBRE_TYPE,Err)
!   CALL cmfe_Field_MeshDecompositionSet(EquationsSetFieldDarcy,Decomposition,Err)        
!   CALL cmfe_Field_GeometricFieldSet(EquationsSetFieldDarcy,GeometricFieldSolid,Err)
!   CALL cmfe_Field_NumberOfVariablesSet(EquationsSetFieldDarcy,FieldFibreSolidNumberOfVariables,Err)
!   CALL cmfe_Field_NumberOfComponentsSet(EquationsSetFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,FieldFibreSolidNumberOfComponents,Err)  
!   CALL cmfe_Field_ComponentMeshComponentSet(EquationsSetFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,1,SolidMeshComponenetNumber,Err)
!   CALL cmfe_Field_ComponentMeshComponentSet(EquationsSetFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,2,SolidMeshComponenetNumber,Err)
!   CALL cmfe_Field_ComponentMeshComponentSet(EquationsSetFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,3,SolidMeshComponenetNumber,Err)
!   CALL cmfe_Field_CreateFinish(EquationsSetFieldDarcy,Err)

    !These two lines have to go into each example file:
    CALL cmfe_EquationsSet_Initialise(EquationsSetDarcy(icomp),Err)
    CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberDarcy,Region,GeometricField, &
      & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMFE_EQUATIONS_SET_DARCY_EQUATION_TYPE, &
      & CMFE_EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE],EquationsSetFieldDarcyUserNumber,EquationsSetFieldDarcy(icomp), &
      & EquationsSetDarcy(icomp),Err)

!Create materials_field for coupling coefficients acccording to create_equations_set_field
!  for first material field (auto-)create it, subsequently just pass it in as we do for the
!  shared dependent field between daryc and elaasticy


    CALL cmfe_EquationsSet_CreateFinish(EquationsSetDarcy(icomp),Err)

    CALL cmfe_Field_ParameterSetUpdateConstant(EquationsSetFieldDarcy(icomp),CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,1,icomp,Err)
    CALL cmfe_Field_ParameterSetUpdateConstant(EquationsSetFieldDarcy(icomp),CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,2,Ncompartments,Err)
  ENDDO
  write (*,*) "Equations set fields now created."
  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS
! 
!     !Initialise dependent field (velocity components)
!     DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!       CALL cmfe_Field_ComponentValuesInitialise(DependentFieldDarcy(icomp),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
!         & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
!     ENDDO
!   ENDDO

   CALL cmfe_Field_Initialise(DependentFieldDarcy,Err)
 

    CALL cmfe_Field_CreateStart(DependentFieldUserNumberDarcy,Region,DependentFieldDarcy,Err)

    CALL cmfe_Field_TypeSet(DependentFieldDarcy,CMFE_FIELD_GENERAL_TYPE,Err)  
    CALL cmfe_Field_MeshDecompositionSet(DependentFieldDarcy,Decomposition,Err)
    CALL cmfe_Field_GeometricFieldSet(DependentFieldDarcy,GeometricField,Err) 
    CALL cmfe_Field_DependentTypeSet(DependentFieldDarcy,CMFE_FIELD_DEPENDENT_TYPE,Err) 

    CALL cmfe_Field_NumberOfVariablesSet(DependentFieldDarcy,2*Ncompartments,Err) 
    !create two variables for each compartment
    ALLOCATE(VariableTypes(2*Ncompartments))
    DO num_var=1,Ncompartments
       VariableTypes(2*num_var-1)=CMFE_FIELD_U_VARIABLE_TYPE+(CMFE_FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
       VariableTypes(2*num_var)=CMFE_FIELD_DELUDELN_VARIABLE_TYPE+(CMFE_FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
    ENDDO
    CALL cmfe_Field_VariableTypesSet(DependentFieldDarcy,VariableTypes,Err) 

    DO icompartment=1,2*Ncompartments
      !set dimension type
      CALL cmfe_Field_DimensionSet(DependentFieldDarcy,VariableTypes(icompartment), &
         & CMFE_FIELD_VECTOR_DIMENSION_TYPE,Err)
      CALL cmfe_Field_NumberOfComponentsSet(DependentFieldDarcy,VariableTypes(icompartment),NUMBER_OF_DIMENSIONS+1,Err)
    ENDDO
    DO icompartment=1,2*Ncompartments
     DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDarcy,VariableTypes(icompartment),COMPONENT_NUMBER, & 
         & MESH_COMPONENT_NUMBER_VELOCITY,Err)
     ENDDO
    ENDDO
    DO icompartment=1,2*Ncompartments
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDarcy,VariableTypes(icompartment),NUMBER_OF_DIMENSIONS+1, & 
         & MESH_COMPONENT_NUMBER_PRESSURE,Err)
    ENDDO
    CALL cmfe_Field_CreateFinish(DependentFieldDarcy,Err)
    
    DO icomp=1,Ncompartments
     
      CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetDarcy(icomp),DependentFieldUserNumberDarcy,DependentFieldDarcy,Err)
      CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetDarcy(icomp),Err)

    ENDDO



  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  DO icomp = 1,Ncompartments
    !Create the equations set materials field variables for Static Darcy
    MaterialsFieldUserNumberDarcy = 400+icomp
    CALL cmfe_Field_Initialise(MaterialsFieldDarcy(icomp),Err)
    CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetDarcy(icomp),MaterialsFieldUserNumberDarcy, & 
      & MaterialsFieldDarcy(icomp),Err)
    !Finish the equations set materials field variables
    CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetDarcy(icomp),Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDarcy(icomp),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDarcy(icomp),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !EQUATIONS

  DO icomp = 1,Ncompartments
    !Create the equations set equations
    CALL cmfe_Equations_Initialise(EquationsDarcy(icomp),Err)
    CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetDarcy(icomp),EquationsDarcy(icomp),Err)
    !Set the equations matrices sparsity type
    CALL cmfe_Equations_SparsityTypeSet(EquationsDarcy(icomp),CMFE_EQUATIONS_SPARSE_MATRICES,Err)
    !Set the equations set output
    CALL cmfe_Equations_OutputTypeSet(EquationsDarcy(icomp),EQUATIONS_DARCY_OUTPUT,Err)
    !Finish the equations set equations
    CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetDarcy(icomp),Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS,CMFE_PROBLEM_DARCY_EQUATION_TYPE, &
    & CMFE_PROBLEM_STANDARD_DARCY_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(LinearSolverDarcy,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
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
  CALL cmfe_Solver_Initialise(LinearSolverDarcy,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsDarcy,Err)

  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the Darcy solver equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverDarcyUserNumber,LinearSolverDarcy,Err)
  CALL cmfe_Solver_SolverEquationsGet(LinearSolverDarcy,SolverEquationsDarcy,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsDarcy,CMFE_SOLVER_SPARSE_MATRICES,Err)

  DO icomp=1,Ncompartments
    CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsDarcy,EquationsSetDarcy(icomp), &
      & EquationsSetIndices(icomp),Err) !extend EquationsSetDarcy(icomp),EquationsSetIndices(icomp) to be arrays
  ENDDO

  !
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsDarcy,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)
  DO icomp = 1,Ncompartments

!     !--- BCs on normal velocity only
!     CONDITION = CMFE_BOUNDARY_CONDITION_MOVED_WALL
! 
!     IF( CM%D==2_CMISSIntg ) THEN
!       DO NODE_NUMBER=1_CMISSIntg,NUMBER_OF_NODES_GEOMETRY
!         COORD_X = CM%N(NODE_NUMBER,1_CMISSIntg)
!         COORD_Y = CM%N(NODE_NUMBER,2_CMISSIntg)
! 
!         IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) ) THEN
!           !x-velocity
!           VALUE = 1.0_CMISSRP
!           CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_NO_GLOBAL_DERIV, & 
!             & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
!           !x-velocity
!           VALUE = 1.0_CMISSRP
!           CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_NO_GLOBAL_DERIV, & 
!             & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
!           !y-velocity
!           VALUE = 2.0_CMISSRP
!           CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_NO_GLOBAL_DERIV, & 
!             & NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!           !y-velocity
!           VALUE = 2.0_CMISSRP
!           CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_NO_GLOBAL_DERIV, & 
!             & NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!       END DO
!     ELSE IF( CM%D==3_CMISSIntg ) THEN
!       DO NODE_NUMBER=1_CMISSIntg,NUMBER_OF_NODES_GEOMETRY
!         COORD_X = CM%N(NODE_NUMBER,1_CMISSIntg)
!         COORD_Y = CM%N(NODE_NUMBER,2_CMISSIntg)
!         COORD_Z = CM%N(NODE_NUMBER,3_CMISSIntg)
! 
!         IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) ) THEN
!           !x-velocity
!           VALUE = 2.0_CMISSRP
!           CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_NO_GLOBAL_DERIV, & 
!             & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
!           !x-velocity
!           VALUE = 1.0_CMISSRP
!           CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_NO_GLOBAL_DERIV, & 
!             & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
!           !y-velocity
!           VALUE = 1.0_CMISSRP
!           CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_NO_GLOBAL_DERIV, & 
!             & NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!           !y-velocity
!           VALUE = 2.0_CMISSRP
!           CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_NO_GLOBAL_DERIV, & 
!             & NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_Z-DOMAIN_Z1) < GEOMETRY_TOLERANCE) ) THEN
!           !z-velocity
!           VALUE = 1.0_CMISSRP
!           CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_NO_GLOBAL_DERIV, & 
!             & NODE_NUMBER,3_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) THEN
!           !z-velocity
!           VALUE = 1.0_CMISSRP
!           CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_NO_GLOBAL_DERIV, & 
!             & NODE_NUMBER,3_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!       END DO
!     END IF

  ENDDO
  !Finish the creation of the equations set boundary conditions for Darcy
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsDarcy,Err)

  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn off PETSc error handling
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
    CALL cmfe_Fields_NodesExport(Fields,"StaticDarcy","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"StaticDarcy","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF

  !Finialise CMISS
!   CALL cmfe_Finalise(Err)

  IF (ALLOCATED(EquationsSetFieldDarcy)) DEALLOCATE(EquationsSetFieldDarcy)
  IF (ALLOCATED(EquationsSetDarcy)) DEALLOCATE(EquationsSetDarcy)
  IF (ALLOCATED(MaterialsFieldDarcy)) DEALLOCATE(MaterialsFieldDarcy)
  IF (ALLOCATED(EquationsSetIndices)) DEALLOCATE(EquationsSetIndices)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM DARCYSTATICEXAMPLE
