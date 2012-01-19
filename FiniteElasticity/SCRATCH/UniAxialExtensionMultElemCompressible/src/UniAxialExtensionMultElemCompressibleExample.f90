!> \file
!> \authors Christian Michler
!> \brief This is an example program to solve a Compressible Finite Elastiticity equation using openCMISS calls and FLUID_MECHANICS_IO_ROUTINES.
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

!> \example FiniteElasticity/SCRATCH/UniAxialExtensionMultElemCompressible/src/UniAxialExtensionMultElemCompressibleExample.f90
!! Example program to solve FiniteElasticity equations using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/SCRATCH/UniAxialExtensionMultElemCompressible/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/SCRATCH/UniAxialExtensionMultElemCompressible/build-intel'>Linux GNU Build</a>
!!
!<

! ! 
! !  This example considers a weakly compressible Finite Elasticity problem
! ! 

!> Main program

PROGRAM FINITEELASTICITYDARCYEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OPENCMISS
  USE FLUID_MECHANICS_IO_ROUTINES
  USE MPI

#ifdef WIN32
  USE IFQWINCMISS
#endif

  !
  !================================================================================================================================
  !

  !PROGRAM VARIABLES AND TYPES

  IMPLICIT NONE

  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=13

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NonlinearSolverSolidUserNumber=1

  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_GEOMETRY
!   INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SOLID
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_GEOMETRY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_GEOMETRY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_GEOMETRY
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ALL_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE

  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: LINEAR_SOLVER_SOLID_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: NONLINEAR_SOLVER_SOLID_OUTPUT_TYPE

  REAL(CMISSDP) :: COORD_X, COORD_Y, COORD_Z
  REAL(CMISSDP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2
  REAL(CMISSDP) :: GEOMETRY_TOLERANCE

  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_SOLID_DIRECT_FLAG

  !CMISS variables

  !Regions
  TYPE(CMISSRegionType) :: Region
  TYPE(CMISSRegionType) :: WorldRegion
  !Coordinate systems
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  !Basis
  TYPE(CMISSBasisType) :: BasisGeometry
  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Elements
  TYPE(CMISSMeshElementsType) :: MeshElementsGeometry
  !Meshes
  TYPE(CMISSMeshType) :: Mesh
  !Decompositions
  TYPE(CMISSDecompositionType) :: Decomposition
  !Fields
  TYPE(CMISSFieldsType) :: Fields
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: NonlinearSolverSolid
  TYPE(CMISSSolverType) :: LinearSolverSolid

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: Err


  INTEGER(CMISSIntg) :: DIAG_LEVEL_LIST(5)
!   CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(8) !,TIMING_ROUTINE_LIST(1)
  CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(1) !,TIMING_ROUTINE_LIST(1)

  
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !

  !Program variables and types (finite elasticity part)

  !Test program parameters

  INTEGER(CMISSIntg) :: BASIS_NUMBER_SOLID
! 
  INTEGER(CMISSIntg) :: TotalNumberOfSolidNodes
  INTEGER(CMISSIntg) :: SolidMeshComponenetNumber
! 
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometrySolidUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometrySolidNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometrySolidNumberOfComponents=3
! 
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreSolidUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreSolidNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreSolidNumberOfComponents=3
! 
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialSolidUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialSolidNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialSolidNumberOfComponents=3
! 
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidNumberOfVariables=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidNumberOfComponents=3
! 
  INTEGER(CMISSIntg), PARAMETER :: EquationSetSolidUserNumber=1
! 
!   !Program types
! 
! 
!   !Program variables
! 
!   !CMISS variables
! 
  TYPE(CMISSBasisType) :: BasisSolid
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsSolid
  TYPE(CMISSEquationsType) :: EquationsSolid
  TYPE(CMISSEquationsSetType) :: EquationsSetSolid
  TYPE(CMISSFieldType) :: GeometricFieldSolid,FibreFieldSolid,MaterialFieldSolid,DependentFieldSolid,EquationsSetField
  TYPE(CMISSSolverEquationsType) :: SolverEquationsSolid
  TYPE(CMISSMeshElementsType) :: MeshElementsSolid

  !End - Program variables and types (finite elasticity part)

  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !


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

  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  !
  !================================================================================================================================
  !

  !PROBLEM CONTROL PANEL

  !Import cmHeart mesh information
  CALL FLUID_MECHANICS_IO_READ_CMHEART(CM,Err)  
  BASIS_NUMBER_GEOMETRY=CM%ID_M
  BASIS_NUMBER_SOLID=CM%ID_P
  NUMBER_OF_DIMENSIONS=CM%D
  BASIS_TYPE=CM%IT_T
  BASIS_XI_INTERPOLATION_GEOMETRY=CM%IT_M
!   BASIS_XI_INTERPOLATION_SOLID=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  NUMBER_OF_NODES_GEOMETRY=CM%N_M
  TOTAL_NUMBER_OF_NODES=CM%N_T
  TOTAL_NUMBER_OF_ELEMENTS=CM%E_T
  NUMBER_OF_ELEMENT_NODES_GEOMETRY=CM%EN_M
  !Set domain dimensions
  DOMAIN_X1 =  0.0_CMISSDP
  DOMAIN_X2 =  1.0_CMISSDP
  DOMAIN_Y1 =  0.0_CMISSDP
  DOMAIN_Y2 =  1.0_CMISSDP
  DOMAIN_Z1 =  0.0_CMISSDP
  DOMAIN_Z2 =  1.0_CMISSDP
  !Set geometric tolerance
  GEOMETRY_TOLERANCE = 1.0E-12_CMISSDP
  !Set material parameters
  !Set number of Gauss points (Mind that also material field may be interpolated)
  BASIS_XI_GAUSS_GEOMETRY=3 !4
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_SOLID_OUTPUT_TYPE=CMISS_SOLVER_NO_OUTPUT
  NONLINEAR_SOLVER_SOLID_OUTPUT_TYPE=CMISS_SOLVER_PROGRESS_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  !Set solver parameters
  LINEAR_SOLVER_SOLID_DIRECT_FLAG=.TRUE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E5_CMISSDP !default: 1.0E5
  MAXIMUM_ITERATIONS=10000_CMISSIntg !default: 100000
  RESTART_VALUE=30_CMISSIntg !default: 30
  LINESEARCH_ALPHA=1.0_CMISSDP


  !
  !================================================================================================================================
  !

  !Set diagnostics

  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5

  DIAG_ROUTINE_LIST(1)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"

  !CMISS_ALL_DIAG_TYPE/CMISS_IN_DIAG_TYPE/CMISS_FROM_DIAG_TYPE
!   CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISS_ALL_TIMING_TYPE/CMISS_IN_TIMING_TYPE/CMISS_FROM_TIMING_TYPE
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !REGION

  !Start the creation of a new region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system as defined above
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !BASES

  !Start the creation of new bases: Geometry
  MESH_NUMBER_OF_COMPONENTS=1
  CALL CMISSBasis_Initialise(BasisGeometry,Err)
  CALL CMISSBasis_CreateStart(BASIS_NUMBER_GEOMETRY,BasisGeometry,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasis_TypeSet(BasisGeometry,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL CMISSBasis_NumberOfXiSet(BasisGeometry,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  CALL CMISSBasis_InterpolationXiSet(BasisGeometry,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY, & 
    & BASIS_XI_INTERPOLATION_GEOMETRY/),Err)                         
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisGeometry,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY, &
    & BASIS_XI_GAUSS_GEOMETRY/),Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(BasisGeometry,Err)

  !Start the creation of the solid basis
  CALL CMISSBasis_Initialise(BasisSolid,Err)
  CALL CMISSBasis_CreateStart(BASIS_NUMBER_SOLID,BasisSolid,Err) 

  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasis_TypeSet(BasisSolid,BASIS_TYPE,Err)
!   CALL CMISSBasis_TypeSet(BasisSolid,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)

  CALL CMISSBasis_NumberOfXiSet(BasisSolid,NUMBER_OF_DIMENSIONS,Err)

  CALL CMISSBasis_InterpolationXiSet(BasisSolid,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY, & 
    & BASIS_XI_INTERPOLATION_GEOMETRY/),Err)                         
!   CALL CMISSBasis_InterpolationXiSet(BasisSolid,(/CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
!     & CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION/),Err)

  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisSolid,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY, &
    & BASIS_XI_GAUSS_GEOMETRY/),Err)

  CALL CMISSBasis_CreateFinish(BasisSolid,Err)

  !
  !================================================================================================================================
  !

  !MESH

  TotalNumberOfSolidNodes = NUMBER_OF_NODES_GEOMETRY
  TOTAL_NUMBER_OF_ALL_NODES = TOTAL_NUMBER_OF_NODES + TotalNumberOfSolidNodes

  !Start the creation of mesh nodes
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(Region,TOTAL_NUMBER_OF_ALL_NODES,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMesh_NumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)
  !
  CALL CMISSMeshElements_Initialise(MeshElementsGeometry,Err)
  CALL CMISSMeshElements_Initialise(MeshElementsSolid,Err)
  MESH_COMPONENT_NUMBER_GEOMETRY=1
  SolidMeshComponenetNumber = 1
  !Specify spatial mesh component
  CALL CMISSMeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_GEOMETRY,BasisGeometry,MeshElementsGeometry,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL CMISSMeshElements_NodesSet(MeshElementsGeometry,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(MeshElementsGeometry,Err)

  !Specify Solid Mesh Component
  MeshElementsSolid=MeshElementsGeometry
  SolidMeshComponenetNumber=MESH_COMPONENT_NUMBER_GEOMETRY

  !Finish the creation of the mesh
  CALL CMISSMesh_CreateFinish(Mesh,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition:
  !All mesh components share the same decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,DomainUserNumber,Err)
  ! ??? Above, this should be: 'NumberOfDomains' rather than 'DomainUserNumber' ???
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field type
  CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL CMISSField_ScalingTypeSet(GeometricField,CMISS_FIELD_NO_SCALING,Err)
  !Set the mesh component to be used by the field components.

  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
  ENDDO

  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
        & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSField_ParameterSetUpdateStart(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create a decomposition

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSField_Initialise(GeometricFieldSolid,Err)
  CALL CMISSField_CreateStart(FieldGeometrySolidUserNumber,Region,GeometricFieldSolid,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricFieldSolid,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricFieldSolid,CMISS_FIELD_GEOMETRIC_TYPE,Err)  
  CALL CMISSField_NumberOfVariablesSet(GeometricFieldSolid,FieldGeometrySolidNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometrySolidNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,SolidMeshComponenetNumber,Err)
  CALL CMISSField_CreateFinish(GeometricFieldSolid,Err)

!---
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSField_ParameterSetUpdateNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
        & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSField_ParameterSetUpdateStart(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
!---

  !Create a fibre field and attach it to the geometric field  
  CALL CMISSField_Initialise(FibreFieldSolid,Err)
  CALL CMISSField_CreateStart(FieldFibreSolidUserNumber,Region,FibreFieldSolid,Err)
  CALL CMISSField_TypeSet(FibreFieldSolid,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreFieldSolid,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(FibreFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreFieldSolid,FieldFibreSolidNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,FieldFibreSolidNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,SolidMeshComponenetNumber,Err)
  CALL CMISSField_CreateFinish(FibreFieldSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations_set
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetSolid,Err)
  CALL CMISSEquationsSet_CreateStart(EquationSetSolidUserNumber,Region,FibreFieldSolid,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSetSolid,Err)
!   CALL CMISSEquationsSet_SpecificationSet(EquationsSetSolid,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
! !     & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_NO_SUBTYPE,Err)
! !     & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,Err)
!     & ,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create a material field and attach it to the geometric field  
  CALL CMISSField_Initialise(MaterialFieldSolid,Err)
  !
  CALL CMISSField_CreateStart(FieldMaterialSolidUserNumber,Region,MaterialFieldSolid,Err)
  !
  CALL CMISSField_TypeSet(MaterialFieldSolid,CMISS_FIELD_MATERIAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(MaterialFieldSolid,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(MaterialFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialFieldSolid,FieldMaterialSolidNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,FieldMaterialSolidNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,SolidMeshComponenetNumber,Err)
  !
  CALL CMISSField_CreateFinish(MaterialFieldSolid,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & 2.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
    & 6.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3, &
    & 1.0e7_CMISSDP,Err)
!   CALL CMISSField_ComponentValuesInitialise(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,1.0e3_CMISSDP,Err)

  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetSolid,FieldMaterialSolidUserNumber,MaterialFieldSolid,Err)  
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetSolid,Err)


  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create a dependent field with two variables and three components
  CALL CMISSField_Initialise(DependentFieldSolid,Err)
  !
  CALL CMISSField_CreateStart(FieldDependentSolidUserNumber,Region,DependentFieldSolid,Err)
  !
  CALL CMISSField_TypeSet(DependentFieldSolid,CMISS_FIELD_GENERAL_TYPE,Err)  
  CALL CMISSField_MeshDecompositionSet(DependentFieldSolid,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentFieldSolid,GeometricFieldSolid,Err) 
  CALL CMISSField_DependentTypeSet(DependentFieldSolid,CMISS_FIELD_DEPENDENT_TYPE,Err) 
  CALL CMISSField_NumberOfVariablesSet(DependentFieldSolid,FieldDependentSolidNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE, &
    & FieldDependentSolidNumberOfComponents,Err)
  !
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,SolidMeshComponenetNumber,Err)  
!   CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3,SolidMeshComponenetNumber,Err)  
!   CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
!     & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  !
  CALL CMISSField_CreateFinish(DependentFieldSolid,Err)  
  !
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetSolid,FieldDependentSolidUserNumber,DependentFieldSolid,Err) 
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetSolid,Err)


  !
  !================================================================================================================================
  !

  !EQUATIONS

  !Create the equations set equations
  CALL CMISSEquations_Initialise(EquationsSolid,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetSolid,EquationsSolid,Err)
  CALL CMISSEquations_SparsityTypeSet(EquationsSolid,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(EquationsSolid,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetSolid,Err)   


  !
  !================================================================================================================================
  !

  !Initialise dependent field from undeformed geometry and displacement bcs 
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
!   CALL CMISSField_ComponentValuesInitialise(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,-8.0_CMISSDP,Err)


  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a finite elasticity problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_ELASTICITY_CLASS,CMISS_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMISS_PROBLEM_NO_SUBTYPE,Err)

  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)


  !
  !================================================================================================================================
  !

  !SOLVERS

!---Customary: begin

!   !Start the creation of the problem solvers
!   CALL CMISSSolver_Initialise(NonlinearSolverSolid,Err)
!   CALL CMISSSolver_Initialise(LinearSolverSolid,Err)
!   CALL CMISSProblem_SolversCreateStart(Problem,Err)
! 
!   !Get the finite elasticity solver
!   CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,NonlinearSolverSolidUserNumber,NonlinearSolverSolid,Err)
!   CALL CMISSSolver_OutputTypeSet(NonlinearSolverSolid,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
!   CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolverSolid,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
! 
!   CALL CMISSSolverNonLinearTypeSet(NonlinearSolverSolid,CMISS_SOLVER_NONLINEAR_NEWTON,Err)
!   CALL CMISSSolver_LibraryTypeSet(NonlinearSolverSolid,CMISS_SOLVER_PETSC_LIBRARY,Err)
! 
!   CALL CMISSSolver_NewtonLinearSolverGet(NonlinearSolverSolid,LinearSolverSolid,Err)
!   CALL CMISSSolver_LinearTypeSet(LinearSolverSolid,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
! 
!   !Finish the creation of the problem solver
!   CALL CMISSProblem_SolversCreateFinish(Problem,Err)

!---Customary: end


!---tob

  !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(NonlinearSolverSolid,Err)
  CALL CMISSSolver_Initialise(LinearSolverSolid,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  !Get the nonlinear static solver
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,NonlinearSolverSolidUserNumber,NonlinearSolverSolid,Err)
  !Set the nonlinear Jacobian type
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolverSolid,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolverSolid,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(NonlinearSolverSolid,NONLINEAR_SOLVER_SOLID_OUTPUT_TYPE,Err)
  !Set the solver settings
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(NonlinearSolverSolid,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(NonlinearSolverSolid,RELATIVE_TOLERANCE,Err)
  !Get the nonlinear linear solver
  CALL CMISSSolver_NewtonLinearSolverGet(NonlinearSolverSolid,LinearSolverSolid,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(LinearSolverSolid,LINEAR_SOLVER_SOLID_OUTPUT_TYPE,Err)


  !Set the solver settings
  IF(LINEAR_SOLVER_SOLID_DIRECT_FLAG) THEN
    CALL CMISSSolver_LinearTypeSet(LinearSolverSolid,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL CMISSSolver_LibraryTypeSet(LinearSolverSolid,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL CMISSSolver_LinearTypeSet(LinearSolverSolid,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolverSolid,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolver_LinearIterativeDivergenceToleranceSet(LinearSolverSolid,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(LinearSolverSolid,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(LinearSolverSolid,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeGMRESRestartSet(LinearSolverSolid,RESTART_VALUE,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)
!---toe


  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(NonlinearSolverSolid,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsSolid,Err)

  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !
  !Get the finite elasticity solver equations
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,NonlinearSolverSolidUserNumber,NonlinearSolverSolid,Err)
  CALL CMISSSolver_SolverEquationsGet(NonlinearSolverSolid,SolverEquationsSolid,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsSolid,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsSolid,EquationsSetSolid,EquationsSetIndex,Err)
  !
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)


  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsSolid,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditionsSolid,Err)

  !--- BCs on displacement only
  CONDITION = CMISS_BOUNDARY_CONDITION_DIRICHLET

  DO NODE_NUMBER=1_CMISSIntg,NUMBER_OF_NODES_GEOMETRY  !What if different number of nodes geometry and velocity ?
    COORD_X = CM%N(NODE_NUMBER,1_CMISSIntg)
    COORD_Y = CM%N(NODE_NUMBER,2_CMISSIntg)
    COORD_Z = CM%N(NODE_NUMBER,3_CMISSIntg)

    IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) ) THEN
      !x-displacement
      VALUE = 0.0_CMISSDP
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
        & NODE_NUMBER,1_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    END IF
    !
    IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
      !x-displacement
      VALUE = 1.05_CMISSDP  ! * WIDTH
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
        & NODE_NUMBER,1_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    END IF
    !
    IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
      !y-displacement
      VALUE = 0.0_CMISSDP
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
        & NODE_NUMBER,2_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    END IF
    !
    IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
      !y-displacement
      ! do nothing ???
    END IF
    !
    IF( (ABS(COORD_Z-DOMAIN_Z1) < GEOMETRY_TOLERANCE) ) THEN
      !z-displacement
      VALUE = 0.0_CMISSDP
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
        & NODE_NUMBER,3_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    END IF
    !
    IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) THEN
      !z-displacement
      ! do nothing ???
    END IF
  END DO

  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsSolid,Err)


  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn of PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL CMISSProblem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"


  !
  !================================================================================================================================
  !

  !OUTPUT

  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"CompressibleFiniteElasticity","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"CompressibleFiniteElasticity","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF


  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM FINITEELASTICITYDARCYEXAMPLE
