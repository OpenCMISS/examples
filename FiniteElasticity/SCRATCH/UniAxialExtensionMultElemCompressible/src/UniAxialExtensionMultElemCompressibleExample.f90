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

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

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
!   BASIS_XI_INTERPOLATION_SOLID=CMISSBasisLinearLagrangeInterpolation
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
  LINEAR_SOLVER_SOLID_OUTPUT_TYPE=CMISSSolverNoOutput
  NONLINEAR_SOLVER_SOLID_OUTPUT_TYPE=CMISSSolverProgressOutput
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

  !CMISSAllDiagType/CMISSInDiagType/CMISSFromDiagType
!   CALL CMISSDiagnosticsSetOn(CMISSInDiagType,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISSAllTimingType/CMISSInTimingType/CMISSFromTimingType
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !REGION

  !Start the creation of a new region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system as defined above
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !BASES

  !Start the creation of new bases: Geometry
  MESH_NUMBER_OF_COMPONENTS=1
  CALL CMISSBasisTypeInitialise(BasisGeometry,Err)
  CALL CMISSBasisCreateStart(BASIS_NUMBER_GEOMETRY,BasisGeometry,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasisTypeSet(BasisGeometry,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL CMISSBasisNumberOfXiSet(BasisGeometry,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  CALL CMISSBasisInterpolationXiSet(BasisGeometry,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY, & 
    & BASIS_XI_INTERPOLATION_GEOMETRY/),Err)                         
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisGeometry,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY, &
    & BASIS_XI_GAUSS_GEOMETRY/),Err)
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(BasisGeometry,Err)

  !Start the creation of the solid basis
  CALL CMISSBasisTypeInitialise(BasisSolid,Err)
  CALL CMISSBasisCreateStart(BASIS_NUMBER_SOLID,BasisSolid,Err) 

  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasisTypeSet(BasisSolid,BASIS_TYPE,Err)
!   CALL CMISSBasisTypeSet(BasisSolid,CMISSBasisLagrangeHermiteTPType,Err)

  CALL CMISSBasisNumberOfXiSet(BasisSolid,NUMBER_OF_DIMENSIONS,Err)

  CALL CMISSBasisInterpolationXiSet(BasisSolid,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY, & 
    & BASIS_XI_INTERPOLATION_GEOMETRY/),Err)                         
!   CALL CMISSBasisInterpolationXiSet(BasisSolid,(/CMISSBasisLinearLagrangeInterpolation, &
!     & CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation/),Err)

  CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisSolid,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY, &
    & BASIS_XI_GAUSS_GEOMETRY/),Err)

  CALL CMISSBasisCreateFinish(BasisSolid,Err)

  !
  !================================================================================================================================
  !

  !MESH

  TotalNumberOfSolidNodes = NUMBER_OF_NODES_GEOMETRY
  TOTAL_NUMBER_OF_ALL_NODES = TOTAL_NUMBER_OF_NODES + TotalNumberOfSolidNodes

  !Start the creation of mesh nodes
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSNodesCreateStart(Region,TOTAL_NUMBER_OF_ALL_NODES,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMeshNumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMeshNumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)
  !
  CALL CMISSMeshElementsTypeInitialise(MeshElementsGeometry,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsSolid,Err)
  MESH_COMPONENT_NUMBER_GEOMETRY=1
  SolidMeshComponenetNumber = 1
  !Specify spatial mesh component
  CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_GEOMETRY,BasisGeometry,MeshElementsGeometry,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL CMISSMeshElementsNodesSet(MeshElementsGeometry,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
  ENDDO
  CALL CMISSMeshElementsCreateFinish(MeshElementsGeometry,Err)

  !Specify Solid Mesh Component
  MeshElementsSolid=MeshElementsGeometry
  SolidMeshComponenetNumber=MESH_COMPONENT_NUMBER_GEOMETRY

  !Finish the creation of the mesh
  CALL CMISSMeshCreateFinish(Mesh,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition:
  !All mesh components share the same decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,DomainUserNumber,Err)
  ! ??? Above, this should be: 'NumberOfDomains' rather than 'DomainUserNumber' ???
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
      & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
  ENDDO

  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
        & CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSFieldParameterSetUpdateStart(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create a decomposition

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSFieldTypeInitialise(GeometricFieldSolid,Err)
  CALL CMISSFieldCreateStart(FieldGeometrySolidUserNumber,Region,GeometricFieldSolid,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricFieldSolid,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricFieldSolid,CMISSFieldGeometricType,Err)  
  CALL CMISSFieldNumberOfVariablesSet(GeometricFieldSolid,FieldGeometrySolidNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricFieldSolid,CMISSFieldUVariableType,FieldGeometrySolidNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,1,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,2,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,3,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldCreateFinish(GeometricFieldSolid,Err)

!---
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSFieldParameterSetUpdateNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
        & CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSFieldParameterSetUpdateStart(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
!---

  !Create a fibre field and attach it to the geometric field  
  CALL CMISSFieldTypeInitialise(FibreFieldSolid,Err)
  CALL CMISSFieldCreateStart(FieldFibreSolidUserNumber,Region,FibreFieldSolid,Err)
  CALL CMISSFieldTypeSet(FibreFieldSolid,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreFieldSolid,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(FibreFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreFieldSolid,FieldFibreSolidNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreFieldSolid,CMISSFieldUVariableType,FieldFibreSolidNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,1,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,2,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,3,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldCreateFinish(FibreFieldSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations_set
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSetSolid,Err)
  CALL CMISSEquationsSetCreateStart(EquationSetSolidUserNumber,Region,FibreFieldSolid,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetCompressibleFiniteElasticitySubtype, &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSetSolid,Err)
!   CALL CMISSEquationsSetSpecificationSet(EquationsSetSolid,CMISSEquationsSetElasticityClass, &
! !     & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetNoSubtype,Err)
! !     & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetMooneyRivlinSubtype,Err)
!     & ,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create a material field and attach it to the geometric field  
  CALL CMISSFieldTypeInitialise(MaterialFieldSolid,Err)
  !
  CALL CMISSFieldCreateStart(FieldMaterialSolidUserNumber,Region,MaterialFieldSolid,Err)
  !
  CALL CMISSFieldTypeSet(MaterialFieldSolid,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialFieldSolid,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(MaterialFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialFieldSolid,FieldMaterialSolidNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialFieldSolid,CMISSFieldUVariableType,FieldMaterialSolidNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldSolid,CMISSFieldUVariableType,1,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldSolid,CMISSFieldUVariableType,2,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldSolid,CMISSFieldUVariableType,3,SolidMeshComponenetNumber,Err)
  !
  CALL CMISSFieldCreateFinish(MaterialFieldSolid,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,6.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,1.0e7_CMISSDP,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,1.0e3_CMISSDP,Err)

  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetSolid,FieldMaterialSolidUserNumber,MaterialFieldSolid,Err)  
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetSolid,Err)


  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create a dependent field with two variables and three components
  CALL CMISSFieldTypeInitialise(DependentFieldSolid,Err)
  !
  CALL CMISSFieldCreateStart(FieldDependentSolidUserNumber,Region,DependentFieldSolid,Err)
  !
  CALL CMISSFieldTypeSet(DependentFieldSolid,CMISSFieldGeneralType,Err)  
  CALL CMISSFieldMeshDecompositionSet(DependentFieldSolid,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(DependentFieldSolid,GeometricFieldSolid,Err) 
  CALL CMISSFieldDependentTypeSet(DependentFieldSolid,CMISSFieldDependentType,Err) 
  CALL CMISSFieldNumberOfVariablesSet(DependentFieldSolid,FieldDependentSolidNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldUVariableType,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,FieldDependentSolidNumberOfComponents,Err)
  !
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,1,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,2,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,3,SolidMeshComponenetNumber,Err)  
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,4,CMISSFieldElementBasedInterpolation,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,1,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,2,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,3,SolidMeshComponenetNumber,Err)  
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4, &
!     & CMISSFieldElementBasedInterpolation,Err)
  !
  CALL CMISSFieldCreateFinish(DependentFieldSolid,Err)  
  !
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetSolid,FieldDependentSolidUserNumber,DependentFieldSolid,Err) 
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetSolid,Err)


  !
  !================================================================================================================================
  !

  !EQUATIONS

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsSolid,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetSolid,EquationsSolid,Err)
  CALL CMISSEquationsSparsityTypeSet(EquationsSolid,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(EquationsSolid,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetSolid,Err)   


  !
  !================================================================================================================================
  !

  !Initialise dependent field from undeformed geometry and displacement bcs 
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)
!   CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,-8.0_CMISSDP,Err)


  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a finite elasticity problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemElasticityClass,CMISSProblemFiniteElasticityType, &
    & CMISSProblemNoSubtype,Err)

  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)


  !
  !================================================================================================================================
  !

  !SOLVERS

!---Customary: begin

!   !Start the creation of the problem solvers
!   CALL CMISSSolverTypeInitialise(NonlinearSolverSolid,Err)
!   CALL CMISSSolverTypeInitialise(LinearSolverSolid,Err)
!   CALL CMISSProblemSolversCreateStart(Problem,Err)
! 
!   !Get the finite elasticity solver
!   CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,NonlinearSolverSolidUserNumber,NonlinearSolverSolid,Err)
!   CALL CMISSSolverOutputTypeSet(NonlinearSolverSolid,CMISSSolverProgressOutput,Err)
!   CALL CMISSSolverNewtonJacobianCalculationTypeSet(NonlinearSolverSolid,CMISSSolverNewtonJacobianFDCalculated,Err)
! 
!   CALL CMISSSolverNonLinearTypeSet(NonlinearSolverSolid,CMISSSolverNonlinearNewton,Err)
!   CALL CMISSSolverLibraryTypeSet(NonlinearSolverSolid,CMISSSolverPETScLibrary,Err)
! 
!   CALL CMISSSolverNewtonLinearSolverGet(NonlinearSolverSolid,LinearSolverSolid,Err)
!   CALL CMISSSolverLinearTypeSet(LinearSolverSolid,CMISSSolverLinearDirectSolveType,Err)
! 
!   !Finish the creation of the problem solver
!   CALL CMISSProblemSolversCreateFinish(Problem,Err)

!---Customary: end


!---tob

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(NonlinearSolverSolid,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverSolid,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  !Get the nonlinear static solver
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,NonlinearSolverSolidUserNumber,NonlinearSolverSolid,Err)
  !Set the nonlinear Jacobian type
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(NonlinearSolverSolid,CMISSSolverNewtonJacobianAnalyticCalculated,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(NonlinearSolverSolid,CMISSSolverNewtonJacobianFDCalculated,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(NonlinearSolverSolid,NONLINEAR_SOLVER_SOLID_OUTPUT_TYPE,Err)
  !Set the solver settings
  CALL CMISSSolverNewtonAbsoluteToleranceSet(NonlinearSolverSolid,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolverNewtonRelativeToleranceSet(NonlinearSolverSolid,RELATIVE_TOLERANCE,Err)
  !Get the nonlinear linear solver
  CALL CMISSSolverNewtonLinearSolverGet(NonlinearSolverSolid,LinearSolverSolid,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(LinearSolverSolid,LINEAR_SOLVER_SOLID_OUTPUT_TYPE,Err)


  !Set the solver settings
  IF(LINEAR_SOLVER_SOLID_DIRECT_FLAG) THEN
    CALL CMISSSolverLinearTypeSet(LinearSolverSolid,CMISSSolverLinearDirectSolveType,Err)
    CALL CMISSSolverLibraryTypeSet(LinearSolverSolid,CMISSSolverMUMPSLibrary,Err)
  ELSE
    CALL CMISSSolverLinearTypeSet(LinearSolverSolid,CMISSSolverLinearIterativeSolveType,Err)
    CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolverSolid,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolverLinearIterativeDivergenceToleranceSet(LinearSolverSolid,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeRelativeToleranceSet(LinearSolverSolid,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(LinearSolverSolid,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeGMRESRestartSet(LinearSolverSolid,RESTART_VALUE,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)
!---toe


  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL CMISSSolverTypeInitialise(NonlinearSolverSolid,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsSolid,Err)

  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !
  !Get the finite elasticity solver equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,NonlinearSolverSolidUserNumber,NonlinearSolverSolid,Err)
  CALL CMISSSolverSolverEquationsGet(NonlinearSolverSolid,SolverEquationsSolid,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsSolid,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsSolid,EquationsSetSolid,EquationsSetIndex,Err)
  !
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)


  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsSolid,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditionsSolid,Err)

  !--- BCs on displacement only
  CONDITION = CMISSBoundaryConditionDirichlet

  DO NODE_NUMBER=1_CMISSIntg,NUMBER_OF_NODES_GEOMETRY  !What if different number of nodes geometry and velocity ?
    COORD_X = CM%N(NODE_NUMBER,1_CMISSIntg)
    COORD_Y = CM%N(NODE_NUMBER,2_CMISSIntg)
    COORD_Z = CM%N(NODE_NUMBER,3_CMISSIntg)

    IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) ) THEN
      !x-displacement
      VALUE = 0.0_CMISSDP
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,1, &
        & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
    END IF
    !
    IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
      !x-displacement
      VALUE = 1.05_CMISSDP  ! * WIDTH
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,1, &
        & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
    END IF
    !
    IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
      !y-displacement
      VALUE = 0.0_CMISSDP
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,1, &
        & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
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
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,1, &
        & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
    END IF
    !
    IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) THEN
      !z-displacement
      ! do nothing ???
    END IF
  END DO

  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsSolid,Err)


  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn of PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL CMISSProblemSolve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"


  !
  !================================================================================================================================
  !

  !OUTPUT

  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"CompressibleFiniteElasticity","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"CompressibleFiniteElasticity","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF


  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM FINITEELASTICITYDARCYEXAMPLE
