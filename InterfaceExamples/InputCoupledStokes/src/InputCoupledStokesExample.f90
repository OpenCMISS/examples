!> \file
!> \author Chris Bradley
!> \brief This is an example program which solves a weakly coupled Stokes equation in two regions using OpenCMISS calls.
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

!> \example InterfaceExamples/CoupledStokes/src/CoupledStoesExample.f90
!! Example program which sets up a field in two regions using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/InterfaceExamples/CoupledStokes/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/InterfaceExamples/CoupledStokes/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM COUPLEDSTOKES

  USE OPENCMISS
  USE FLUID_MECHANICS_IO_ROUTINES
  
#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystem1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystem2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: Region1UserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: Region2UserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: Basis1UserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: Basis2UserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: InterfaceBasisUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: InterfaceGeneratedMeshUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: Mesh1UserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: Mesh2UserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: InterfaceMeshUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: Decomposition1UserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: Decomposition2UserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: InterfaceDecompositionUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: GeometricField1UserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: GeometricField2UserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: InterfaceGeometricFieldUserNumber=19
  INTEGER(CMISSIntg), PARAMETER :: EquationsSet1UserNumber=20
  INTEGER(CMISSIntg), PARAMETER :: EquationsSet2UserNumber=21
  INTEGER(CMISSIntg), PARAMETER :: DependentField1UserNumber=22
  INTEGER(CMISSIntg), PARAMETER :: DependentField2UserNumber=23
  INTEGER(CMISSIntg), PARAMETER :: InterfaceUserNumber=24
  INTEGER(CMISSIntg), PARAMETER :: InterfaceConditionUserNumber=25
  INTEGER(CMISSIntg), PARAMETER :: LagrangeFieldUserNumber=26
  INTEGER(CMISSIntg), PARAMETER :: CoupledProblemUserNumber=27
  INTEGER(CMISSIntg), PARAMETER :: InterfaceMappingBasisUserNumber=28
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetField1UserNumber=40
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetField2UserNumber=41

  INTEGER(CMISSIntg), PARAMETER :: MaterialsComponentUserNumber1=1 !MU
  INTEGER(CMISSIntg), PARAMETER :: MaterialsComponentUserNumber2=2 !RHO
 
  !Program types

  TYPE(EXPORT_CONTAINER):: CM1,CM2,CM3
  TYPE(COUPLING_PARAMETERS):: CMX
  TYPE(BOUNDARY_PARAMETERS):: BC
  
  !Program variables
  INTEGER(CMISSIntg) :: EquationsSet1Index,EquationsSet2Index
  INTEGER(CMISSIntg) :: BoundaryNodeDomain
  INTEGER(CMISSIntg) :: InterfaceConditionIndex
  INTEGER(CMISSIntg) :: Mesh1Index,Mesh2Index
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,ic_idx
  REAL(CMISSDP) :: VALUE

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS1,NUMBER_OF_DIMENSIONS2,NUMBER_OF_DIMENSIONS_INTERFACE
  
  INTEGER(CMISSIntg) :: BASIS_TYPE1,BASIS_TYPE2,BASIS_TYPE_INTERFACE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_SPACE1,BASIS_NUMBER_SPACE2,BASIS_NUMBER_SPACE_INTERFACE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_VELOCITY1,BASIS_NUMBER_VELOCITY2
  INTEGER(CMISSIntg) :: BASIS_NUMBER_PRESSURE1,BASIS_NUMBER_PRESSURE2
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_SPACE1,BASIS_XI_GAUSS_SPACE2,BASIS_XI_GAUSS_INTERFACE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_VELOCITY1,BASIS_XI_GAUSS_VELOCITY2
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_PRESSURE1,BASIS_XI_GAUSS_PRESSURE2
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SPACE1,BASIS_XI_INTERPOLATION_SPACE2,BASIS_XI_INTERPOLATION_INTERFACE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_VELOCITY1,BASIS_XI_INTERPOLATION_VELOCITY2
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_PRESSURE1,BASIS_XI_INTERPOLATION_PRESSURE2
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS1,MESH_NUMBER_OF_COMPONENTS2,MESH_NUMBER_OF_COMPONENTS_INTERFACE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_SPACE1,MESH_COMPONENT_NUMBER_SPACE2,MESH_COMPONENT_NUMBER_INTERFACE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_VELOCITY1,MESH_COMPONENT_NUMBER_VELOCITY2
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_PRESSURE1,MESH_COMPONENT_NUMBER_PRESSURE2
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_SPACE1,NUMBER_OF_NODES_SPACE2,NUMBER_OF_NODES_INTERFACE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_VELOCITY1,NUMBER_OF_NODES_VELOCITY2
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_PRESSURE1,NUMBER_OF_NODES_PRESSURE2
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_SPACE1,NUMBER_OF_ELEMENT_NODES_SPACE2,NUMBER_OF_ELEMENT_NODES_INTERFACE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_VELOCITY1,NUMBER_OF_ELEMENT_NODES_VELOCITY2
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_PRESSURE1,NUMBER_OF_ELEMENT_NODES_PRESSURE2
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES1,TOTAL_NUMBER_OF_NODES2,TOTAL_NUMBER_OF_NODES_INTERFACE
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS1,TOTAL_NUMBER_OF_ELEMENTS2,TOTAL_NUMBER_OF_ELEMENTS_INTERFACE
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER,COMPONENT_NUMBER,NODE_NUMBER,NODE_COUNTER,CONDITION
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES1,NUMBER_OF_FIXED_WALL_NODES2
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES1,NUMBER_OF_INLET_WALL_NODES2

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES1
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES1
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS1(3)
  LOGICAL :: FIXED_WALL_NODES1_FLAG=.TRUE.
  LOGICAL :: INLET_WALL_NODES1_FLAG=.TRUE.

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES2
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES2
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS2(3)
  LOGICAL :: FIXED_WALL_NODES2_FLAG=.TRUE.
  LOGICAL :: INLET_WALL_NODES2_FLAG=.TRUE.


  REAL(CMISSDP) :: RHO_PARAM=1.0_CMISSDP
  REAL(CMISSDP) :: MU_PARAM=1.0_CMISSDP

  !CMISS variables

  TYPE(CMISSBasisType) :: BasisSpace1,BasisSpace2
  TYPE(CMISSBasisType) :: BasisVelocity1,BasisVelocity2
  TYPE(CMISSBasisType) :: BasisPressure1,BasisPressure2
  TYPE(CMISSBasisType) :: InterfaceBasis,InterfaceMappingBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem1,CoordinateSystem2,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition1,Decomposition2,InterfaceDecomposition
  TYPE(CMISSEquationsType) :: Equations1,Equations2
  TYPE(CMISSEquationsSetType) :: EquationsSet1,EquationsSet2
  TYPE(CMISSFieldType) :: GeometricField1,GeometricField2
  TYPE(CMISSFieldType) :: InterfaceGeometricField,LagrangeField
  TYPE(CMISSFieldType) :: DependentField1,DependentField2
  TYPE(CMISSFieldType) :: EquationsSetField1,EquationsSetField2
  TYPE(CMISSFieldType) :: MaterialsField1,MaterialsField2

  TYPE(CMISSFieldsType) :: Fields1,Fields2,InterfaceFields
  TYPE(CMISSInterfaceType) :: Interface
  TYPE(CMISSInterfaceConditionType) :: InterfaceCondition
  TYPE(CMISSInterfaceEquationsType) :: InterfaceEquations
  TYPE(CMISSInterfaceMeshConnectivityType) :: InterfaceMeshConnectivity
  !Nodes
  TYPE(CMISSNodesType) :: Nodes1
  TYPE(CMISSNodesType) :: Nodes2
  TYPE(CMISSNodesType) :: InterfaceNodes
  !Elements
  TYPE(CMISSMeshElementsType) :: MeshElementsSpace1,MeshElementsVelocity1,MeshElementsPressure1
  TYPE(CMISSMeshElementsType) :: MeshElementsSpace2,MeshElementsVelocity2,MeshElementsPressure2
  TYPE(CMISSMeshElementsType) :: InterfaceMeshElements
  TYPE(CMISSMeshType) :: Mesh1,Mesh2,InterfaceMesh
  TYPE(CMISSProblemType) :: CoupledProblem
  TYPE(CMISSRegionType) :: Region1,Region2,WorldRegion
  TYPE(CMISSSolverType) :: CoupledSolver
  TYPE(CMISSSolverEquationsType) :: CoupledSolverEquations
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
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

  !
  !================================================================================================================================
  !

  !INITIALISE OPENCMISS


  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Set error handling mode
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
 
  !Set diganostics for testing
  !CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["SOLVER_MAPPING_CALCULATE         ", &
  !  & "SOLVER_MATRIX_STRUCTURE_CALCULATE"],Err)
  
  !
  !================================================================================================================================
  !

  !CHECK COMPUTATIONAL NODE

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !
  !================================================================================================================================
  !

  !PROBLEM CONTROL PANEL

  !Import cmHeart mesh information
  CALL FLUID_MECHANICS_IO_READ_CMHEART(CM1,CM2,CM3,CMX,BC,Err)

  !Information for mesh 1
  BASIS_NUMBER_SPACE1=CM1%ID_M
  BASIS_NUMBER_VELOCITY1=CM1%ID_V
  BASIS_NUMBER_PRESSURE1=CM1%ID_P
  NUMBER_OF_DIMENSIONS1=CM1%D
  BASIS_TYPE1=CM1%IT_T
  BASIS_XI_INTERPOLATION_SPACE1=CM1%IT_M
  BASIS_XI_INTERPOLATION_VELOCITY1=CM1%IT_V
  BASIS_XI_INTERPOLATION_PRESSURE1=CM1%IT_P
  NUMBER_OF_NODES_SPACE1=CM1%N_M
  NUMBER_OF_NODES_VELOCITY1=CM1%N_V
  NUMBER_OF_NODES_PRESSURE1=CM1%N_P
  TOTAL_NUMBER_OF_NODES1=CM1%N_T
  TOTAL_NUMBER_OF_ELEMENTS1=CM1%E_T
  NUMBER_OF_ELEMENT_NODES_SPACE1=CM1%EN_M
  NUMBER_OF_ELEMENT_NODES_VELOCITY1=CM1%EN_V
  NUMBER_OF_ELEMENT_NODES_PRESSURE1=CM1%EN_P
  !Information for mesh 2
  BASIS_NUMBER_SPACE2=CM2%ID_M+10
  BASIS_NUMBER_VELOCITY2=CM2%ID_V+10
  BASIS_NUMBER_PRESSURE2=CM2%ID_P+10
  NUMBER_OF_DIMENSIONS2=CM2%D
  BASIS_TYPE2=CM2%IT_T
  BASIS_XI_INTERPOLATION_SPACE2=CM2%IT_M
  BASIS_XI_INTERPOLATION_VELOCITY2=CM2%IT_V
  BASIS_XI_INTERPOLATION_PRESSURE2=CM2%IT_P
  NUMBER_OF_NODES_SPACE2=CM2%N_M
  NUMBER_OF_NODES_VELOCITY2=CM2%N_V
  NUMBER_OF_NODES_PRESSURE2=CM2%N_P
  TOTAL_NUMBER_OF_NODES2=CM2%N_T
  TOTAL_NUMBER_OF_ELEMENTS2=CM2%E_T
  NUMBER_OF_ELEMENT_NODES_SPACE2=CM2%EN_M
  NUMBER_OF_ELEMENT_NODES_VELOCITY2=CM2%EN_V
  NUMBER_OF_ELEMENT_NODES_PRESSURE2=CM2%EN_P
  !Set interpolation parameters
  BASIS_XI_GAUSS_SPACE1=3
  BASIS_XI_GAUSS_VELOCITY1=3
  BASIS_XI_GAUSS_PRESSURE1=3
  BASIS_XI_GAUSS_SPACE2=3
  BASIS_XI_GAUSS_VELOCITY2=3
  BASIS_XI_GAUSS_PRESSURE2=3
  BASIS_NUMBER_SPACE_INTERFACE=CM3%ID_M
  NUMBER_OF_DIMENSIONS_INTERFACE=CM3%D
  BASIS_TYPE_INTERFACE=CM3%IT_T
  BASIS_XI_INTERPOLATION_INTERFACE=CM3%IT_M
  BASIS_XI_INTERPOLATION_VELOCITY2=CM2%IT_V
  BASIS_XI_INTERPOLATION_PRESSURE2=CM2%IT_P
  NUMBER_OF_NODES_INTERFACE=CM3%N_M
  NUMBER_OF_NODES_VELOCITY2=CM2%N_V
  NUMBER_OF_NODES_PRESSURE2=CM2%N_P
  TOTAL_NUMBER_OF_NODES_INTERFACE=CM3%N_T
  TOTAL_NUMBER_OF_ELEMENTS_INTERFACE=CM3%E_T
  NUMBER_OF_ELEMENT_NODES_INTERFACE=CM3%EN_M
  NUMBER_OF_ELEMENT_NODES_VELOCITY2=CM2%EN_V
  NUMBER_OF_ELEMENT_NODES_PRESSURE2=CM2%EN_P
  !Set interpolation parameters
  BASIS_XI_GAUSS_INTERFACE=3
  BASIS_XI_GAUSS_VELOCITY1=3
  BASIS_XI_GAUSS_PRESSURE1=3
  BASIS_XI_GAUSS_INTERFACE=3
  BASIS_XI_GAUSS_VELOCITY2=3
  BASIS_XI_GAUSS_PRESSURE2=3

  IF(FIXED_WALL_NODES1_FLAG) THEN
    NUMBER_OF_FIXED_WALL_NODES1=80
    ALLOCATE(FIXED_WALL_NODES1(NUMBER_OF_FIXED_WALL_NODES1))
    FIXED_WALL_NODES1=(/1,2,3,4,5,7,9,10,11,12,13,14,17,20,24,28,29,30,31,32,33,34,35,37,39, & 
    & 41,44,46,47,48,50,51,52,53,54,57,60,64,65,66,67,68,70,72,74,76,77,78,79,80,83,86, & 
    & 89,90,91,92,93,94,95,97,99,101,102,103,104,105,106,107,108,111,114,115,116,117,118, & 
    & 120,122,123,124,125/)
  ENDIF
  IF(INLET_WALL_NODES1_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES1=9
    ALLOCATE(INLET_WALL_NODES1(NUMBER_OF_INLET_WALL_NODES1))
    INLET_WALL_NODES1=(/67,74,76,79,87,88,102,103,113/)
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS1(1)=1.0_CMISSDP
    BOUNDARY_CONDITIONS1(2)=0.0_CMISSDP
    BOUNDARY_CONDITIONS1(3)=0.0_CMISSDP
  ENDIF
  IF(FIXED_WALL_NODES2_FLAG) THEN
    NUMBER_OF_FIXED_WALL_NODES2=80
    ALLOCATE(FIXED_WALL_NODES2(NUMBER_OF_FIXED_WALL_NODES2))
    FIXED_WALL_NODES2=(/1,2,3,4,5,7,9,10,11,12,13,14,17,20,24,28,29,30,31,32,33,34,35,37,39, & 
    & 41,44,46,47,48,50,51,52,53,54,57,60,64,65,66,67,68,70,72,74,76,77,78,79,80,83,86, & 
    & 89,90,91,92,93,94,95,97,99,101,102,103,104,105,106,107,108,111,114,115,116,117,118, & 
    & 120,122,123,124,125/)
  ENDIF
  IF(INLET_WALL_NODES2_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES2=9
    ALLOCATE(INLET_WALL_NODES2(NUMBER_OF_INLET_WALL_NODES1))
    INLET_WALL_NODES2=(/67,74,76,79,87,88,102,103,113/)
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS2(1)=1.0_CMISSDP
    BOUNDARY_CONDITIONS2(2)=0.0_CMISSDP
    BOUNDARY_CONDITIONS2(3)=0.0_CMISSDP
  ENDIF



  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  
  !Start the creation of a new RC coordinate system for the first region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM(1) << == '
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem1,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystem1UserNumber,CoordinateSystem1,Err)
  !Set the coordinate system dimension
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem1,NUMBER_OF_DIMENSIONS1,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem1,Err)

  !Start the creation of a new RC coordinate system for the second region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM(2) << == '
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem2,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystem2UserNumber,CoordinateSystem2,Err)
  !Set the coordinate system dimension
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem2,NUMBER_OF_DIMENSIONS2,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem2,Err)

  !
  !================================================================================================================================
  !

  !REGION
  
  !Start the creation of the first region
  PRINT *, ' == >> CREATING REGION(1) << == '
  CALL CMISSRegion_Initialise(Region1,Err)
  CALL CMISSRegion_CreateStart(Region1UserNumber,WorldRegion,Region1,Err)
  CALL CMISSRegion_LabelSet(Region1,"Region1",Err)
  !Set the regions coordinate system as defined above
  CALL CMISSRegion_CoordinateSystemSet(Region1,CoordinateSystem1,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region1,Err)

  !Start the creation of the second region
  PRINT *, ' == >> CREATING REGION(2) << == '
  CALL CMISSRegion_Initialise(Region2,Err)
  CALL CMISSRegion_CreateStart(Region2UserNumber,WorldRegion,Region2,Err)
  CALL CMISSRegion_LabelSet(Region2,"Region2",Err)
  !Set the regions coordinate system as defined above
  CALL CMISSRegion_CoordinateSystemSet(Region2,CoordinateSystem2,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region2,Err)

  !
  !================================================================================================================================
  !

  !BASES


  !Start the creation of a bI/tri-linear-Lagrange basis
  PRINT *, ' == >> CREATING BASIS(1) << == '
  !Start the creation of new bases
  MESH_NUMBER_OF_COMPONENTS1=1
  CALL CMISSBasis_Initialise(BasisSpace1,Err)
  CALL CMISSBasis_CreateStart(BASIS_NUMBER_SPACE1,BasisSpace1,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasis_TypeSet(BasisSpace1,BASIS_TYPE1,Err)
  !Set the basis xi number
  CALL CMISSBasis_NumberOfXiSet(BasisSpace1,NUMBER_OF_DIMENSIONS1,Err)
  !Set the basis xi interpolation and number of Gauss points
  IF(NUMBER_OF_DIMENSIONS1==2.AND.NUMBER_OF_DIMENSIONS2==2) THEN
    CALL CMISSBasis_InterpolationXiSet(BasisSpace1,(/BASIS_XI_INTERPOLATION_SPACE1,BASIS_XI_INTERPOLATION_SPACE1/),Err)
    IF(BASIS_TYPE1/=CMISS_BASIS_SIMPLEX_TYPE) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisSpace1,(/BASIS_XI_GAUSS_SPACE1,BASIS_XI_GAUSS_SPACE1/),Err)
    ELSE
      CALL CMISSBasis_QuadratureOrderSet(BasisSpace1,BASIS_XI_GAUSS_SPACE1+1,Err)
    ENDIF
  ELSE IF(NUMBER_OF_DIMENSIONS1==3.AND.NUMBER_OF_DIMENSIONS2==3) THEN
    CALL CMISSBasis_InterpolationXiSet(BasisSpace1,(/BASIS_XI_INTERPOLATION_SPACE1,BASIS_XI_INTERPOLATION_SPACE1, & 
      & BASIS_XI_INTERPOLATION_SPACE1/),Err)                         
    IF(BASIS_TYPE1/=CMISS_BASIS_SIMPLEX_TYPE) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisSpace1,(/BASIS_XI_GAUSS_SPACE1,BASIS_XI_GAUSS_SPACE1, &
        & BASIS_XI_GAUSS_SPACE1/),Err)
    ELSE
      CALL CMISSBasis_QuadratureOrderSet(BasisSpace1,BASIS_XI_GAUSS_SPACE1+1,Err)
    ENDIF
  ELSE
    CALL HANDLE_ERROR("Dimension coupling error.")
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(BasisSpace1,Err)
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_VELOCITY1==BASIS_XI_INTERPOLATION_SPACE1) THEN
    BasisVelocity1=BasisSpace1
  ELSE
    MESH_NUMBER_OF_COMPONENTS1=MESH_NUMBER_OF_COMPONENTS1+1
    !Initialise a new velocity basis
    CALL CMISSBasis_Initialise(BasisVelocity1,Err)
    !Start the creation of a basis
    CALL CMISSBasis_CreateStart(BASIS_NUMBER_VELOCITY1,BasisVelocity1,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasis_TypeSet(BasisVelocity1,BASIS_TYPE1,Err)
    !Set the basis xi number
    CALL CMISSBasis_NumberOfXiSet(BasisVelocity1,NUMBER_OF_DIMENSIONS1,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS1==2.AND.NUMBER_OF_DIMENSIONS2==2) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisVelocity1,(/BASIS_XI_INTERPOLATION_VELOCITY1,BASIS_XI_INTERPOLATION_VELOCITY1/),Err)
      IF(BASIS_TYPE1/=CMISS_BASIS_SIMPLEX_TYPE) THEN
        CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisVelocity1,(/BASIS_XI_GAUSS_VELOCITY1,BASIS_XI_GAUSS_VELOCITY1/),Err)
      ELSE
        CALL CMISSBasis_QuadratureOrderSet(BasisVelocity1,BASIS_XI_GAUSS_VELOCITY1+1,Err)
      ENDIF
    ELSE IF(NUMBER_OF_DIMENSIONS1==3.AND.NUMBER_OF_DIMENSIONS2==3) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisVelocity1,(/BASIS_XI_INTERPOLATION_VELOCITY1,BASIS_XI_INTERPOLATION_VELOCITY1, & 
        & BASIS_XI_INTERPOLATION_VELOCITY1/),Err)                         
      IF(BASIS_TYPE1/=CMISS_BASIS_SIMPLEX_TYPE) THEN
        CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisVelocity1,(/BASIS_XI_GAUSS_VELOCITY1,BASIS_XI_GAUSS_VELOCITY1, &
          & BASIS_XI_GAUSS_VELOCITY1/),Err)
      ELSE
        CALL CMISSBasis_QuadratureOrderSet(BasisVelocity1,BASIS_XI_GAUSS_VELOCITY1+1,Err)
      ENDIF
    ELSE
      CALL HANDLE_ERROR("Dimension coupling error.")
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasis_CreateFinish(BasisVelocity1,Err)
  ENDIF
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_PRESSURE1==BASIS_XI_INTERPOLATION_SPACE1) THEN
    BasisPressure1=BasisSpace1
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE1==BASIS_XI_INTERPOLATION_VELOCITY1) THEN
    BasisPressure1=BasisVelocity1
  ELSE
    MESH_NUMBER_OF_COMPONENTS1=MESH_NUMBER_OF_COMPONENTS1+1
    !Initialise a new pressure basis
    CALL CMISSBasis_Initialise(BasisPressure1,Err)
    !Start the creation of a basis
    CALL CMISSBasis_CreateStart(BASIS_NUMBER_PRESSURE1,BasisPressure1,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasis_TypeSet(BasisPressure1,BASIS_TYPE1,Err)
    !Set the basis xi number
    CALL CMISSBasis_NumberOfXiSet(BasisPressure1,NUMBER_OF_DIMENSIONS1,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS1==2.AND.NUMBER_OF_DIMENSIONS2==2) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisPressure1,(/BASIS_XI_INTERPOLATION_PRESSURE1,BASIS_XI_INTERPOLATION_PRESSURE1/),Err)
      IF(BASIS_TYPE1/=CMISS_BASIS_SIMPLEX_TYPE) THEN
        CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisPressure1,(/BASIS_XI_GAUSS_PRESSURE1,BASIS_XI_GAUSS_PRESSURE1/),Err)
      ELSE
        CALL CMISSBasis_QuadratureOrderSet(BasisPressure1,BASIS_XI_GAUSS_PRESSURE1+1,Err)
      ENDIF
    ELSE IF(NUMBER_OF_DIMENSIONS1==3.AND.NUMBER_OF_DIMENSIONS2==3) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisPressure1,(/BASIS_XI_INTERPOLATION_PRESSURE1,BASIS_XI_INTERPOLATION_PRESSURE1, & 
        & BASIS_XI_INTERPOLATION_PRESSURE1/),Err)                         
      IF(BASIS_TYPE1/=CMISS_BASIS_SIMPLEX_TYPE) THEN
        CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisPressure1,(/BASIS_XI_GAUSS_PRESSURE1,BASIS_XI_GAUSS_PRESSURE1, &
          & BASIS_XI_GAUSS_PRESSURE1/),Err)
      ELSE
        CALL CMISSBasis_QuadratureOrderSet(BasisPressure1,BASIS_XI_GAUSS_PRESSURE1+1,Err)
      ENDIF
    ELSE
      CALL HANDLE_ERROR("Dimension coupling error.")
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasis_CreateFinish(BasisPressure1,Err)
  ENDIF


  !Start the creation of a bI/tri-XXX-Lagrange basis
  PRINT *, ' == >> CREATING BASIS(2) << == '
  !Start the creation of new bases
  MESH_NUMBER_OF_COMPONENTS2=1
  CALL CMISSBasis_Initialise(BasisSpace2,Err)
  CALL CMISSBasis_CreateStart(BASIS_NUMBER_SPACE2,BasisSpace2,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasis_TypeSet(BasisSpace2,BASIS_TYPE2,Err)
  !Set the basis xi number
  CALL CMISSBasis_NumberOfXiSet(BasisSpace2,NUMBER_OF_DIMENSIONS2,Err)
  !Set the basis xi interpolation and number of Gauss points
  IF(NUMBER_OF_DIMENSIONS1==2.AND.NUMBER_OF_DIMENSIONS2==2) THEN
    CALL CMISSBasis_InterpolationXiSet(BasisSpace2,(/BASIS_XI_INTERPOLATION_SPACE2,BASIS_XI_INTERPOLATION_SPACE2/),Err)
    IF(BASIS_TYPE2/=CMISS_BASIS_SIMPLEX_TYPE) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisSpace2,(/BASIS_XI_GAUSS_SPACE2,BASIS_XI_GAUSS_SPACE2/),Err)
    ELSE
      CALL CMISSBasis_QuadratureOrderSet(BasisSpace2,BASIS_XI_GAUSS_SPACE2+1,Err)
    ENDIF
  ELSE IF(NUMBER_OF_DIMENSIONS1==3.AND.NUMBER_OF_DIMENSIONS2==3) THEN
    CALL CMISSBasis_InterpolationXiSet(BasisSpace2,(/BASIS_XI_INTERPOLATION_SPACE2,BASIS_XI_INTERPOLATION_SPACE2, & 
      & BASIS_XI_INTERPOLATION_SPACE2/),Err)                         
    IF(BASIS_TYPE2/=CMISS_BASIS_SIMPLEX_TYPE) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisSpace2,(/BASIS_XI_GAUSS_SPACE2,BASIS_XI_GAUSS_SPACE2, &
        & BASIS_XI_GAUSS_SPACE2/),Err)
    ELSE
      CALL CMISSBasis_QuadratureOrderSet(BasisSpace2,BASIS_XI_GAUSS_SPACE2+1,Err)
    ENDIF
  ELSE
    CALL HANDLE_ERROR("Dimension coupling error.")
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(BasisSpace2,Err)
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_VELOCITY2==BASIS_XI_INTERPOLATION_SPACE2) THEN
    BasisVelocity2=BasisSpace2
  ELSE
    MESH_NUMBER_OF_COMPONENTS2=MESH_NUMBER_OF_COMPONENTS2+1
    !Initialise a new velocity basis
    CALL CMISSBasis_Initialise(BasisVelocity2,Err)
    !Start the creation of a basis
    CALL CMISSBasis_CreateStart(BASIS_NUMBER_VELOCITY2,BasisVelocity2,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasis_TypeSet(BasisVelocity2,BASIS_TYPE2,Err)
    !Set the basis xi number
    CALL CMISSBasis_NumberOfXiSet(BasisVelocity2,NUMBER_OF_DIMENSIONS2,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS1==2.AND.NUMBER_OF_DIMENSIONS2==2) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisVelocity2,(/BASIS_XI_INTERPOLATION_VELOCITY2,BASIS_XI_INTERPOLATION_VELOCITY2/),Err)
      IF(BASIS_TYPE2/=CMISS_BASIS_SIMPLEX_TYPE) THEN
        CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisVelocity2,(/BASIS_XI_GAUSS_VELOCITY2,BASIS_XI_GAUSS_VELOCITY2/),Err)
      ELSE
        CALL CMISSBasis_QuadratureOrderSet(BasisVelocity2,BASIS_XI_GAUSS_VELOCITY2+1,Err)
      ENDIF
    ELSE IF(NUMBER_OF_DIMENSIONS1==3.AND.NUMBER_OF_DIMENSIONS2==3) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisVelocity2,(/BASIS_XI_INTERPOLATION_VELOCITY2,BASIS_XI_INTERPOLATION_VELOCITY2, & 
        & BASIS_XI_INTERPOLATION_VELOCITY2/),Err)                         
      IF(BASIS_TYPE2/=CMISS_BASIS_SIMPLEX_TYPE) THEN
        CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisVelocity2,(/BASIS_XI_GAUSS_VELOCITY2,BASIS_XI_GAUSS_VELOCITY2, &
          & BASIS_XI_GAUSS_VELOCITY2/),Err)
      ELSE
        CALL CMISSBasis_QuadratureOrderSet(BasisVelocity2,BASIS_XI_GAUSS_VELOCITY2+1,Err)
      ENDIF
    ELSE
      CALL HANDLE_ERROR("Dimension coupling error.")
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasis_CreateFinish(BasisVelocity2,Err)
  ENDIF
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_PRESSURE2==BASIS_XI_INTERPOLATION_SPACE2) THEN
    BasisPressure2=BasisSpace2
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE2==BASIS_XI_INTERPOLATION_VELOCITY2) THEN
    BasisPressure2=BasisVelocity2
  ELSE
    MESH_NUMBER_OF_COMPONENTS2=MESH_NUMBER_OF_COMPONENTS2+1
    !Initialise a new pressure basis
    CALL CMISSBasis_Initialise(BasisPressure2,Err)
    !Start the creation of a basis
    CALL CMISSBasis_CreateStart(BASIS_NUMBER_PRESSURE2,BasisPressure2,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasis_TypeSet(BasisPressure2,BASIS_TYPE2,Err)
    !Set the basis xi number
    CALL CMISSBasis_NumberOfXiSet(BasisPressure2,NUMBER_OF_DIMENSIONS2,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS1==2.AND.NUMBER_OF_DIMENSIONS2==2) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisPressure2,(/BASIS_XI_INTERPOLATION_PRESSURE2,BASIS_XI_INTERPOLATION_PRESSURE2/),Err)
      IF(BASIS_TYPE2/=CMISS_BASIS_SIMPLEX_TYPE) THEN
        CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisPressure2,(/BASIS_XI_GAUSS_PRESSURE2,BASIS_XI_GAUSS_PRESSURE2/),Err)
      ELSE
        CALL CMISSBasis_QuadratureOrderSet(BasisPressure2,BASIS_XI_GAUSS_PRESSURE2+1,Err)
      ENDIF
    ELSE IF(NUMBER_OF_DIMENSIONS1==3.AND.NUMBER_OF_DIMENSIONS2==3) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisPressure2,(/BASIS_XI_INTERPOLATION_PRESSURE2,BASIS_XI_INTERPOLATION_PRESSURE2, & 
        & BASIS_XI_INTERPOLATION_PRESSURE2/),Err)                         
      IF(BASIS_TYPE2/=CMISS_BASIS_SIMPLEX_TYPE) THEN
        CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisPressure2,(/BASIS_XI_GAUSS_PRESSURE2,BASIS_XI_GAUSS_PRESSURE2, &
          & BASIS_XI_GAUSS_PRESSURE2/),Err)
      ELSE
        CALL CMISSBasis_QuadratureOrderSet(BasisPressure2,BASIS_XI_GAUSS_PRESSURE2+1,Err)
      ENDIF
    ELSE
      CALL HANDLE_ERROR("Dimension coupling error.")
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasis_CreateFinish(BasisPressure2,Err)
  ENDIF


  !
  !================================================================================================================================
  !

  !MESH

  
  !Start the creation of a generated mesh in the first region
  PRINT *, ' == >> CREATING MESH(1) FROM INPUT DATA << == '
 !Start the creation of mesh nodes
  CALL CMISSNodes_Initialise(Nodes1,Err)
  CALL CMISSMesh_Initialise(Mesh1,Err)
  CALL CMISSNodes_CreateStart(Region1,TOTAL_NUMBER_OF_NODES1,Nodes1,Err)
  CALL CMISSNodes_CreateFinish(Nodes1,Err)
  !Start the creation of the mesh
  CALL CMISSMesh_CreateStart(Mesh1UserNumber,Region1,NUMBER_OF_DIMENSIONS1,Mesh1,Err)
  !Set number of mesh elements
  CALL CMISSMesh_NumberOfElementsSet(Mesh1,TOTAL_NUMBER_OF_ELEMENTS1,Err)
  !Set number of mesh components
  CALL CMISSMesh_NumberOfComponentsSet(Mesh1,MESH_NUMBER_OF_COMPONENTS1,Err)
  !Specify spatial mesh component
  CALL CMISSMeshElements_Initialise(MeshElementsSpace1,Err)
  CALL CMISSMeshElements_Initialise(MeshElementsVelocity1,Err)
  CALL CMISSMeshElements_Initialise(MeshElementsPressure1,Err)
  MESH_COMPONENT_NUMBER_SPACE1=1
  MESH_COMPONENT_NUMBER_VELOCITY1=1
  MESH_COMPONENT_NUMBER_PRESSURE1=1
  CALL CMISSMeshElements_CreateStart(Mesh1,MESH_COMPONENT_NUMBER_SPACE1,BasisSpace1,MeshElementsSpace1,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS1
    CALL CMISSMeshElements_NodesSet(MeshElementsSpace1,ELEMENT_NUMBER,CM1%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_SPACE1),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(MeshElementsSpace1,Err)
  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY1==BASIS_XI_INTERPOLATION_SPACE1) THEN
    MeshElementsVelocity1=MeshElementsSpace1
  ELSE
    MESH_COMPONENT_NUMBER_VELOCITY1=MESH_COMPONENT_NUMBER_SPACE1+1
    CALL CMISSMeshElements_CreateStart(Mesh1,MESH_COMPONENT_NUMBER_VELOCITY1,BasisVelocity1,MeshElementsVelocity1,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS1
      CALL CMISSMeshElements_NodesSet(MeshElementsVelocity1,ELEMENT_NUMBER,CM1%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_VELOCITY1),Err)
    ENDDO
    CALL CMISSMeshElements_CreateFinish(MeshElementsVelocity1,Err)
  ENDIF
  !Specify pressure mesh component
  IF(BASIS_XI_INTERPOLATION_PRESSURE1==BASIS_XI_INTERPOLATION_SPACE1) THEN
    MeshElementsPressure1=MeshElementsSpace1
    MESH_COMPONENT_NUMBER_PRESSURE1=MESH_COMPONENT_NUMBER_SPACE1
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE1==BASIS_XI_INTERPOLATION_VELOCITY1) THEN
    MeshElementsPressure1=MeshElementsVelocity1
    MESH_COMPONENT_NUMBER_PRESSURE1=MESH_COMPONENT_NUMBER_VELOCITY1
  ELSE
    MESH_COMPONENT_NUMBER_PRESSURE1=MESH_COMPONENT_NUMBER_VELOCITY1+1
    CALL CMISSMeshElements_CreateStart(Mesh1,MESH_COMPONENT_NUMBER_PRESSURE1,BasisPressure1,MeshElementsPressure1,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS1
      CALL CMISSMeshElements_NodesSet(MeshElementsPressure1,ELEMENT_NUMBER,CM1%P(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_PRESSURE1),Err)
    ENDDO
    CALL CMISSMeshElements_CreateFinish(MeshElementsPressure1,Err)
  ENDIF
  !Finish the creation of the mesh
  CALL CMISSMesh_CreateFinish(Mesh1,Err)

  !Start the creation of a generated mesh in the second region
  PRINT *, ' == >> CREATING MESH(2) FROM INPUT DATA << == '
  !Start the creation of mesh nodes
  CALL CMISSNodes_Initialise(Nodes2,Err)
  CALL CMISSMesh_Initialise(Mesh2,Err)
  CALL CMISSNodes_CreateStart(Region2,TOTAL_NUMBER_OF_NODES2,Nodes2,Err)
  CALL CMISSNodes_CreateFinish(Nodes2,Err)
  !Start the creation of the mesh
  CALL CMISSMesh_CreateStart(Mesh2UserNumber,Region2,NUMBER_OF_DIMENSIONS2,Mesh2,Err)
  !Set number of mesh elements
  CALL CMISSMesh_NumberOfElementsSet(Mesh2,TOTAL_NUMBER_OF_ELEMENTS2,Err)
  !Set number of mesh components
  CALL CMISSMesh_NumberOfComponentsSet(Mesh2,MESH_NUMBER_OF_COMPONENTS2,Err)
  !Specify spatial mesh component
  CALL CMISSMeshElements_Initialise(MeshElementsSpace2,Err)
  CALL CMISSMeshElements_Initialise(MeshElementsVelocity2,Err)
  CALL CMISSMeshElements_Initialise(MeshElementsPressure2,Err)
  MESH_COMPONENT_NUMBER_SPACE2=1
  MESH_COMPONENT_NUMBER_VELOCITY2=1
  MESH_COMPONENT_NUMBER_PRESSURE2=1
  CALL CMISSMeshElements_CreateStart(Mesh2,MESH_COMPONENT_NUMBER_SPACE2,BasisSpace2,MeshElementsSpace2,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS2
    CALL CMISSMeshElements_NodesSet(MeshElementsSpace2,ELEMENT_NUMBER,CM2%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_SPACE2),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(MeshElementsSpace2,Err)
  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY2==BASIS_XI_INTERPOLATION_SPACE2) THEN
    MeshElementsVelocity2=MeshElementsSpace2
  ELSE
    MESH_COMPONENT_NUMBER_VELOCITY2=MESH_COMPONENT_NUMBER_SPACE2+1
    CALL CMISSMeshElements_CreateStart(Mesh2,MESH_COMPONENT_NUMBER_VELOCITY2,BasisVelocity2,MeshElementsVelocity2,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS2
      CALL CMISSMeshElements_NodesSet(MeshElementsVelocity2,ELEMENT_NUMBER,CM2%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_VELOCITY2),Err)
    ENDDO
    CALL CMISSMeshElements_CreateFinish(MeshElementsVelocity2,Err)
  ENDIF
  !Specify pressure mesh component
  IF(BASIS_XI_INTERPOLATION_PRESSURE2==BASIS_XI_INTERPOLATION_SPACE2) THEN
    MeshElementsPressure2=MeshElementsSpace2
    MESH_COMPONENT_NUMBER_PRESSURE2=MESH_COMPONENT_NUMBER_SPACE2
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE2==BASIS_XI_INTERPOLATION_VELOCITY2) THEN
    MeshElementsPressure2=MeshElementsVelocity2
    MESH_COMPONENT_NUMBER_PRESSURE2=MESH_COMPONENT_NUMBER_VELOCITY2
  ELSE
    MESH_COMPONENT_NUMBER_PRESSURE2=MESH_COMPONENT_NUMBER_VELOCITY2+1
    CALL CMISSMeshElements_CreateStart(Mesh2,MESH_COMPONENT_NUMBER_PRESSURE2,BasisPressure2,MeshElementsPressure2,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS2
      CALL CMISSMeshElements_NodesSet(MeshElementsPressure2,ELEMENT_NUMBER,CM2%P(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_PRESSURE2),Err)
    ENDDO
    CALL CMISSMeshElements_CreateFinish(MeshElementsPressure2,Err)
  ENDIF
  !Finish the creation of the mesh
  CALL CMISSMesh_CreateFinish(Mesh2,Err)


  !
  !================================================================================================================================
  !

  !INTERFACE DEFINITION

  !Create an interface between the two meshes
  PRINT *, ' == >> CREATING INTERFACE << == '
  CALL CMISSInterface_Initialise(Interface,Err)
  CALL CMISSInterface_CreateStart(InterfaceUserNumber,WorldRegion,Interface,Err)
  CALL CMISSInterface_LabelSet(Interface,"Interface",Err)
  !Add in the two meshes
  CALL CMISSInterface_MeshAdd(Interface,Mesh1,Mesh1Index,Err)
  CALL CMISSInterface_MeshAdd(Interface,Mesh2,Mesh2Index,Err)
  !Finish creating the interface
  CALL CMISSInterface_CreateFinish(Interface,Err)

  !Start the creation of a (bi)-linear-Lagrange basis
  PRINT *, ' == >> CREATING INTERFACE BASIS << == '
  CALL CMISSBasis_Initialise(InterfaceBasis,Err)
  CALL CMISSBasis_CreateStart(InterfaceBasisUserNumber,InterfaceBasis,Err)
  CALL CMISSBasis_TypeSet(InterfaceBasis,BASIS_TYPE_INTERFACE,Err)
  CALL CMISSBasis_NumberOfXiSet(InterfaceBasis,NUMBER_OF_DIMENSIONS_INTERFACE,Err)
  !Set the basis xi interpolation and number of Gauss points
  IF(NUMBER_OF_DIMENSIONS1==3.AND.NUMBER_OF_DIMENSIONS2==3.AND.NUMBER_OF_DIMENSIONS_INTERFACE==2) THEN
    CALL CMISSBasis_InterpolationXiSet(InterfaceBasis,(/BASIS_XI_INTERPOLATION_INTERFACE,BASIS_XI_INTERPOLATION_INTERFACE/),Err)
    IF(BASIS_TYPE_INTERFACE/=CMISS_BASIS_SIMPLEX_TYPE) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(InterfaceBasis,(/BASIS_XI_GAUSS_INTERFACE,BASIS_XI_GAUSS_INTERFACE/),Err)
    ELSE
      CALL CMISSBasis_QuadratureOrderSet(InterfaceBasis,BASIS_XI_GAUSS_INTERFACE+1,Err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(InterfaceBasis,Err)


  !
  !================================================================================================================================
  !

  !INTERFACE MAPPING

  !Start the creation of a (bi)-linear-Lagrange basis
  PRINT *, ' == >> CREATING INTERFACE MAPPING BASIS << == '
  CALL CMISSBasis_Initialise(InterfaceMappingBasis,Err)
  CALL CMISSBasis_CreateStart(InterfaceMappingBasisUserNumber,InterfaceMappingBasis,Err)
  CALL CMISSBasis_TypeSet(InterfaceMappingBasis,BASIS_TYPE_INTERFACE,Err)
  CALL CMISSBasis_NumberOfXiSet(InterfaceMappingBasis,NUMBER_OF_DIMENSIONS_INTERFACE,Err)
  IF(NUMBER_OF_DIMENSIONS1==3.AND.NUMBER_OF_DIMENSIONS2==3.AND.NUMBER_OF_DIMENSIONS_INTERFACE==2) THEN
    CALL CMISSBasis_InterpolationXiSet(InterfaceMappingBasis,(/BASIS_XI_INTERPOLATION_INTERFACE, &
      & BASIS_XI_INTERPOLATION_INTERFACE/),Err)
    IF(BASIS_TYPE_INTERFACE/=CMISS_BASIS_SIMPLEX_TYPE) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(InterfaceMappingBasis,(/BASIS_XI_GAUSS_INTERFACE,BASIS_XI_GAUSS_INTERFACE/),Err)
    ELSE
      CALL CMISSBasis_QuadratureOrderSet(InterfaceMappingBasis,BASIS_XI_GAUSS_INTERFACE+1,Err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(InterfaceMappingBasis,Err)

  !
  !================================================================================================================================
  !

  !INTERFACE MESH
  
  !Start the creation of a generated mesh for the interface
  PRINT *, ' == >> CREATING INTERFACE MESH FROM INPUT DATA << == '
  !Start the creation of mesh nodes
  CALL CMISSNodes_Initialise(InterfaceNodes,Err)
  CALL CMISSMesh_Initialise(InterfaceMesh,Err)
  CALL CMISSNodes_CreateStart(Interface,TOTAL_NUMBER_OF_NODES_INTERFACE,InterfaceNodes,Err)
  CALL CMISSNodes_CreateFinish(InterfaceNodes,Err)
  !Start the creation of the mesh
  CALL CMISSMesh_CreateStart(InterfaceMeshUserNumber,Interface,NUMBER_OF_DIMENSIONS_INTERFACE,InterfaceMesh,Err)
  !Set number of mesh elements
  CALL CMISSMesh_NumberOfElementsSet(InterfaceMesh,TOTAL_NUMBER_OF_ELEMENTS_INTERFACE,Err)
  MESH_NUMBER_OF_COMPONENTS_INTERFACE=1
  !Set number of mesh components
  CALL CMISSMesh_NumberOfComponentsSet(InterfaceMesh,MESH_NUMBER_OF_COMPONENTS_INTERFACE,Err)
  !Specify spatial mesh component
  CALL CMISSMeshElements_Initialise(InterfaceMeshElements,Err)
  MESH_COMPONENT_NUMBER_INTERFACE=1
  CALL CMISSMeshElements_CreateStart(InterfaceMesh,MESH_COMPONENT_NUMBER_INTERFACE,InterfaceBasis,InterfaceMeshElements,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS_INTERFACE
    CALL CMISSMeshElements_NodesSet(InterfaceMeshElements,ELEMENT_NUMBER,CM3%M(ELEMENT_NUMBER, &
      & 1:NUMBER_OF_ELEMENT_NODES_INTERFACE),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(InterfaceMeshElements,Err)
  !Finish the creation of the mesh
  CALL CMISSMesh_CreateFinish(InterfaceMesh,Err)


  !
  !================================================================================================================================
  !

  !INTERFACE CONNECTIVITY

  !Couple the interface meshes
  PRINT *, ' == >> CREATING INTERFACE MESHES CONNECTIVITY << == '
  CALL CMISSInterfaceMeshConnectivity_Initialise(InterfaceMeshConnectivity,Err)
  CALL CMISSInterfaceMeshConnectivity_CreateStart(Interface,InterfaceMesh,InterfaceMeshConnectivity,Err)
  CALL CMISSInterfaceMeshConnectivity_SetBasis(InterfaceMeshConnectivity,InterfaceMappingBasis,Err)

  DO ic_idx=1,CMX%NUMBER_OF_COUPLINGS
    !Map the interface element to the elements in mesh 1
    CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity,CMX%INTERFACE_ELEMENT_NUMBER(ic_idx), &
      & CMX%MESH1_ID,CMX%MESH1_ELEMENT_NUMBER(ic_idx),Err)
    !Map the interface element to the elements in mesh 2
    CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity,CMX%INTERFACE_ELEMENT_NUMBER(ic_idx), &
      & CMX%MESH2_ID,CMX%MESH2_ELEMENT_NUMBER(ic_idx),Err)
  ENDDO !ic_idx

  DO ic_idx=1,CMX%NUMBER_OF_COUPLINGS
    !Define xi mapping in mesh 1
    CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,CMX%INTERFACE_ELEMENT_NUMBER(ic_idx), & 
      & CMX%MESH1_ID,CMX%MESH1_ELEMENT_NUMBER(ic_idx),CMX%INTERFACE_ELEMENT_LOCAL_NODE(ic_idx),1, &
      & CMX%MESH1_ELEMENT_XI(ic_idx,1:3),Err)
    !Define xi mapping in mesh 2
    CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity,CMX%INTERFACE_ELEMENT_NUMBER(ic_idx), & 
      & CMX%MESH2_ID,CMX%MESH2_ELEMENT_NUMBER(ic_idx),CMX%INTERFACE_ELEMENT_LOCAL_NODE(ic_idx),1, &
      & CMX%MESH2_ELEMENT_XI(ic_idx,1:3),Err)
  ENDDO !ic_idx


  CALL CMISSInterfaceMeshConnectivity_CreateFinish(InterfaceMeshConnectivity,Err)


  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD & DECOMPOSITION

  !Create a decomposition for mesh1
  PRINT *, ' == >> CREATING MESH(1) DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(Decomposition1,Err)
  CALL CMISSDecomposition_CreateStart(Decomposition1UserNumber,Mesh1,Decomposition1,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition1,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition1,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition1,Err)

  !Create a decomposition for mesh2
  PRINT *, ' == >> CREATING MESH(2) DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(Decomposition2,Err)
  CALL CMISSDecomposition_CreateStart(Decomposition2UserNumber,Mesh2,Decomposition2,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition2,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition2,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition2,Err)
  
  !Create a decomposition for the interface mesh
  PRINT *, ' == >> CREATING INTERFACE DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(InterfaceDecomposition,Err)
  CALL CMISSDecomposition_CreateStart(InterfaceDecompositionUserNumber,InterfaceMesh,InterfaceDecomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(InterfaceDecomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(InterfaceDecomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(InterfaceDecomposition,Err)

  !Start to create a default (geometric) field on the first region
  PRINT *, ' == >> CREATING MESH(1) GEOMETRIC FIELD << == '
  CALL CMISSField_Initialise(GeometricField1,Err)
  CALL CMISSField_CreateStart(GeometricField1UserNumber,Region1,GeometricField1,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField1,Decomposition1,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_OF_DIMENSIONS1==3.AND.NUMBER_OF_DIMENSIONS2==3) THEN
    CALL CMISSField_ComponentMeshComponentSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the first field
  CALL CMISSField_CreateFinish(GeometricField1,Err)

  !Start to create a default (geometric) field on the second region
  PRINT *, ' == >> CREATING MESH(2) GEOMETRIC FIELD << == '
  CALL CMISSField_Initialise(GeometricField2,Err)
  CALL CMISSField_CreateStart(GeometricField2UserNumber,Region2,GeometricField2,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField2,Decomposition2,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_OF_DIMENSIONS1==3.AND.NUMBER_OF_DIMENSIONS2==3) THEN
    CALL CMISSField_ComponentMeshComponentSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the second field
  CALL CMISSField_CreateFinish(GeometricField2,Err)

  !Update the geometric field parameters for the first field
  DO NODE_NUMBER=1,NUMBER_OF_NODES_SPACE1
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS1
      VALUE=CM1%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, & 
        & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSField_ParameterSetUpdateStart(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Update the geometric field parameters for the second field
  DO NODE_NUMBER=1,NUMBER_OF_NODES_SPACE2
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS2
      VALUE=CM2%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, & 
        & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSField_ParameterSetUpdateStart(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)


  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

   !Create the equations set for the first region
  PRINT *, ' == >> CREATING EQUATION SET(1) << == '
  CALL CMISSField_Initialise(EquationsSetField1,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSet1,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSet1UserNumber,Region1,GeometricField1,CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & CMISS_EQUATIONS_SET_STOKES_EQUATION_TYPE,CMISS_EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EquationsSetField1UserNumber,&
    & EquationsSetField1,EquationsSet1,Err)
  !Set the equations set to be a standard Stokes problem
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet1,Err)

  !Create the equations set for the second region
  PRINT *, ' == >> CREATING EQUATION SET(2) << == '
  CALL CMISSField_Initialise(EquationsSetField2,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSet2,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSet2UserNumber,Region2,GeometricField2,CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & CMISS_EQUATIONS_SET_STOKES_EQUATION_TYPE,CMISS_EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EquationsSetField2UserNumber,&
    & EquationsSetField2,EquationsSet2,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet2,Err)

  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for the first equations set
  PRINT *, ' == >> CREATING DEPENDENT FIELD(1) << == '
  CALL CMISSField_Initialise(DependentField1,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet1,DependentField1UserNumber,DependentField1,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet1,Err)

  !Create the equations set dependent field variables for the second equations set
  PRINT *, ' == >> CREATING DEPENDENT FIELD(2) << == '
  CALL CMISSField_Initialise(DependentField2,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet2,DependentField2UserNumber,DependentField2,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet2,Err)


  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for Field 1
  CALL CMISSField_Initialise(MaterialsField1,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet1,MaterialsComponentUserNumber1, & 
    & MaterialsField1,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet1,Err)

  !Create the equations set materials field variables for Field 2
  CALL CMISSField_Initialise(MaterialsField2,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet2,MaterialsComponentUserNumber2, & 
    & MaterialsField2,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet2,Err)

  
  !Initialise same material parameters for each field
  CALL CMISSField_ComponentValuesInitialise(MaterialsField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialsComponentUserNumber1,MU_PARAM,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialsComponentUserNumber2,RHO_PARAM,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialsComponentUserNumber1,MU_PARAM,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialsComponentUserNumber2,RHO_PARAM,Err)


  !
  !================================================================================================================================
  !

  !EQUATIONS

  !Create the equations set equations for the first equations set
  PRINT *, ' == >> CREATING EQUATIONS(1) << == '
  CALL CMISSEquations_Initialise(Equations1,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet1,Equations1,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations1,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  !CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet1,Err)

  !Create the equations set equations for the second equations set
  PRINT *, ' == >> CREATING EQUATIONS(2) << == '
  CALL CMISSEquations_Initialise(Equations2,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet2,Equations2,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations2,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  !CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet2,Err)


  !
  !================================================================================================================================
  !


  !INTERFACE GEOMETRIC FIELD

  !Start to create a default (geometric) field on the Interface
  PRINT *, ' == >> CREATING INTERFACE GEOMETRIC FIELD << == '
  CALL CMISSField_Initialise(InterfaceGeometricField,Err)
  CALL CMISSField_CreateStart(InterfaceGeometricFieldUserNumber,Interface,InterfaceGeometricField,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(InterfaceGeometricField,InterfaceDecomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(InterfaceGeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(InterfaceGeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_OF_DIMENSIONS1==3.AND.NUMBER_OF_DIMENSIONS2==3.AND.NUMBER_OF_DIMENSIONS_INTERFACE==2) THEN
    CALL CMISSField_ComponentMeshComponentSet(InterfaceGeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the first field
  CALL CMISSField_CreateFinish(InterfaceGeometricField,Err)

 !Update the geometric field parameters for the interface field
  DO NODE_NUMBER=1,NUMBER_OF_NODES_INTERFACE
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS_INTERFACE+1
      VALUE=CM3%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSField_ParameterSetUpdateNode(InterfaceGeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, & 
        & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSField_ParameterSetUpdateStart(InterfaceGeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(InterfaceGeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Create an interface condition between the two meshes
  PRINT *, ' == >> CREATING INTERFACE CONDITIONS << == '
  CALL CMISSInterfaceCondition_Initialise(InterfaceCondition,Err)
  CALL CMISSInterfaceCondition_CreateStart(InterfaceConditionUserNumber,Interface,InterfaceGeometricField, &
    & InterfaceCondition,Err)
  !Specify the method for the interface condition
  CALL CMISSInterfaceCondition_MethodSet(InterfaceCondition,CMISS_INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,Err)

! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!

  !Specify the type of interface condition operator
  CALL CMISSInterfaceCondition_OperatorSet(InterfaceCondition,CMISS_INTERFACE_CONDITION_FIELD_CONTINUITY_OPERATOR,Err)

! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!

  !Add in the dependent variables from the equations sets
  CALL CMISSInterfaceCondition_DependentVariableAdd(InterfaceCondition,Mesh1Index,EquationsSet1, &
    & CMISS_FIELD_U_VARIABLE_TYPE,Err)
  CALL CMISSInterfaceCondition_DependentVariableAdd(InterfaceCondition,Mesh2Index,EquationsSet2, &
    & CMISS_FIELD_U_VARIABLE_TYPE,Err)
  !Finish creating the interface condition
  CALL CMISSInterfaceCondition_CreateFinish(InterfaceCondition,Err)

  !Create the Lagrange multipliers field
  PRINT *, ' == >> CREATING INTERFACE LAGRANGE FIELD << == '
  CALL CMISSField_Initialise(LagrangeField,Err)
  CALL CMISSInterfaceCondition_LagrangeFieldCreateStart(InterfaceCondition,LagrangeFieldUserNumber,LagrangeField,Err)
  !Finish the Lagrange multipliers field
  CALL CMISSInterfaceCondition_LagrangeFieldCreateFinish(InterfaceCondition,Err)

  !Create the interface condition equations
  PRINT *, ' == >> CREATING INTERFACE EQUATIONS << == '
  CALL CMISSInterfaceEquations_Initialise(InterfaceEquations,Err)
  CALL CMISSInterfaceCondition_EquationsCreateStart(InterfaceCondition,InterfaceEquations,Err)
  !Set the interface equations sparsity
  CALL CMISSInterfaceEquations_SparsitySet(InterfaceEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the interface equations output
  CALL CMISSInterfaceEquations_OutputTypeSet(InterfaceEquations,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !Finish creating the interface equations
  CALL CMISSInterfaceCondition_EquationsCreateFinish(InterfaceCondition,Err)

  !
  !================================================================================================================================
  !

  !PROBLEMS
  
  !Start the creation of a coupled problem.
  PRINT *, ' == >> CREATING PROBLEM << == '
  CALL CMISSProblem_Initialise(CoupledProblem,Err)
  CALL CMISSProblem_CreateStart(CoupledProblemUserNumber,CoupledProblem,Err)
  !Set the problem to be a standard Stokes problem
  CALL CMISSProblem_SpecificationSet(CoupledProblem,CMISS_PROBLEM_FLUID_MECHANICS_CLASS, &
    & CMISS_PROBLEM_STOKES_EQUATION_TYPE,CMISS_PROBLEM_STATIC_STOKES_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(CoupledProblem,Err)

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem control loop for the coupled problem
  PRINT *, ' == >> CREATING PROBLEM CONTROL LOOP << == '
  CALL CMISSProblem_ControlLoopCreateStart(CoupledProblem,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(CoupledProblem,Err)
 
  !Start the creation of the problem solver for the coupled problem
  PRINT *, ' == >> CREATING PROBLEM SOLVERS << == '
  CALL CMISSSolver_Initialise(CoupledSolver,Err)
  CALL CMISSProblem_SolversCreateStart(CoupledProblem,Err)
  CALL CMISSProblem_SolverGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE,1,CoupledSolver,Err)
  !CALL CMISSSolver_OutputTypeSet(CoupledSolver,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(CoupledSolver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(CoupledSolver,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(CoupledSolver,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  CALL CMISSSolver_OutputTypeSet(CoupledSolver,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  CALL CMISSSolver_LinearTypeSet(CoupledSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSSolver_LibraryTypeSet(CoupledSolver,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(CoupledProblem,Err)


  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations for the coupled problem
  PRINT *, ' == >> CREATING PROBLEM SOLVER EQUATIONS << == '
  CALL CMISSSolver_Initialise(CoupledSolver,Err)
  CALL CMISSSolverEquations_Initialise(CoupledSolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(CoupledProblem,Err)
  !Get the solve equations
  CALL CMISSProblem_SolverGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE,1,CoupledSolver,Err)
  CALL CMISSSolver_SolverEquationsGet(CoupledSolver,CoupledSolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquations_SparsityTypeSet(CoupledSolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !CALL CMISSSolverEquations_SparsityTypeSet(CoupledSolverEquations,CMISS_SOLVER_FULL_MATRICES,Err)  
  !Add in the first equations set
  CALL CMISSSolverEquations_EquationsSetAdd(CoupledSolverEquations,EquationsSet1,EquationsSet1Index,Err)
  !Add in the second equations set
  CALL CMISSSolverEquations_EquationsSetAdd(CoupledSolverEquations,EquationsSet2,EquationsSet2Index,Err)

! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!

  !Add in the interface condition
! ! !   CALL CMISSSolverEquations_InterfaceConditionAdd(CoupledSolverEquations,InterfaceCondition,InterfaceConditionIndex,Err)

! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!

  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(CoupledProblem,Err)

  !
  !================================================================================================================================
  !


  !BOUNDARY CONDITIONS

  !Start the creation of the equations set boundary conditions for the first equations set
  PRINT *, ' == >> CREATING BOUNDARY CONDITIONS << == '
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(CoupledSolverEquations,BoundaryConditions,Err)

  IF(FIXED_WALL_NODES1_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_FIXED_WALL_NODES1
      NODE_NUMBER=FIXED_WALL_NODES1(NODE_COUNTER)
      CONDITION=CMISS_BOUNDARY_CONDITION_FIXED_WALL
      CALL CMISSDecomposition_NodeDomainGet(Decomposition1,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS1
          VALUE=0.0_CMISSDP
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_NO_GLOBAL_DERIV,1,&
            & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !Set velocity boundary conditions
  IF(INLET_WALL_NODES1_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES1
      NODE_NUMBER=INLET_WALL_NODES1(NODE_COUNTER)
      CONDITION=CMISS_BOUNDARY_CONDITION_FIXED_INLET
      CALL CMISSDecomposition_NodeDomainGet(Decomposition1,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS1
          VALUE=BOUNDARY_CONDITIONS1(COMPONENT_NUMBER)
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_NO_GLOBAL_DERIV,1,&
            & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  
  !Set boundary conditions for second dependent field
  IF(FIXED_WALL_NODES2_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_FIXED_WALL_NODES2
      NODE_NUMBER=FIXED_WALL_NODES2(NODE_COUNTER)
      CONDITION=CMISS_BOUNDARY_CONDITION_FIXED_WALL
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS2
          VALUE=0.0_CMISSDP
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
            & CMISS_NO_GLOBAL_DERIV,&
            & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !Set velocity boundary conditions
  IF(INLET_WALL_NODES2_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES2
      NODE_NUMBER=INLET_WALL_NODES2(NODE_COUNTER)
      CONDITION=CMISS_BOUNDARY_CONDITION_FIXED_INLET
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS2
          VALUE=BOUNDARY_CONDITIONS2(COMPONENT_NUMBER)
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1, &
            & CMISS_NO_GLOBAL_DERIV,&
            & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(CoupledSolverEquations,Err)

  !
  !================================================================================================================================
  !

  !RUN SOLVERS


! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!

! ! !   !Solve the problem
! ! !   PRINT *, ' == >> SOLVING PROBLEM << == '
! ! !   CALL CMISSProblem_Solve(CoupledProblem,Err)

! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!
! ! ! THIS NEEDS TO BE CHANGES!!!

  !Export the fields
  PRINT *, ' == >> EXPORTING FIELDS << == '
  CALL CMISSFields_Initialise(Fields1,Err)
  CALL CMISSFields_Create(Region1,Fields1,Err)
  CALL CMISSFields_NodesExport(Fields1,"CoupledStokes_1","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields1,"CoupledStokes_1","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields1,Err)
  CALL CMISSFields_Initialise(Fields2,Err)
  CALL CMISSFields_Create(Region2,Fields2,Err)
  CALL CMISSFields_NodesExport(Fields2,"CoupledStokes_2","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields2,"CoupledStokes_2","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields2,Err)
  CALL CMISSFields_Initialise(InterfaceFields,Err)
  CALL CMISSFields_Create(INTERFACE,InterfaceFields,Err)
  CALL CMISSFields_NodesExport(InterfaceFields,"CoupledStokes_Interface","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(InterfaceFields,"CoupledStokes_Interface","FORTRAN",Err)
  CALL CMISSFields_Finalise(InterfaceFields,Err)
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP
 
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
     
END PROGRAM COUPLEDSTOKES
