!> \file
!> $Id: LinearElasticityExample.f90 20 2009-02-15 13:26:52Z cpb $
!> \author Chris Bradley
!> \brief This is an example program to solve a linear elasticity equation using openCMISS calls.
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

!> \example LinearElasticity/src/LinearElasticityExample.f90
!! Example program to solve a linear elasticity equation using openCMISS calls.
!<

!> Main program
PROGRAM 3DLagrangeBasis

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=10

  !Program types
  !Element Node Types for prescribing Nodal positions
  !ELEMENTS(:) 					#TYPE(ELEMENT_TYPE)
  !  MESH_COMP(:) 				#TYPE(ELEMENT_NODES_TYPE)
  !    ELEMENT					#TYPE(MESH_ELEMENTS_TYPE)
  !    ELEM_NODE_NUMBERS(:)		#INTEGER(CMISSIntg),ALLOCATABLE
  !    DERIV(8)					#TYPE(NODE_COORDINATES_TYPE)
  !      NODE_COORDINATES(:)	#REAL(CMISSDP),ALLOCATABLE 
  TYPE NODE_COORDINATES_TYPE
    REAL(CMISSDP),ALLOCATABLE :: NODE_COORDINATES(:)
  END TYPE NODE_COORDINATES_TYPE
  TYPE ELEMENT_NODES_TYPE
    TYPE(CMISSMeshElementsType) :: ELEMENT
    INTEGER(CMISSIntg),ALLOCATABLE :: ELEM_NODE_NUMBERS(:)
    TYPE(NODE_COORDINATES_TYPE) :: DERIV(8)
  END TYPE ELEMENT_NODES_TYPE
  TYPE ELEMENT_TYPE
    TYPE(ELEMENT_NODES_TYPE),POINTER :: MESH_COMP(:)
  END TYPE ELEMENT_TYPE
  TYPE(ELEMENT_TYPE), POINTER :: ELEMENTS(:)

  !Boundary Condition Types
  !DISP_BC 								#TYPE(BC_MESH_COMP)
  !  MESH_COMP(:) 						#TYPE(BC_ELEMENT_NODES_TYPE)
  !    NUMBER_OF_BC_NODES_IN_DERIV(8)	#INTEGER(CMISSIntg)
  !    DERIV(8)							#TYPE(BC_NODE_COORDINATES_TYPE)
  !      NODE_NUMBER(:)					#INTEGER(CMISSIntg),ALLOCATABLE
  !      NODE_COORDINATES(:)			#REAL(CMISSDP),ALLOCATABLE 
  TYPE BC_NODE_COORDINATES_TYPE
    INTEGER(CMISSIntg),ALLOCATABLE :: NODE_NUMBER(:)
    REAL(CMISSDP),ALLOCATABLE :: NODE_COORDINATES(:)
  END TYPE BC_NODE_COORDINATES_TYPE
  TYPE BC_ELEMENT_NODES_TYPE
    INTEGER(CMISSIntg) :: NUMBER_OF_BC_NODES_IN_DERIV(8)
    TYPE(BC_NODE_COORDINATES_TYPE) :: DERIV(8)
  END TYPE BC_ELEMENT_NODES_TYPE
  TYPE BC_MESH_COMP
    TYPE(BC_ELEMENT_NODES_TYPE),POINTER :: MESH_COMP(:)
  END TYPE BC_MESH_COMP
  TYPE(BC_MESH_COMP), POINTER :: DISP_BC,FORCE_BC

  !Program variables
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NUMBER_COMPUTATIONAL_NODES,MY_COMPUTATIONAL_NODE_NUMBER
  INTEGER(CMISSIntg) :: CS_USER_NUMBER, REGION_USER_NUMBER, BASIS_USER_NUMBER, MESH_USER_NUMBER, DECOMPOSITION_USER_NUMBER
  INTEGER(CMISSIntg) :: FIELD_USER_NUMBER,EQUATION_SET_USER_NUMBER, PROBLEM_USER_NUMBER
  INTEGER(CMISSIntg) :: NUMBER_OF_MESH_COMPS,NUMBER_OF_XI,NUMBER_OF_NODES,NUMBER_OF_ELEMENTS,NUMBER_OF_DERIV
  INTEGER(CMISSIntg) :: NUMBER_OF_FIELD_VARIABLES,NUMBER_OF_BC_NODES,NUMBER_OF_FIELD_COMPS
  INTEGER(CMISSIntg) :: MESH_COMP_IDX,FIELD_COMP_IDX,FIELD_DERIVATIVE_IDX,SOLVER_IDX,EQUATION_SET_IDX,np,xi,ne,nu
  INTEGER(CMISSIntg) :: NUMBER_OF_GLOBAL_DEPENDENT_DOFS
  INTEGER(CMISSIntg),ALLOCATABLE :: XI_INTERPOLATION(:,:),NUMBER_OF_GAUSS_POINTS(:),MESH_COMP_NUMBER_OF_ELEMENT_NODES(:)

  REAL(CMISSDP) :: l,w,h
  REAL(CMISSDP), POINTER :: FIELD_DATA(:)
  REAL(CMISSDP),ALLOCATABLE :: MATE_PARA(:)

  LOGICAL :: EXPORT_FIELD

  CHARACTER(LEN=255) :: ERROR

  !CMISS variables
  TYPE(CMISSBasisType),ALLOCATABLE :: Basis_xyz(:)
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

  !Generic CMISS variables 
  INTEGER(CMISSIntg) :: Err

  !WIN32 Variables
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

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


  !=================================================================================================================================

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !=BROADCAST PARAMETERS TO COMPUTATIONAL NODES====================================================================================
  NUMBER_GLOBAL_X_ELEMENTS=2
  NUMBER_GLOBAL_Y_ELEMENTS=2
  NUMBER_GLOBAL_Z_ELEMENTS=0
  NUMBER_OF_DOMAINS=1
  NUMBER_OF_XI = 3
  NUMBER_OF_ELEMENTS = 1
  NUMBER_OF_NODES = 8
    
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !=CREATE COORDINATE SYSTEM=======================================================================================================
  !Start the creation of a new 3D RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemTypeSet(CoordinateSystem,CMISSCOORDINATERECTANGULARCARTESIANTYPE,Err)
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NUMBER_OF_XI,Err)
  CALL CMISSCoordinateSystemOriginSet(CoordinateSystem,(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/),Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !=CREATE REGION==================================================================================================================
  !Create Region and set CS to newly created 3D RC CS
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !=CREATE BASIS ==================================================================================================================
  ALLOCATE(Basis_xyz(NUMBER_OF_XI),STAT=Err)
  IF(Err/=0) CALL FLAG_ERROR("Could not allocate Basis_xyz",Err,ERROR,*999)
  ALLOCATE(XI_INTERPOLATION(NUMBER_OF_XI,NUMBER_OF_XI),STAT=Err)
  IF(Err/=0) CALL FLAG_ERROR("Could not allocate XI_INTERPOLATION",Err,ERROR,*999)
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe Interpolation in each direction
  !NOTE if you change interpolation you need to change Boundary Conditions
  XI_INTERPOLATION(1,:) = &
    & (/CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation/)
  XI_INTERPOLATION(2,:) = &
    & (/CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation/)
  XI_INTERPOLATION(3,:) = &
    & (/CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation/)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ALLOCATE(NUMBER_OF_GAUSS_POINTS(NUMBER_OF_XI),STAT=Err)
  IF(Err/=0) CALL FLAG_ERROR("Could not allocate NUMBER_OF_GAUSS_POINTS",Err,ERROR,*999)
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe Number of Gauss Points
  !NOTE:: Num of Gauss points must be the same across X,Y & Z coordinates and be sufficient to accurately integrate the hightest order interpolation being used
  NUMBER_OF_GAUSS_POINTS(:) = (/4,4,4/)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  DO xi=1,NUMBER_OF_XI
    BASIS_USER_NUMBER = xi
    CALL CMISSBasisTypeInitialise(Basis_xyz(xi),Err)
    CALL CMISSBasisCreateStart(BASIS_USER_NUMBER,Basis_xyz(xi),Err)
    CALL CMISSBasisTypeSet(Basis_xyz(xi),CMISSBasisLagrangeHermiteTPType,Err)
    CALL CMISSBasisNumberOfXiSet(Basis_xyz(xi),NUMBER_OF_XI,Err)
    CALL CMISSBasisInterpolationXiSet(Basis_xyz(xi),XI_INTERPOLATION(xi,:),Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis_xyz(xi),NUMBER_OF_GAUSS_POINTS,Err)
    CALL CMISSBasisCreateFinish(Basis_xyz(xi),Err)
  ENDDO !xi

  !=MANUAL MESH CREATION===========================================================================================================

  !=NODE CREATION==================================================================================================================
  !Create nodes in REGION and set initial coordinates to 0,0,0
  CALL CMISSNodesTypeInitialise(NODES,Err)
  CALL CMISSNodesCreateStart(REGION,NUMBER_OF_NODES,NODES,Err)
  CALL CMISSNodesCreateFinish(NODES,Err)
  
  !=CREATE MESH====================================================================================================================
  !Create a mesh with xi number of coordinates
  MESH_USER_NUMBER = 1
  NUMBER_OF_MESH_COMPS = NUMBER_OF_XI
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSMeshCreateStart(MESH_USER_NUMBER,REGION,NUMBER_OF_XI,MESH,Err)
  CALL CMISSMeshNumberOfElementsSet(MESH,NUMBER_OF_ELEMENTS,Err)
  CALL CMISSMeshNumberOfComponentsSet(MESH,NUMBER_OF_MESH_COMPS,Err)

  ALLOCATE(ELEMENTS(NUMBER_OF_ELEMENTS),STAT=Err)
  IF(Err/=0) CALL FLAG_ERROR("Could not allocate MATE_PARA",Err,ERROR,*999)
  DO ne=1,NUMBER_OF_ELEMENTS
    ALLOCATE(ELEMENTS(ne)%MESH_COMP(NUMBER_OF_MESH_COMPS),STAT=Err)
    IF(Err/=0) CALL FLAG_ERROR("Could not allocate ELEMENT(ne)%MESH_COMP",Err,ERROR,*999)
  ENDDO !ne
  ALLOCATE(MESH_COMP_NUMBER_OF_ELEMENT_NODES(NUMBER_OF_MESH_COMPS),STAT=Err)
  IF(Err/=0) CALL FLAG_ERROR("Could not allocate MESH_COMP_NUMBER_OF_ELEMENT_NODES",Err,ERROR,*999)
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX=xi
    IF (XI_INTERPOLATION(1,1)<=3) THEN
      MESH_COMP_NUMBER_OF_ELEMENT_NODES(MESH_COMP_IDX)=PRODUCT(XI_INTERPOLATION(xi,:)+1)
      NUMBER_OF_DERIV = 1
    ELSEIF (XI_INTERPOLATION(1,1)==4) THEN
      MESH_COMP_NUMBER_OF_ELEMENT_NODES(MESH_COMP_IDX)=(2**NUMBER_OF_XI)
      NUMBER_OF_DERIV = 2**NUMBER_OF_XI
    ELSE
      CALL FLAG_ERROR("This example is only setup for Lagrange and cubic hermite bases",Err,ERROR,*999)
    ENDIF
  ENDDO !xi
  DO ne=1,NUMBER_OF_ELEMENTS
    DO xi=1,NUMBER_OF_XI
      MESH_COMP_IDX=xi
      ALLOCATE(ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEM_NODE_NUMBERS(MESH_COMP_NUMBER_OF_ELEMENT_NODES(MESH_COMP_IDX)),STAT=Err)
      IF(Err/=0) CALL FLAG_ERROR("Could not allocate ELEMENT(ne)%MESH_COMP(MESH_COMP_IDX)%ELEM_NODE_NUMBERS",Err,ERROR,*999)
      ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEM_NODE_NUMBERS = 0
    ENDDO !xi
  ENDDO !ne
!  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  !Prescribe Element node numbers
  ELEMENTS(1)%MESH_COMP(1)%ELEM_NODE_NUMBERS(:) = (/1,2,3,4,5,6,7,8/)
  ELEMENTS(1)%MESH_COMP(2)%ELEM_NODE_NUMBERS(:) = (/1,2,3,4,5,6,7,8/)
  ELEMENTS(1)%MESH_COMP(3)%ELEM_NODE_NUMBERS(:) = (/1,2,3,4,5,6,7,8/)
!  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  !!TODO:: When diagnostics are on - MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET trys to output MESH%TOPOLOGY(1)%PTR%ELEMENTS%ELEMENTS(1)%MESH_ELEMENT_NODES which are only allocated when the MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH command is given
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX=xi
    DO ne=1,NUMBER_OF_ELEMENTS
      CALL CMISSMeshElementsCreateStart(MESH,MESH_COMP_IDX,Basis_xyz(xi),ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEMENT,Err)
      CALL CMISSMeshElementsNodesSet(ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEMENT,ne, & 
        & ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEM_NODE_NUMBERS(:),Err)
      CALL CMISSMeshElementsCreateFinish(ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEMENT,Err)
    ENDDO !ne
  ENDDO !xi
  CALL CMISSMeshCreateFinish(MESH,Err)

  !=CREATE DECOMPOSITION===========================================================================================================
  !Create mesh decomposition dividing mesh into number_of_domains for parallel solving
  DECOMPOSITION_USER_NUMBER = 1
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DECOMPOSITION_USER_NUMBER,MESH,DECOMPOSITION,Err)
  CALL CMISSDecompositionTypeSet(DECOMPOSITION,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(DECOMPOSITION,NUMBER_OF_DOMAINS,Err)
  CALL CMISSDecompositionCreateFinish(DECOMPOSITION,Err)

  !=CREATE GEOMETRIC FIELD=========================================================================================================
  !Start to create a default (geometric) field on the region
  FIELD_USER_NUMBER = 1
  NUMBER_OF_FIELD_VARIABLES = 1 !Geometric Field Coordinates
  NUMBER_OF_FIELD_COMPS = NUMBER_OF_XI
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(FIELD_USER_NUMBER,REGION,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,DECOMPOSITION,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,NUMBER_OF_FIELD_VARIABLES,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,NUMBER_OF_FIELD_COMPS,Err)
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX = xi
    FIELD_COMP_IDX = xi
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,FIELD_COMP_IDX,MESH_COMP_IDX,Err)
  ENDDO !xi
  CALL CMISSFieldCreateFinish(GeometricField,Err)
  !Initialize node coordinates
  DO ne=1,NUMBER_OF_ELEMENTS
    DO xi=1,NUMBER_OF_XI
      MESH_COMP_IDX=xi
      DO nu=1,NUMBER_OF_DERIV
        ALLOCATE(ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES&
          &(MESH_COMP_NUMBER_OF_ELEMENT_NODES(MESH_COMP_IDX)),STAT=Err)
        IF(Err/=0) CALL FLAG_ERROR("Could not allocate ELEMENT(ne)%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES", & 
          & Err,ERROR,*999)
        ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES = 0.0_CMISSDP
      ENDDO
    ENDDO
  ENDDO
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe node coordinates
  l = 120.0_CMISSDP
  w = 160.0_CMISSDP
  h = 10.0_CMISSDP
  ELEMENTS(1)%MESH_COMP(1)%DERIV(1)%NODE_COORDINATES = (/0.0_CMISSDP,l,0.0_CMISSDP,l,0.0_CMISSDP,l,0.0_CMISSDP,l/)
  ELEMENTS(1)%MESH_COMP(2)%DERIV(1)%NODE_COORDINATES = (/0.0_CMISSDP,0.0_CMISSDP,w,w,0.0_CMISSDP,0.0_CMISSDP,w,w/)
  ELEMENTS(1)%MESH_COMP(3)%DERIV(1)%NODE_COORDINATES = (/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,h,h,h,h/)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX = xi
    FIELD_COMP_IDX = xi
    DO ne=1,NUMBER_OF_ELEMENTS
      DO nu=1,NUMBER_OF_DERIV
        DO np=1,MESH_COMP_NUMBER_OF_ELEMENT_NODES(MESH_COMP_IDX) !Can also loop using nu
           CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,nu, &
            & ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEM_NODE_NUMBERS(np),FIELD_COMP_IDX, &
            & ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES(np),Err)
        ENDDO !np
      ENDDO !nu
    ENDDO !ne
  ENDDO !xi

  !=CREATE DEPENDENT FIELD=========================================================================================================
  !Create a dependent field
  FIELD_USER_NUMBER = 2
  NUMBER_OF_FIELD_VARIABLES = 2 !Dependent Field Displacement Coordinates & Force 
  NUMBER_OF_FIELD_COMPS = NUMBER_OF_XI
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSFieldCreateStart(FIELD_USER_NUMBER,REGION,DependentField,Err)
  CALL CMISSFieldTypeSet(DependentField,CMISSFieldGeneralType,Err)
  CALL CMISSFieldMeshDecompositionSet(DependentField,DECOMPOSITION,Err)
  CALL CMISSFieldGeometricFieldSet(DependentField,GeometricField,Err)
  CALL CMISSFieldDependentTypeSet(DependentField,CMISSFieldDependentType,Err)
  CALL CMISSFieldNumberOfVariablesSet(DependentField,NUMBER_OF_FIELD_VARIABLES,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldUVariableType,NUMBER_OF_FIELD_COMPS,Err)
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX = xi
    FIELD_COMP_IDX = xi
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,FIELD_COMP_IDX,MESH_COMP_IDX,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,FIELD_COMP_IDX,MESH_COMP_IDX,Err)
  ENDDO !xi
  CALL CMISSFieldCreateFinish(DependentField,Err)

  !=CREATE MATERIAL FIELD==========================================================================================================
  !!TODO:: Set Material Field Interpolation to constant element based interpolation when field i/o and cmgui allows for this
  !Create a material field for a general 2D isotropic material
  FIELD_USER_NUMBER = 3
  NUMBER_OF_FIELD_VARIABLES = 1
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe material properties E1,E2,E3 & v13,v23,v12
  NUMBER_OF_FIELD_COMPS = 6 !Young's Modulus & Poisson's Ratio
  ALLOCATE(MATE_PARA(NUMBER_OF_FIELD_COMPS),STAT=Err)
  IF(Err/=0) CALL FLAG_ERROR("Could not allocate MATE_PARA",Err,ERROR,*999)
  MATE_PARA = (/30E6_CMISSDP,30E6_CMISSDP,30E6_CMISSDP,0.25_CMISSDP,0.25_CMISSDP,0.25_CMISSDP/)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  CALL CMISSFieldTypeInitialise(MaterialField,Err)
  CALL CMISSFieldCreateStart(FIELD_USER_NUMBER,REGION,MaterialField,Err)
  CALL CMISSFieldTypeSet(MaterialField,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialField,DECOMPOSITION,Err)
  CALL CMISSFieldGeometricFieldSet(MaterialField,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialField,NUMBER_OF_FIELD_VARIABLES,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialField,CMISSFieldUVariableType,NUMBER_OF_FIELD_COMPS,Err)
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX = xi
    DO FIELD_COMP_IDX=1,NUMBER_OF_FIELD_COMPS
      CALL CMISSFieldComponentMeshComponentSet(MaterialField,CMISSFieldUVariableType,FIELD_COMP_IDX,MESH_COMP_IDX,Err)
    ENDDO !FIELD_COMP_IDX
  ENDDO !xi
  CALL CMISSFieldCreateFinish(MaterialField,Err)
  FIELD_DERIVATIVE_IDX = 1
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX = xi
    DO ne=1,NUMBER_OF_ELEMENTS
      DO np=1,MESH_COMP_NUMBER_OF_ELEMENT_NODES(MESH_COMP_IDX) !Can also loop using nu
        DO FIELD_COMP_IDX=1,NUMBER_OF_FIELD_COMPS
          CALL CMISSFieldParameterSetUpdateNode(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
            & FIELD_DERIVATIVE_IDX,ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEM_NODE_NUMBERS(np),FIELD_COMP_IDX, &
            & MATE_PARA(FIELD_COMP_IDX),Err)
        ENDDO !FIELD_COMP_IDX
      ENDDO !np
    ENDDO !ne
  ENDDO !xi

  !=CREATE EQUATIONS SET===========================================================================================================
  !Create a Elasticity Class, Linear Elasticity type, no subtype, EquationsSet
  EQUATION_SET_USER_NUMBER = 1
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSEquationsSetCreateStart(EQUATION_SET_USER_NUMBER,REGION,GeometricField,EquationsSet,Err)
  !Set the equations set to be a Elasticity Class, Linear Elasticity type, no subtype, EquationsSet
  CALL CMISSEquationsSetSpecificationSet(EquationsSet,CMISSEquationsSetElasticityClass,CMISSEquationsSetLinearElasticityType, &
    & CMISSEquationsSetThreeDimensionalSubtype,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)
  !Create the equations set dependent field variables
  FIELD_USER_NUMBER = 2 !Dependent Field
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,FIELD_USER_NUMBER,DependentField,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)
  !Create the equations set material field variables
  FIELD_USER_NUMBER = 3 !Material Field
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,FIELD_USER_NUMBER,MaterialField,Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  !=CREATE EQUATIONS SET EQUATION==================================================================================================
  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquationsSparsityTypeSet(EQUATIONS,CMISSEquationsFullMatrices,Err)
                                              !CMISSEquationsSparseMatrices=1 !<Use sparse matrices for the equations.
                                              !CMISSEquationsFullMatrices=2 !<Use fully populated matrices for the equations. 
  CALL CMISSEquationsOutputTypeSet(EQUATIONS,CMISSEquationsTimingOutput,Err)
                                            !CMISSEquationsNoOutput !<No output from the equations.
                                            !CMISSEquationsTimingOutput !<Timing information output.
                                            !CMISSEquationsMatrixOutput !<All below and equation matrices output.
                                            !CMISSEquationsElementMatrixOutput !<All below and element matrices output.
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)

  !=PRESCRIBE BOUNDARY CONDITIONS==================================================================================================
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSet,BoundaryConditions,Err)
  !Allocate Memory for Displacement & Force BC
  NULLIFY(DISP_BC,FORCE_BC)
  ALLOCATE(DISP_BC,STAT=Err)
  IF(Err/=0) CALL FLAG_ERROR("Could not allocate DISP_BC",Err,ERROR,*999)
  ALLOCATE(FORCE_BC,STAT=Err)
  IF(Err/=0) CALL FLAG_ERROR("Could not allocate FORCE_BC",Err,ERROR,*999)
  ALLOCATE(DISP_BC%MESH_COMP(MESH_COMP_IDX),STAT=Err)
  IF(Err/=0) CALL FLAG_ERROR("Could not allocate DISP_BC%MESH_COMP",Err,ERROR,*999)
  ALLOCATE(FORCE_BC%MESH_COMP(MESH_COMP_IDX),STAT=Err)
  IF(Err/=0) CALL FLAG_ERROR("Could not allocate FORCE_BC%MESH_COMP",Err,ERROR,*999)
  !Initialize Number of Displacement & Force BC
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX=xi
    DISP_BC%MESH_COMP(MESH_COMP_IDX)%NUMBER_OF_BC_NODES_IN_DERIV=0
    FORCE_BC%MESH_COMP(MESH_COMP_IDX)%NUMBER_OF_BC_NODES_IN_DERIV=0
  ENDDO
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe Number of BC & thier DERIV
  DISP_BC%MESH_COMP(1)%NUMBER_OF_BC_NODES_IN_DERIV(1)=4
  DISP_BC%MESH_COMP(2)%NUMBER_OF_BC_NODES_IN_DERIV=DISP_BC%MESH_COMP(1)%NUMBER_OF_BC_NODES_IN_DERIV
  DISP_BC%MESH_COMP(3)%NUMBER_OF_BC_NODES_IN_DERIV=DISP_BC%MESH_COMP(1)%NUMBER_OF_BC_NODES_IN_DERIV
  FORCE_BC%MESH_COMP(1)%NUMBER_OF_BC_NODES_IN_DERIV(1)=4
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX=xi
    DO nu=1,NUMBER_OF_DERIV
      !Displacement BC
      NUMBER_OF_BC_NODES = DISP_BC%MESH_COMP(MESH_COMP_IDX)%NUMBER_OF_BC_NODES_IN_DERIV(nu)
      IF(NUMBER_OF_BC_NODES/=0) THEN
        ALLOCATE(DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER(NUMBER_OF_BC_NODES),STAT=Err)
        IF(Err/=0) CALL FLAG_ERROR("Could not allocate DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER",Err,ERROR,*999)
        ALLOCATE(DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES(NUMBER_OF_BC_NODES),STAT=Err)
        IF(Err/=0) CALL FLAG_ERROR("Could not allocate DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES",Err,ERROR,*999)
        !Initialize Displacement BC Derivative_Node_Numbers & Node Coordinates
        DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER=0
        DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES=0.0_CMISSDP
      ENDIF
      !Force BC
      NUMBER_OF_BC_NODES = FORCE_BC%MESH_COMP(MESH_COMP_IDX)%NUMBER_OF_BC_NODES_IN_DERIV(nu)
      IF(NUMBER_OF_BC_NODES/=0) THEN
        ALLOCATE(FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER(NUMBER_OF_BC_NODES),STAT=Err)
        IF(Err/=0) CALL FLAG_ERROR("Could not allocate FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER",Err,ERROR,*999)
        ALLOCATE(FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES(NUMBER_OF_BC_NODES),STAT=Err)
        IF(Err/=0) CALL FLAG_ERROR("Could not allocate FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES", &
          & Err,ERROR,*999)
        !Initialize Displacement BC Derivative_Node_Numbers & Node Coordinates
        FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER=0
        FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES=0.0_CMISSDP
      ENDIF
    ENDDO !nu
  ENDDO !xi
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe Displacement BC Derivative Node Numbers & Node Coordinates
  DISP_BC%MESH_COMP(1)%DERIV(1)%NODE_NUMBER(:)=(/1,3,5,7/)
  DISP_BC%MESH_COMP(2)=DISP_BC%MESH_COMP(1)
  DISP_BC%MESH_COMP(3)=DISP_BC%MESH_COMP(1)
  !Prescribe Force BC Derivative Node Numbers & Node Coordinates
  FORCE_BC%MESH_COMP(1)%DERIV(1)%NODE_NUMBER(:)=(/2,4,6,8/)
  FORCE_BC%MESH_COMP(1)%DERIV(1)%NODE_COORDINATES(:)=400.0_CMISSDP
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX=xi
    FIELD_COMP_IDX=xi
    DO nu=1,NUMBER_OF_DERIV
      !Displacement BC
      NUMBER_OF_BC_NODES = DISP_BC%MESH_COMP(MESH_COMP_IDX)%NUMBER_OF_BC_NODES_IN_DERIV(nu)
      IF(NUMBER_OF_BC_NODES/=0) THEN
        DO np=1,NUMBER_OF_BC_NODES
          CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,nu, &
            & DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER(np),FIELD_COMP_IDX,CMISSBoundaryConditionFixed, &
            & DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES(np),Err)
        ENDDO !np
      ENDIF
      !Force BC
      NUMBER_OF_BC_NODES = FORCE_BC%MESH_COMP(MESH_COMP_IDX)%NUMBER_OF_BC_NODES_IN_DERIV(nu)
      IF(NUMBER_OF_BC_NODES/=0) THEN
        DO np=1,NUMBER_OF_BC_NODES
          CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldDelUDelNVariableType,nu, &
            & FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER(np),FIELD_COMP_IDX,CMISSBoundaryConditionFixed, &
            & FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES(np),Err)
        ENDDO !np
      ENDIF
    ENDDO !nu
  ENDDO !xi
  CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSet,Err)
  
  !=CREATE PROBLEM=================================================================================================================
  !Create the problem
  PROBLEM_USER_NUMBER=1
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(PROBLEM_USER_NUMBER,PROBLEM,Err)
  !Set the problem to be a elasticity class, linear elasticity type with no subtype.
  CALL CMISSProblemSpecificationSet(PROBLEM,CMISSProblemElasticityClass,CMISSProblemLinearElasticityType, &
    & CMISSProblemNoSubtype,Err)
  CALL CMISSProblemCreateFinish(PROBLEM,Err)

  !=CREATE PROBLEM CONTROL LOOP====================================================================================================
  !Create the problem control loop
   CALL CMISSProblemControlLoopCreateStart(PROBLEM,Err)
   CALL CMISSProblemControlLoopCreateFinish(PROBLEM,Err)

  !=CREATE PROBLEM SOLVER==========================================================================================================
  !Start the creation of the problem solvers
  SOLVER_IDX = 1
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolversCreateStart(PROBLEM,Err)
  CALL CMISSProblemSolverGet(PROBLEM,CMISSControlLoopNode,SOLVER_IDX,SOLVER,Err)
  CALL CMISSSolverLibraryTypeSet(SOLVER,CMISSSolverPETScLibrary,Err)
                                       !CMISSSolverCMISSLibrary     !<CMISS (internal) solver library.
                                       !CMISSSolverPETScLibrary     !<PETSc solver library.
                                       !CMISSSolverMUMPSLibrary     !<MUMPS solver library.
                                       !CMISSSolverSuperLULibrary   !<SuperLU solver library.
                                       !CMISSSolverSpoolesLULibrary !<SPOOLES solver library.
                                       !CMISSSolverUMFPACKLibrary   !<UMFPACK solver library.
                                       !CMISSSolverLUSOLLibrary     !<LUSOL solver library.
                                       !CMISSSolverESSLLibrary      !<ESSL solver library.
                                       !CMISSSolverLAPACKLibrary    !<LAPACK solver library.
                                       !CMISSSolverTAOLibrary       !<TAO solver library.
                                       !CMISSSolverHypreLibrary     !<Hypre solver library.
  CALL CMISSSolverLinearTypeSet(SOLVER,CMISSSolverLinearDirectSolveType,Err)
                                      !CMISSSolverLinearDirectSolveType    !<Direct linear solver type.
                                      !CMISSSolverLinearIterativeSolveType !<Iterative linear solver type.
  !CALL CMISSSolverLinearDirectTypeSet(SOLVER,CMISSSolverDirectLU,Err)
                                            !CMISSSolverDirectLU       !<LU direct linear solver.
                                            !CMISSSolverDirectCholesky !<Cholesky direct linear solver.
                                            !CMISSSolverDirectSVD      !<SVD direct linear solver.
  CALL CMISSSolverOutputTypeSet(SOLVER,CMISSSolverSolverMatrixOutput,Err)
                                      !CMISSSolverNoOutput !<No output from the solver routines. \see OPENCMISS_SolverOutputTypes,OPENCMISS
                                      !CMISSSolverProgressOutput !<Progress output from solver routines.
                                      !CMISSSolverTimingOutput !<Timing output from the solver routines plus below.
                                      !CMISSSolverSolverOutput !<Solver specific output from the solver routines plus below.
                                      !CMISSSolverSolverMatrixOutput !<Solver matrices output from the solver routines plus below.
  CALL CMISSProblemSolversCreateFinish(PROBLEM,Err)

  !=CREATE PROBLEM SOLVER EQUATIONS================================================================================================
  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(PROBLEM,Err)
  CALL CMISSProblemSolverGet(PROBLEM,CMISSControlLoopNode,SOLVER_IDX,SOLVER,Err)
  CALL CMISSSolverSolverEquationsGet(SOLVER,SolverEquations,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
                                                          !CMISSSolverEquationsSparseMatrices !<Use sparse solver matrices.
                                                          !CMISSSolverEquationsFullMatrices !<Use fully populated solver matrices.
  EQUATION_SET_IDX = 1 !Initialize index of the equations set that has been added 
                        !(Variable is returned from PROBLEM_SolverEquations_EquationsSet_ADD)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EQUATION_SET_IDX,Err)
  CALL CMISSProblemSolverEquationsCreateFinish(PROBLEM,Err)

  !=SOLVE PROBLEM==================================================================================================================
  !Solve the problem
  CALL CMISSProblemSolve(PROBLEM,Err)

  !=OUTPUT SOLUTION================================================================================================================
  !!TODO:: Output reaction forces in ipnode files
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"LinearElasticity3DLagrangeBasisExample","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"LinearElasticity3DLagrangeBasisExample","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
  ENDIF

  !================================================================================================================================

  !Finialise OpenCMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
999 WRITE(*,'(A)') ERROR
  STOP 1

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets the Error string specified by a character string and flags an Error 
  SUBROUTINE FLAG_ERROR(STRING,Err,ERROR,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The Error condition string
    INTEGER(CMISSIntg), INTENT(OUT) :: Err !<The Error code
    CHARACTER(LEN=255), INTENT(OUT) :: ERROR !<The Error string
    !Local variables
    INTEGER(CMISSIntg) :: STRING_LENGTH

    IF(Err==0) Err=1
    STRING_LENGTH=LEN_TRIM(STRING)
    ERROR=STRING(1:STRING_LENGTH)

    RETURN 1
  END SUBROUTINE FLAG_ERROR

END PROGRAM 3DLagrangeBasis


