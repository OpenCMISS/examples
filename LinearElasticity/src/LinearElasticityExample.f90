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
PROGRAM LINEARELASTICITYEXAMPLE

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CMISS
  USE CMISS_MPI
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE EQUATIONS_SET_ROUTINES
  USE FIELD_ROUTINES
  USE FIELD_IO_ROUTINES
  USE LINEAR_ELASTICITY_ROUTINES  
  USE GENERATED_MESH_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MESH_ROUTINES
  USE MPI
  USE NODE_ROUTINES  
  USE PROBLEM_CONSTANTS
  USE PROBLEM_ROUTINES
  USE REGION_ROUTINES
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters
!#### Index: ne
!###  Description:
!###    Index label for a element.
!#### Index: ng
!###  Description:
!###    Index label for a gauss point.
!#### Index: ni
!###  Description:
!###    Index label for a xi direction.
!#### Index: nk
!###  Description:
!###    Index label for a derivative with respect to the global directions.
!#### Index: nn
!###  Description:
!###    Index for a local node within an element.
!#### Index: np
!###  Description:
!###    Index for a node.
!#### Index: ns
!###  Description:
!###    Index for a element parameter within an element.
!#### Index: nu
!###  Description:
!###    Index for a partial derivative.

  !Program types
  TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
  TYPE(REGION_TYPE), POINTER :: REGION,WORLD_REGION
  TYPE(BASIS_TYPE), POINTER :: BASIS
  TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
  TYPE(NODES_TYPE), POINTER :: NODES
  TYPE(MESH_TYPE), POINTER :: MESH
  TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
  TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
  TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD,DEPENDENT_FIELD,MATERIAL_FIELD
  TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
  TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
  TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
  TYPE(SOLVER_TYPE), POINTER :: SOLVER
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
  TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
  TYPE(VARYING_STRING) :: FILE,METHOD

  !> A type to hold the boundary conditions
  TYPE BC_SET_TYPE
    INTEGER(INTG) :: NUM_BC_NODES(3) !<The number of Nodal Components with Boundary Conditions applied for the field
    INTEGER(INTG), ALLOCATABLE :: BC_INFO(:,:) !<BC_NODE(np,i=1..3). The np'th node and i'th component number for the BC
    REAL(DP), ALLOCATABLE :: BC_VALUE(:) !<BC_VALUE(np,i=1..3). The np'th node and i'th component number for the BC Value
  END TYPE BC_SET_TYPE
  TYPE(BC_SET_TYPE) :: DISP_BC,FORCE_BC

  !Program variables
  INTEGER(INTG) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(INTG) :: NUMBER_OF_DOMAINS
  INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES
  INTEGER(INTG) :: MY_COMPUTATIONAL_NODE_NUMBER
  INTEGER(INTG) :: MPI_IERROR
  INTEGER(INTG) :: User_Num_CS, User_Num_REG, User_Num_BASIS, User_Num_MESH, User_Num_DECOMP
  INTEGER(INTG) :: GEOMETRIC_FIELD_USER_NUM,DEPENDENT_FIELD_USER_NUM,MATERIAL_FIELD_USER_NUM
  INTEGER(INTG) :: User_Num_EQUATION_SET, User_Num_PROBLEM
  INTEGER(INTG) :: Num_Coordinates, Num_xi_coordinates, Num_Nodes, Num_Mesh_dim, Num_mesh_components, Num_Elem
  INTEGER(INTG) :: Mesh_component, Num_Field_Variables, Num_Field_Components,Field_Variable
  INTEGER(INTG) :: Field_Derivative, Global_Element
  INTEGER(INTG) :: Solver_Idx,Equations_Set_Idx,BC_idx,Variable_idx,Component_idx,np
  REAL(DP) :: l,w,h,E(3),v(3),E_v(6)
  REAL(DP),ALLOCATABLE :: Node_Coordinates(:,:)
  LOGICAL :: Export_Field,USE_GENERATED_MESHES

  !Timing Variables
  REAL(SP) :: START_USER_TIME(1),STOP_USER_TIME(1),START_SYSTEM_TIME(1),STOP_SYSTEM_TIME(1)

  !Generic CMISS variables
  INTEGER(INTG) :: ERR
  TYPE(VARYING_STRING) :: ERROR

  INTEGER(INTG) :: DIAG_LEVEL_LIST(5)
  CHARACTER(LEN=MAXSTRLEN) :: DIAG_ROUTINE_LIST(3),TIMING_ROUTINE_LIST(1)

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

  !================================================================================================================================

  !Intialise cmiss
  NULLIFY(WORLD_REGION)
  CALL CMISS_INITIALISE(WORLD_REGION,ERR,ERROR,*999)

  !Set all diganostic levels on for testing
  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5
  !DIAG_ROUTINE_LIST(1)="MATRIX_VALUES_ADD_DP2"
  !DIAG_ROUTINE_LIST(2)="EQUATIONS_MATRIX_STRUCTURE_CALCULATE"
  !DIAG_ROUTINE_LIST(3)="EQUATIONS_MAPPING_CALCULATE"

  !CALL DIAGNOSTICS_SET_ON(ALL_DIAG_TYPE,DIAG_LEVEL_LIST,"LinearElasticityExample",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
  !CALL DIAGNOSTICS_SET_ON(IN_DIAG_TYPE,DIAG_LEVEL_LIST,"",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
  !CALL DIAGNOSTICS_SET_ON(IN_DIAG_TYPE,DIAG_LEVEL_LIST,"",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

  !Calculate the start times
  CALL CPU_TIMER(USER_CPU,START_USER_TIME,ERR,ERROR,*999)
  CALL CPU_TIMER(SYSTEM_CPU,START_SYSTEM_TIME,ERR,ERROR,*999)

  !=BROADCAST PARAMETERS TO COMPUTATIONAL NODES====================================================================================
  !Get the number of computational nodes
  NUMBER_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
  IF(ERR/=0) GOTO 999
  !Get my computational node number
  MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
  IF(ERR/=0) GOTO 999
  Number_global_x_elements = 1
  Number_global_y_elements = 1
  Number_global_z_elements = 1
  Number_of_domains = 1
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(Number_global_x_elements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(Number_global_y_elements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(Number_global_z_elements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(Number_of_domains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)

  !=CREATE COORDINATE SYSTEM=======================================================================================================
  !Start the creation of a new 3D RC coordinate system
  User_Num_CS = 1
  Num_Coordinates = 3
  NULLIFY(COORDINATE_SYSTEM)
  CALL COORDINATE_SYSTEM_CREATE_START(User_Num_CS,COORDINATE_SYSTEM,ERR,ERROR,*999)
  CALL COORDINATE_SYSTEM_TYPE_SET(COORDINATE_SYSTEM,COORDINATE_RECTANGULAR_CARTESIAN_TYPE,ERR,ERROR,*999)
  CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,Num_Coordinates,ERR,ERROR,*999)
  CALL COORDINATE_SYSTEM_ORIGIN_SET(COORDINATE_SYSTEM,(/0.0_DP,0.0_DP,0.0_DP/),ERR,ERROR,*999) 
  CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR,ERROR,*999)

  !=CREATE REGION==================================================================================================================
  !Create Region and set CS to newly created 3D RC CS
  User_Num_REG = 1
  NULLIFY(REGION)
  CALL REGION_CREATE_START(User_Num_REG,WORLD_REGION,REGION,ERR,ERROR,*999)
  CALL REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*999)
  CALL REGION_CREATE_FINISH(REGION,ERR,ERROR,*999)

  !=CREATE BASIS ==================================================================================================================
  !Define basis - tri-linear Lagrange  
  User_Num_BASIS = 1
  Num_xi_coordinates = 3
  NULLIFY(BASIS)
  CALL BASIS_CREATE_START(User_Num_BASIS,BASIS,ERR,ERROR,*999) 
  CALL BASIS_TYPE_SET(BASIS,BASIS_LAGRANGE_HERMITE_TP_TYPE,ERR,ERROR,*999)
  CALL BASIS_NUMBER_OF_XI_SET(BASIS,Num_xi_coordinates,ERR,ERROR,*999)
  CALL BASIS_INTERPOLATION_XI_SET(BASIS,(/BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & BASIS_LINEAR_LAGRANGE_INTERPOLATION,BASIS_LINEAR_LAGRANGE_INTERPOLATION/),ERR,ERROR,*999)
  CALL BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET(BASIS,(/3,3,3/),ERR,ERROR,*999)  
  CALL BASIS_CREATE_FINISH(BASIS,ERR,ERROR,*999)

  USE_GENERATED_MESHES = .TRUE.

  IF(USE_GENERATED_MESHES) THEN
    !=AUTOMATIC MESH CREATION======================================================================================================
    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," *** USING GENERATED MESHES  ***",ERR,ERROR,*999)

    !=CREATE GENERATED MESH========================================================================================================
    !Start the creation of a generated mesh in the region
    NULLIFY(GENERATED_MESH)
    CALL GENERATED_MESH_CREATE_START(1,REGION,GENERATED_MESH,ERR,ERROR,*999)
    !Set up a regular 100x100 mesh
    CALL GENERATED_MESH_TYPE_SET(GENERATED_MESH,1,ERR,ERROR,*999)
    CALL GENERATED_MESH_BASIS_SET(GENERATED_MESH,BASIS,ERR,ERROR,*999)
    !Define the mesh on the region
    l = 120.0_DP
    w = 160.0_DP
    h = 10.0_DP
    CALL GENERATED_MESH_EXTENT_SET(GENERATED_MESH,(/l,w,h/),ERR,ERROR,*999)
    CALL GENERATED_MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH,(/NUMBER_GLOBAL_X_ELEMENTS, &
        & NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS/), ERR,ERROR,*999)
    !Finish the creation of a generated mesh in the region
    CALL GENERATED_MESH_CREATE_FINISH(GENERATED_MESH,1,MESH,ERR,ERROR,*999) 

    !=CREATE DECOMPOSITION=========================================================================================================
    !Create a decomposition
    NULLIFY(DECOMPOSITION)
    User_Num_Decomp = 1
    CALL DECOMPOSITION_CREATE_START(User_Num_Decomp,MESH,DECOMPOSITION,ERR,ERROR,*999)
    !Set the decomposition to be a general decomposition with the specified number of domains
    CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)
    CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*999)
    CALL DECOMPOSITION_CREATE_FINISH(MESH,DECOMPOSITION,ERR,ERROR,*999)

    !=CREATE GEOMETRIC FIELD=======================================================================================================
    !Start to create a default (geometric) field on the region
    NULLIFY(GEOMETRIC_FIELD)
    GEOMETRIC_FIELD_USER_NUM = 1  
    CALL FIELD_CREATE_START(GEOMETRIC_FIELD_USER_NUM,REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)
    !Set the decomposition to use
    CALL FIELD_MESH_DECOMPOSITION_SET(GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*999)
    !Set the domain to be used by the field components
    !NB these are needed now as the default mesh component number is 1
    Field_Variable = FIELD_U_VARIABLE_TYPE
    Mesh_component = 1
    Num_Field_Components = 3
    DO Component_idx=1,Num_Field_Components 
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,Field_Variable,Component_idx,Mesh_component,ERR,ERROR,*999)
    ENDDO !Component_idx
    !Finish creating the field
    CALL FIELD_CREATE_FINISH(GEOMETRIC_FIELD,ERR,ERROR,*999)
    !Update the geometric field parameters
    CALL GENERATED_MESH_GEOMETRIC_PARAMETERS_CALCULATE(GEOMETRIC_FIELD,GENERATED_MESH,ERR,ERROR,*999)

  ELSE
    !=MANUAL MESH CREATION=========================================================================================================
    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," *** USING MANUAL MESH CREATION ***",ERR,ERROR,*999)

    !=NODE CREATION================================================================================================================
    !Create 8 nodes in REGION and set initial coordinates to 0,0,0
    Num_Nodes = 8
    NULLIFY(NODES)
    CALL NODES_CREATE_START(REGION,Num_Nodes,NODES,ERR,ERROR,*999)
    CALL NODES_CREATE_FINISH(NODES,ERR,ERROR,*999)
  
    !=CREATE MESH==================================================================================================================
    !Create a general mesh with 1 component
    User_Num_Mesh = 1
    Num_mesh_dim = 3
    Num_mesh_components = 1
    Num_elem = 1
    NULLIFY(MESH)
    CALL MESH_CREATE_START(User_Num_Mesh,REGION,Num_mesh_dim,MESH,ERR,ERROR,*999)    
    CALL MESH_NUMBER_OF_COMPONENTS_SET(MESH,Num_mesh_components,ERR,ERROR,*999) 
    CALL MESH_NUMBER_OF_ELEMENTS_SET(MESH,Num_elem,ERR,ERROR,*999)  
    !Populate single mesh component with 1 element whose nodes are associated with REGION
    Mesh_component = 1
    Global_Element = 1
    NULLIFY(ELEMENTS)
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,Mesh_component,BASIS,ELEMENTS,ERR,ERROR,*999)
    CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(Global_Element,ELEMENTS,(/1,2,3,4,5,6,7,8/),ERR,ERROR,*999)
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH,Mesh_component,ERR,ERROR,*999)
    CALL MESH_CREATE_FINISH(MESH,ERR,ERROR,*999) 
  
    !=CREATE DECOMPOSITION=========================================================================================================
    !Create mesh decomposition dividing mesh into number_of_domains for parallel solving
    User_Num_Decomp = 1
    NULLIFY(DECOMPOSITION)
    CALL DECOMPOSITION_CREATE_START(User_Num_Decomp,MESH,DECOMPOSITION,ERR,ERROR,*999)
    CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)
    CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,number_of_domains,ERR,ERROR,*999)
    CALL DECOMPOSITION_CREATE_FINISH(MESH,DECOMPOSITION,ERR,ERROR,*999)

    !=CREATE GEOMETRIC FIELD=======================================================================================================
    !Start to create a default (geometric) field on the region
    NULLIFY(GEOMETRIC_FIELD)
    GEOMETRIC_FIELD_USER_NUM = 1  
    Num_Field_Variables = 1
    Num_Field_Components = 3
    CALL FIELD_CREATE_START(GEOMETRIC_FIELD_USER_NUM,REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)
    !Set the decomposition to use
    CALL FIELD_MESH_DECOMPOSITION_SET(GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*999)
    CALL FIELD_TYPE_SET(GEOMETRIC_FIELD,FIELD_GEOMETRIC_TYPE,ERR,ERROR,*999)  
    CALL FIELD_NUMBER_OF_VARIABLES_SET(GEOMETRIC_FIELD,Num_Field_Variables,ERR,ERROR,*999)
    CALL FIELD_NUMBER_OF_COMPONENTS_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,Num_Field_Components,ERR,ERROR,*999) 
    !Set the domain to be used by the field components
    !NB these are needed now as the default mesh component number is 1
    Field_Variable = 1
    Mesh_component = 1
    DO Component_idx=1,Num_Field_Components 
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,Field_Variable,Component_idx,Mesh_component,ERR,ERROR,*999)
    ENDDO !Component_idx
    CALL FIELD_CREATE_FINISH(GEOMETRIC_FIELD,ERR,ERROR,*999)
    ALLOCATE(Node_Coordinates(Num_Nodes,3),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Node_Coordinates",ERR,ERROR,*999)
    l = 120.0_DP
    w = 160.0_DP
    h = 10.0_DP
    Node_Coordinates = 0.0_DP
    Node_Coordinates(1:Num_Nodes,1) = (/0.0_DP,l,0.0_DP,l,0.0_DP,l,0.0_DP,l/) ! x coordinates for single 8 node element
    Node_Coordinates(1:Num_Nodes,2) = (/0.0_DP,0.0_DP,w,w,0.0_DP,0.0_DP,w,w/) ! y coordinates for single 8 node element
    Node_Coordinates(1:Num_Nodes,3) = (/0.0_DP,0.0_DP,0.0_DP,0.0_DP,h,h,h,h/) ! z coordinates for single 8 node element
    Field_Variable = 1 ! Displacment is the 1st & only Variable of the GEOMETRIC_FIELD
    Field_Derivative = 1
    DO np=1,Num_Nodes
      DO Component_idx=1,Num_Field_Components
      CALL FIELD_PARAMETER_SET_UPDATE_NODE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,Field_Derivative,np, &
        & Component_idx,Node_Coordinates(np,Component_idx),ERR,ERROR,*999)
      ENDDO !Component_idx
    ENDDO !np

  ENDIF

  !=CREATE DEPENDENT FIELD=========================================================================================================
  !Create a dependent field
  NULLIFY(DEPENDENT_FIELD)
  DEPENDENT_FIELD_USER_NUM = 10
  Mesh_component = 1
  Num_Field_Variables = 2
  Num_Field_Components = 3
  CALL FIELD_CREATE_START(DEPENDENT_FIELD_USER_NUM,REGION,DEPENDENT_FIELD,ERR,ERROR,*999)
  CALL FIELD_TYPE_SET(DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
  CALL FIELD_MESH_DECOMPOSITION_SET(DEPENDENT_FIELD,DECOMPOSITION,ERR,ERROR,*999)
  CALL FIELD_GEOMETRIC_FIELD_SET(DEPENDENT_FIELD,GEOMETRIC_FIELD,ERR,ERROR,*999)
  CALL FIELD_DEPENDENT_TYPE_SET(DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
  CALL FIELD_NUMBER_OF_VARIABLES_SET(DEPENDENT_FIELD,Num_Field_Variables,ERR,ERROR,*999)
  CALL FIELD_NUMBER_OF_COMPONENTS_SET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,Num_Field_Components,ERR,ERROR,*999)
  CALL FIELD_NUMBER_OF_COMPONENTS_SET(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,Num_Field_Components,ERR,ERROR,*999)
  DO Component_idx=1,Num_Field_Components
    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,Component_idx,Mesh_component,ERR,ERROR,*999)
    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,Component_idx,Mesh_component, &
        & ERR,ERROR,*999)
  ENDDO !Component_idx
  CALL FIELD_CREATE_FINISH(DEPENDENT_FIELD,ERR,ERROR,*999) 

  !=CREATE MATERIAL FIELD==========================================================================================================
  !Create a material field for a general 3D orthotropic material
  NULLIFY(MATERIAL_FIELD)
  MATERIAL_FIELD_USER_NUM = 3
  Num_Field_Variables = 1
  Num_Field_Components = 6
  CALL FIELD_CREATE_START(MATERIAL_FIELD_USER_NUM,REGION,MATERIAL_FIELD,ERR,ERROR,*999)
  CALL FIELD_TYPE_SET(MATERIAL_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
  CALL FIELD_MESH_DECOMPOSITION_SET(MATERIAL_FIELD,DECOMPOSITION,ERR,ERROR,*999)        
  CALL FIELD_GEOMETRIC_FIELD_SET(MATERIAL_FIELD,GEOMETRIC_FIELD,ERR,ERROR,*999)
  CALL FIELD_NUMBER_OF_VARIABLES_SET(MATERIAL_FIELD,Num_Field_Variables,ERR,ERROR,*999)
  CALL FIELD_NUMBER_OF_COMPONENTS_SET(MATERIAL_FIELD,FIELD_U_VARIABLE_TYPE,Num_Field_Components,ERR,ERROR,*999)
  DO Variable_idx=1,Num_Field_Variables
    DO Component_idx=1,Num_Field_Components
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(MATERIAL_FIELD,Variable_idx,Component_idx,Mesh_component,ERR,ERROR,*999)
    ENDDO !Component_idx
  ENDDO !Variable_idx
  CALL FIELD_CREATE_FINISH(MATERIAL_FIELD,ERR,ERROR,*999)
  E = (/30E6_DP,30E6_DP,30E6_DP/)
  v = (/0.25_DP,0.25_DP,0.25_DP/)
  E_v = (/E,v/)
  Field_Derivative = 1
  Num_nodes=mesh%TOPOLOGY(Mesh_component)%PTR%NODES%NUMBER_OF_NODES
  DO np=1,Num_nodes
    DO Component_idx=1,Num_Field_Components
      CALL FIELD_PARAMETER_SET_UPDATE_NODE(MATERIAL_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,Field_Derivative,np, &
        & Component_idx,E_v(Component_idx),ERR,ERROR,*999)
    ENDDO !Component_idx
  ENDDO !np

  !=CREATE EQUATION SET============================================================================================================
  !Create a Elasticity Class, Linear Elasticity type, no subtype, equations_set
  User_Num_EQUATION_SET = 1
  NULLIFY(EQUATIONS_SET)
  CALL EQUATIONS_SET_CREATE_START(User_Num_EQUATION_SET,REGION,GEOMETRIC_FIELD,EQUATIONS_SET,ERR,ERROR,*999)
  !Set the equations set to be a Elasticity Class, Linear Elasticity type, no subtype, equations_set
  CALL EQUATIONS_SET_SPECIFICATION_SET(EQUATIONS_SET,EQUATIONS_SET_ELASTICITY_CLASS,EQUATIONS_SET_LINEAR_ELASTICITY_TYPE, &
        EQUATIONS_SET_THREE_DIMENSIONAL_LINEAR_ELASTICITY_SUBTYPE,ERR,ERROR,*999)
  CALL EQUATIONS_SET_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
  !Create the equations set dependent field variables
  CALL EQUATIONS_SET_DEPENDENT_CREATE_START(EQUATIONS_SET,DEPENDENT_FIELD_USER_NUM,DEPENDENT_FIELD,ERR,ERROR,*999)
  CALL EQUATIONS_SET_DEPENDENT_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
  !Create the equations set material field variables
  CALL EQUATIONS_SET_MATERIALS_CREATE_START(EQUATIONS_SET,MATERIAL_FIELD_USER_NUM,MATERIAL_FIELD,ERR,ERROR,*999)  
  CALL EQUATIONS_SET_MATERIALS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
  !Create the equations set equations
  NULLIFY(EQUATIONS)
  CALL EQUATIONS_SET_EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
  CALL EQUATIONS_SPARSITY_TYPE_SET(EQUATIONS,EQUATIONS_SPARSE_MATRICES,ERR,ERROR,*999)
                                            !EQUATIONS_SPARSE_MATRICES=1 !<Use sparse matrices for the equations.
                                            !EQUATIONS_FULL_MATRICES=2 !<Use fully populated matrices for the equations. 
  CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_TIMING_OUTPUT,ERR,ERROR,*999)
                                          !EQUATIONS_ELEMENT_MATRIX_OUTPUT=3 !<All below and element matrices output.
                                          !EQUATIONS_MATRIX_OUTPUT=2 !<All below and equation matrices output.
                                          !EQUATIONS_TIMING_OUTPUT=1 !<Timing information output.
                                          !EQUATIONS_NO_OUTPUT=0 !<No output.
  CALL EQUATIONS_SET_EQUATIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999) 

  !=PRESCRIBE BOUNDARY CONDITIONS==================================================================================================
  !Prescribe Displacment to variable 1 of the Equation set's DEPENDENT_FIELD
  DISP_BC%NUM_BC_NODES = (/4,4,4/)
  IF (SUM(DISP_BC%NUM_BC_NODES) == 0) THEN
    CALL FLAG_ERROR("No DISP Boundary Conditions Applied.",ERR,ERROR,*999)
  ENDIF
  ALLOCATE(DISP_BC%BC_INFO(SUM(DISP_BC%NUM_BC_NODES),2),STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate BC_NODE.",ERR,ERROR,*999)
  ALLOCATE(DISP_BC%BC_VALUE(SUM(DISP_BC%NUM_BC_NODES)),STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate BC_VALUE.",ERR,ERROR,*999)
  IF (DISP_BC%NUM_BC_NODES(1) > 0) THEN 
    DISP_BC%BC_INFO(1:DISP_BC%NUM_BC_NODES(1),1) = (/1,3,5,7/) !Nodes
    DISP_BC%BC_INFO(1:DISP_BC%NUM_BC_NODES(1),2) = 1 !Component Number
    DISP_BC%BC_VALUE(1:DISP_BC%NUM_BC_NODES(1)) = 0.0_DP !Value
  ENDIF
  IF (DISP_BC%NUM_BC_NODES(2) > 0) THEN 
    DISP_BC%BC_INFO(DISP_BC%NUM_BC_NODES(1)+1:SUM(DISP_BC%NUM_BC_NODES(1:2)),1) = (/1,3,5,7/) !Nodes
    DISP_BC%BC_INFO(DISP_BC%NUM_BC_NODES(1)+1:SUM(DISP_BC%NUM_BC_NODES(1:2)),2) = 2 !Component Number
    DISP_BC%BC_VALUE(DISP_BC%NUM_BC_NODES(1)+1:SUM(DISP_BC%NUM_BC_NODES(1:2))) = 0.0_DP !Value
  ENDIF
  IF (DISP_BC%NUM_BC_NODES(3) > 0) THEN 
    DISP_BC%BC_INFO(SUM(DISP_BC%NUM_BC_NODES(1:2))+1:SUM(DISP_BC%NUM_BC_NODES(1:3)),1) = (/1,3,5,7/) !Nodes
    DISP_BC%BC_INFO(SUM(DISP_BC%NUM_BC_NODES(1:2))+1:SUM(DISP_BC%NUM_BC_NODES(1:3)),2) = 3 !Component Number
    DISP_BC%BC_VALUE(SUM(DISP_BC%NUM_BC_NODES(1:2))+1:SUM(DISP_BC%NUM_BC_NODES(1:3))) = 0.0_DP !Value
  ENDIF
  !Prescribe Force BC to variable 2 of the Equation set's DEPENDENT_FIELD
  FORCE_BC%NUM_BC_NODES = (/4,0,0/)
  IF (SUM(FORCE_BC%NUM_BC_NODES) == 0) THEN
    CALL FLAG_ERROR("No force Boundary Conditions Applied.",ERR,ERROR,*999)
  ENDIF
  ALLOCATE(FORCE_BC%BC_INFO(SUM(FORCE_BC%NUM_BC_NODES),2),STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate BC_NODE.",ERR,ERROR,*999)
  ALLOCATE(FORCE_BC%BC_VALUE(SUM(FORCE_BC%NUM_BC_NODES)),STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate BC_VALUE.",ERR,ERROR,*999)
  IF (FORCE_BC%NUM_BC_NODES(1) > 0) THEN 
    FORCE_BC%BC_INFO(1:FORCE_BC%NUM_BC_NODES(1),1) = (/2,4,6,8/) !Nodes
    FORCE_BC%BC_INFO(1:FORCE_BC%NUM_BC_NODES(1),2) = 1 !Component Number
    FORCE_BC%BC_VALUE(1:FORCE_BC%NUM_BC_NODES(1)) = 400.0_DP !Value
  ENDIF
  IF (FORCE_BC%NUM_BC_NODES(2) > 0) THEN 
    FORCE_BC%BC_INFO(FORCE_BC%NUM_BC_NODES(1)+1:SUM(FORCE_BC%NUM_BC_NODES(1:2)),1) = 0 !Nodes
    FORCE_BC%BC_INFO(FORCE_BC%NUM_BC_NODES(1)+1:SUM(FORCE_BC%NUM_BC_NODES(1:2)),2) = 2 !Component Number
    FORCE_BC%BC_VALUE(FORCE_BC%NUM_BC_NODES(1)+1:SUM(FORCE_BC%NUM_BC_NODES(1:2))) = 0.0_DP !Value
  ENDIF
  IF (FORCE_BC%NUM_BC_NODES(3) > 0) THEN 
    FORCE_BC%BC_INFO(SUM(FORCE_BC%NUM_BC_NODES(1:2))+1:SUM(FORCE_BC%NUM_BC_NODES(1:3)),1) = 0 !Nodes
    FORCE_BC%BC_INFO(SUM(FORCE_BC%NUM_BC_NODES(1:2))+1:SUM(FORCE_BC%NUM_BC_NODES(1:3)),2) = 3 !Component Number
    FORCE_BC%BC_VALUE(SUM(FORCE_BC%NUM_BC_NODES(1:2))+1:SUM(FORCE_BC%NUM_BC_NODES(1:3))) = 0.0_DP !Value
  ENDIF

  Field_Derivative = 1 
  NULLIFY(BOUNDARY_CONDITIONS)
  CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
  Field_Variable = 1 ! Displacment is the 1st Variable of the DEPENDENT_FIELD
  DO BC_idx=1,SUM(DISP_BC%NUM_BC_NODES) !Loop over fixed BC nodes
    CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,FIELD_U_VARIABLE_TYPE,Field_Derivative, &
      & DISP_BC%BC_INFO(BC_idx,1),DISP_BC%BC_INFO(BC_idx,2),BOUNDARY_CONDITION_FIXED, &
      & DISP_BC%BC_VALUE(BC_idx),ERR,ERROR,*999)
  ENDDO !BC_idx
  Field_Variable = 2 ! Force is the 2nd Variable of the DEPENDENT_FIELD
  DO BC_idx = 1,SUM(FORCE_BC%NUM_BC_NODES)
!    CALL FIELD_PARAMETER_SET_UPDATE_NODE(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,Field_Derivative,FORCE_BC%BC_INFO(BC_idx,1), &
!      & FORCE_BC%BC_INFO(BC_idx,2),Field_Variable,FORCE_BC%BC_VALUE(BC_idx),ERR,ERROR,*999)
    CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,FIELD_DELUDELN_VARIABLE_TYPE,Field_Derivative, &
      & FORCE_BC%BC_INFO(BC_idx,1),FORCE_BC%BC_INFO(BC_idx,2),BOUNDARY_CONDITION_FIXED, &
      & FORCE_BC%BC_VALUE(BC_idx),ERR,ERROR,*999)
  ENDDO !BC_idx
  CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
  
!  !Prescribe Force BC to variable 2 of the Equation set's DEPENDENT_FIELD
!  FORCE_BC%NUM_BC_NODES = (/4,0,0/)
!  IF (SUM(FORCE_BC%NUM_BC_NODES) == 0) THEN
!    CALL FLAG_ERROR("No force Boundary Conditions Applied.",ERR,ERROR,*999)
!  ENDIF
!  ALLOCATE(FORCE_BC%BC_INFO(SUM(FORCE_BC%NUM_BC_NODES),2),STAT=ERR)
!  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate BC_NODE.",ERR,ERROR,*999)
!  ALLOCATE(FORCE_BC%BC_VALUE(SUM(FORCE_BC%NUM_BC_NODES)),STAT=ERR)
!  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate BC_VALUE.",ERR,ERROR,*999)
!  IF (FORCE_BC%NUM_BC_NODES(1) > 0) THEN 
!    FORCE_BC%BC_INFO(1:FORCE_BC%NUM_BC_NODES(1),1) = (/9,10,11,12/)
!    FORCE_BC%BC_INFO(1:FORCE_BC%NUM_BC_NODES(1),2) = 1
!    FORCE_BC%BC_VALUE(1:FORCE_BC%NUM_BC_NODES(1)) = 400.0_DP
!  ENDIF
!  IF (FORCE_BC%NUM_BC_NODES(2) > 0) THEN 
!    FORCE_BC%BC_INFO(FORCE_BC%NUM_BC_NODES(1)+1:SUM(FORCE_BC%NUM_BC_NODES(1:2)),1) = 0
!    FORCE_BC%BC_INFO(FORCE_BC%NUM_BC_NODES(1)+1:SUM(FORCE_BC%NUM_BC_NODES(1:2)),2) = 2
!    FORCE_BC%BC_VALUE(FORCE_BC%NUM_BC_NODES(1)+1:SUM(FORCE_BC%NUM_BC_NODES(1:2))) = 0.0_DP
!  ENDIF
!  IF (FORCE_BC%NUM_BC_NODES(3) > 0) THEN 
!    FORCE_BC%BC_INFO(SUM(FORCE_BC%NUM_BC_NODES(1:2))+1:SUM(FORCE_BC%NUM_BC_NODES(1:3)),1) = 0
!    FORCE_BC%BC_INFO(SUM(FORCE_BC%NUM_BC_NODES(1:2))+1:SUM(FORCE_BC%NUM_BC_NODES(1:3)),2) = 3
!    FORCE_BC%BC_VALUE(SUM(FORCE_BC%NUM_BC_NODES(1:2))+1:SUM(FORCE_BC%NUM_BC_NODES(1:3))) = 0.0_DP
!  ENDIF
!  Field_Variable = 2 ! Force is the 2nd Variable of the DEPENDENT_FIELD
!  Field_Derivative = 1
!  DO BC_idx = 1,SUM(FORCE_BC%NUM_BC_NODES)
!    CALL FIELD_PARAMETER_SET_UPDATE_NODE(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,Field_Derivative,FORCE_BC%BC_INFO(BC_idx,1), &
!      & FORCE_BC%BC_INFO(BC_idx,2),Field_Variable,FORCE_BC%BC_VALUE(BC_idx),ERR,ERROR,*999)
!  ENDDO !BC_idx

  !=CREATE PROBLEM=================================================================================================================
  !Create the problem
  NULLIFY(PROBLEM)
  User_Num_PROBLEM=1
  CALL PROBLEM_CREATE_START(User_Num_PROBLEM,PROBLEM,ERR,ERROR,*999)
  !Set the problem to be a elasticity class, linear elasticity type with no subtype.
  CALL PROBLEM_SPECIFICATION_SET(PROBLEM,PROBLEM_ELASTICITY_CLASS,PROBLEM_LINEAR_ELASTICITY_TYPE, &
    & PROBLEM_NO_SUBTYPE,ERR,ERROR,*999)
  CALL PROBLEM_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !=CREATE PROBLEM CONTROL LOOP====================================================================================================
  !Create the problem control loop
  CALL PROBLEM_CONTROL_LOOP_CREATE_START(PROBLEM,ERR,ERROR,*999)
  CALL PROBLEM_CONTROL_LOOP_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !=CREATE PROBLEM SOLVER==========================================================================================================
  !Start the creation of the problem solvers
  NULLIFY(SOLVER)
  Solver_Idx = 1
  CALL PROBLEM_SOLVERS_CREATE_START(PROBLEM,ERR,ERROR,*999)
  CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,Solver_Idx,SOLVER,ERR,ERROR,*999)
  !CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
                                     !SOLVER_CMISS_LIBRARY=LIBRARY_CMISS_TYPE !<CMISS (internal) solver library 
                                     !SOLVER_PETSC_LIBRARY=LIBRARY_PETSC_TYPE !<PETSc solver library
  !CALL SOLVER_LINEAR_TYPE_SET(SOLVER,SOLVER_LINEAR_DIRECT_SOLVE_TYPE,ERR,ERROR,*999)
                                    !SOLVER_LINEAR_DIRECT_SOLVE_TYPE=1 !<Direct linear solver type \NOT IMPLEMENTED
                                    !SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE=2 !<Iterative linear solver type
  !CALL SOLVER_LINEAR_DIRECT_TYPE_SET(SOLVER,SOLVER_DIRECT_LU,ERR,ERROR,*999)  
                                           !SOLVER_DIRECT_LU=1 !<LU direct linear solver 
                                           !SOLVER_DIRECT_CHOLESKY=2 !<Cholesky direct linear solver \NOT IMPLEMENTED
                                           !SOLVER_DIRECT_SVD=3 !<SVD direct linear solver \NOT IMPLEMENTED
  CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_SOLVER_OUTPUT,ERR,ERROR,*999)
                                    !SOLVER_MATRIX_OUTPUT=4 !<SolVER matrices output from the solver routines plus below
                                    !SOLVER_SOLVER_OUTPUT=3 !<Solver specific output from the solver routines plus below
                                    !SOLVER_TIMING_OUTPUT=2 !<Timing output from the solver routines plus below
                                    !SOLVER_PROGRESS_OUTPUT=1 !<Progress output from solver routines 
                                    !SOLVER_NO_OUTPUT=0 !<No output from the solver routines 
  CALL PROBLEM_SOLVERS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !=CREATE PROBLEM SOLVER EQUATIONS================================================================================================
  !Create the problem solver equations
  NULLIFY(SOLVER) !WHY NULLIFY HERE AFTER ALREADY NULLIFIED ABOVE????
  NULLIFY(SOLVER_EQUATIONS)
  CALL PROBLEM_SOLVER_EQUATIONS_CREATE_START(PROBLEM,ERR,ERROR,*999)
  CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,Solver_Idx,SOLVER,ERR,ERROR,*999)
  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                                                          !SOLVER_SPARSE_MATRICES=1 !<Use sparse solver matrices
                                                          !SOLVER_FULL_MATRICES=2 !<Use fully populated solver matrices
  Equations_Set_Idx = 0 !Initialize index of the equations set that has been added 
                        !(Variable is returned from PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD)
  CALL SOLVER_EQUATIONS_EQUATIONS_SET_ADD(SOLVER_EQUATIONS,EQUATIONS_SET,Equations_Set_Idx,ERR,ERROR,*999)
  !CALL SOLVER_MATRICES_STATIC_ASSEMBLE(SOLVER,SOLVER_MATRICES_LINEAR_ONLY,ERR,ERROR,*999)
                                         !SOLVER_MATRICES_ALL=1 !<Select all the solver matrices and vectors 
                                         !SOLVER_MATRICES_DYNAMIC_ONLY=2 !<Select only the dynamic solver matrices and vectors 
                                         !SOLVER_MATRICES_LINEAR_ONLY=3 !<Select only the linear solver matrices and vectors 
                                         !SOLVER_MATRICES_NONLINEAR_ONLY=4 !<Select only the nonlinear solver matrices and vectors
                                         !SOLVER_MATRICES_JACOBIAN_ONLY=5 !<Select only the Jacobian solver matrix
                                         !SOLVER_MATRICES_RESIDUAL_ONLY=6 !<Select only the residual solver vector
                                         !SOLVER_MATRICES_RHS_ONLY=7 !<Select only the RHS solver vector
                                         !SOLVER_MATRICES_RHS_RESIDUAL_ONLY=8 !<Select only the residual and RHS solver vectors
  CALL PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !=SOLVE PROBLEM==================================================================================================================
  !Solve the problem
  CALL PROBLEM_SOLVE(PROBLEM,ERR,ERROR,*999)

  !=OUTPUT SOLUTION================================================================================================================
  FILE="LinearElasticityExample"
  METHOD="FORTRAN"
  Export_Field=.TRUE.
  IF(Export_Field) THEN
    CALL FIELD_IO_NODES_EXPORT(REGION%FIELDS, FILE, METHOD, ERR,ERROR,*999)  
    CALL FIELD_IO_ELEMENTS_EXPORT(REGION%FIELDS, FILE, METHOD, ERR,ERROR,*999)
  ENDIF

  !================================================================================================================================

  !Calculate the stop times and write out the elapsed user and system times
  CALL CPU_TIMER(USER_CPU,STOP_USER_TIME,ERR,ERROR,*999)
  CALL CPU_TIMER(SYSTEM_CPU,STOP_SYSTEM_TIME,ERR,ERROR,*999)

  CALL WRITE_STRING_TWO_VALUE(GENERAL_OUTPUT_TYPE,"User time = ",STOP_USER_TIME(1)-START_USER_TIME(1),", System time = ", &
    & STOP_SYSTEM_TIME(1)-START_SYSTEM_TIME(1),ERR,ERROR,*999)

  CALL CMISS_FINALISE(ERR,ERROR,*999)

  WRITE(*,'(A)') "Program successfully completed."

  STOP
999 CALL CMISS_WRITE_ERROR(ERR,ERROR)
  STOP

END PROGRAM LINEARELASTICITYEXAMPLE
