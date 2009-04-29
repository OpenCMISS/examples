!> \file
!> $Id: StokesFlowExample.f90 20 2009-04-08 20:22:52Z cpb $
!> \author Sebastian Krittian
!> \brief This is an example program to solve a Stokes equation using openCMISS calls.
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
!> The Original Code is openCMISS
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

!> Main program

PROGRAM StokesFlow



   USE BASE_ROUTINES				! needed
   USE BASIS_ROUTINES				! needed
   USE CMISS					! needed
   USE CMISS_MPI				! needed with MPI
   USE COMP_ENVIRONMENT				! needed with MPI
   USE CONSTANTS				! basic definitions: pi, e, I, ...
   USE CONTROL_LOOP_ROUTINES			! needed
   USE COORDINATE_ROUTINES			! needed
!   USE DISTRIBUTED_MATRIX_VECTOR		! not needed
   USE DOMAIN_MAPPINGS				! needed
   USE EQUATIONS_ROUTINES			! needed
   USE EQUATIONS_SET_CONSTANTS			! needed
   USE EQUATIONS_SET_ROUTINES			! needed
   USE FIELD_ROUTINES				! needed
   USE FIELD_IO_ROUTINES			! not needed	
!   USE GENERATED_MESH_ROUTINES			! not needed
   USE INPUT_OUTPUT				! needed
   USE ISO_VARYING_STRING			! needed
   USE KINDS					! needed
!   USE LISTS					! needed	
   USE MESH_ROUTINES				! needed
   USE MPI					! needed with MPI
   USE PROBLEM_CONSTANTS			! needed
   USE PROBLEM_ROUTINES				! needed
   USE REGION_ROUTINES				! needed
   USE SOLVER_ROUTINES				! needed
   USE TIMER					! needed
   USE TYPES					! needed
 
#ifdef WIN32
   USE IFQWIN
#endif
 
!use module in order to include cmheart input file for nodes and elements
  USE IMPORT_CMHEART

IMPLICIT NONE
  

  !define CMHeart_ImportFiles
    TYPE(EXPORT_CONTAINER):: CM
 
   !Program types
   
   !Program variables
   
   !Program variables

  INTEGER(INTG) :: NUMBER_OF_DOMAINS
  
  INTEGER(INTG) :: MPI_IERROR

  INTEGER(INTG) :: EQUATIONS_SET_INDEX
  TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOF_MAPPING
  
  TYPE(BASIS_TYPE), POINTER :: BASIS_M,BASIS_V,BASIS_P
  TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
!  TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH		! not needed
  TYPE(MESH_TYPE), POINTER :: MESH
  TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
  TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
  TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
  TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD
  TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
  TYPE(REGION_TYPE), POINTER :: REGION
  TYPE(SOLVER_TYPE), POINTER :: SOLVER
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
   
  LOGICAL :: EXPORT_FIELD,IMPORT_FIELD
  TYPE(VARYING_STRING) :: FILE,METHOD

  REAL(SP) :: START_USER_TIME(1),STOP_USER_TIME(1),START_SYSTEM_TIME(1),STOP_SYSTEM_TIME(1)

   
   INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES
   INTEGER(INTG) :: MY_COMPUTATIONAL_NODE_NUMBER
  
   TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS_M,MESH_ELEMENTS_P,MESH_ELEMENTS_V
   TYPE(NODES_TYPE), POINTER :: NODES
   
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

   !Generic CMISS variables
   
   INTEGER(INTG) :: ERR
   TYPE(VARYING_STRING) :: ERROR
 
   INTEGER(INTG) :: DIAG_LEVEL_LIST(5)
   CHARACTER(LEN=MAXSTRLEN) :: DIAG_ROUTINE_LIST(1),TIMING_ROUTINE_LIST(1)
   
   INTEGER k,l,m,n


   CALL READ_CMHEART_EXE
   CALL RECV_CMHEART_EXE(CM)


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
  
   !Intialise cmiss
   CALL CMISS_INITIALISE(ERR,ERROR,*999)
  
  !Set all diganostic levels on for testing
  !DIAG_LEVEL_LIST(1)=1
  !DIAG_LEVEL_LIST(2)=2
  !DIAG_LEVEL_LIST(3)=3
  !DIAG_LEVEL_LIST(4)=4
  !DIAG_LEVEL_LIST(5)=5
  !DIAG_ROUTINE_LIST(1)=""
  !CALL DIAGNOSTICS_SET_ON(ALL_DIAG_TYPE,DIAG_LEVEL_LIST,"MoreComplexMeshExample",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
  !CALL DIAGNOSTICS_SET_ON(ALL_DIAG_TYPE,DIAG_LEVEL_LIST,"",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
 
  !TIMING_ROUTINE_LIST(1)=""
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)
  
  !Calculate the start times
  CALL CPU_TIMER(USER_CPU,START_USER_TIME,ERR,ERROR,*999)
  CALL CPU_TIMER(SYSTEM_CPU,START_SYSTEM_TIME,ERR,ERROR,*999)
  
  !Get the number of computational nodes
  NUMBER_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
  IF(ERR/=0) GOTO 999
  !Get my computational node number
  MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
  IF(ERR/=0) GOTO 999


  ! HERE IT'S GETTING INTERESTING

  ! SAME FOR STOKES
  !Start the creation of a new RC coordinate system
  NULLIFY(COORDINATE_SYSTEM)
  CALL COORDINATE_SYSTEM_CREATE_START(1,COORDINATE_SYSTEM,ERR,ERROR,*999)
  !Set the coordinate system to be DIMENSIONS
  CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,CM%D,ERR,ERROR,*999)
  !Finish the creation of the coordinate system
  CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR,ERROR,*999)



  ! SAME FOR STOKES
  !Start the creation of a region
  NULLIFY(REGION)
  CALL REGION_CREATE_START(1,REGION,ERR,ERROR,*999)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*999)
  !Finish the creation of the region
  CALL REGION_CREATE_FINISH(REGION,ERR,ERROR,*999)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! START CAREFUL BASIS: 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! -> define *_P and *_V and *_M in input file first

  !Start the creation of a mesh basis
  NULLIFY(BASIS_M)
  CALL BASIS_CREATE_START(CM%ID_M,BASIS_M,ERR,ERROR,*999)  
  !Set the basis to CM%IT_M (1=linear,2=quadratic,3=cubic)
  CALL BASIS_NUMBER_OF_XI_SET(BASIS_M,CM%D,ERR,ERROR,*999)
  IF (CM%D==2) THEN
    CALL BASIS_INTERPOLATION_XI_SET(BASIS_M,(/CM%IT_M,CM%IT_M/),ERR,ERROR,*999)
  ELSE IF (CM%D==3) THEN
    CALL BASIS_INTERPOLATION_XI_SET(BASIS_M,(/CM%IT_M,CM%IT_M,CM%IT_M/),ERR,ERROR,*999)
  ELSE
    GOTO 999
  END IF
  CALL BASIS_CREATE_FINISH(BASIS_M,ERR,ERROR,*999)

  !Start the creation of a velocity basis
  NULLIFY(BASIS_V)
  CALL BASIS_CREATE_START(CM%ID_V,BASIS_V,ERR,ERROR,*999)  
  !Set the basis to CM%IT_V (1=linear,2=quadratic,3=cubic)
  CALL BASIS_NUMBER_OF_XI_SET(BASIS_V,CM%D,ERR,ERROR,*999)
  IF (CM%D==2) THEN
    CALL BASIS_INTERPOLATION_XI_SET(BASIS_V,(/CM%IT_V,CM%IT_V/),ERR,ERROR,*999)
  ELSE IF (CM%D==3) THEN
    CALL BASIS_INTERPOLATION_XI_SET(BASIS_V,(/CM%IT_V,CM%IT_V,CM%IT_V/),ERR,ERROR,*999)
  ELSE
    GOTO 999
  END IF
  CALL BASIS_CREATE_FINISH(BASIS_V,ERR,ERROR,*999)


  !Start the creation of a pressure basis
  NULLIFY(BASIS_P)
  CALL BASIS_CREATE_START(CM%ID_P,BASIS_P,ERR,ERROR,*999)  
  !Set the basis to CM%IT_P (1=linear,2=quadratic,3=cubic)
  CALL BASIS_NUMBER_OF_XI_SET(BASIS_P,CM%D,ERR,ERROR,*999)
  IF (CM%D==2) THEN
    CALL BASIS_INTERPOLATION_XI_SET(BASIS_P,(/CM%IT_P,CM%IT_P/),ERR,ERROR,*999)
  ELSE IF (CM%D==3) THEN
    CALL BASIS_INTERPOLATION_XI_SET(BASIS_P,(/CM%IT_P,CM%IT_P,CM%IT_P/),ERR,ERROR,*999)
  ELSE
    GOTO 999
  END IF
  CALL BASIS_CREATE_FINISH(BASIS_P,ERR,ERROR,*999)


  ! END CAREFUL BASIS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! START CAREFUL MESH
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! -> how to define different mesh components if nodes have to be in one array???
  ! 	e.g. cubic for V and quadratuc for P

  !Create a mesh
  NULLIFY(NODES)
  NULLIFY(MESH)

  NULLIFY(MESH_ELEMENTS_M)
  NULLIFY(MESH_ELEMENTS_P)
  NULLIFY(MESH_ELEMENTS_V)

  CALL NODES_CREATE_START(CM%N_T,REGION,NODES,ERR,ERROR,*999)
  CALL NODES_CREATE_FINISH(REGION,ERR,ERROR,*999)

  

! DEFINE GENERAL MESH ON REGION
    CALL MESH_CREATE_START(1,REGION,CM%D,MESH,ERR,ERROR,*999)
    CALL MESH_NUMBER_OF_ELEMENTS_SET(MESH,CM%E_T,ERR,ERROR,*999)
    CALL MESH_NUMBER_OF_COMPONENTS_SET(MESH,3,ERR,ERROR,*999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TO BE DONE:	cmh_NodesArray anpassen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!DEFINE MESH FOR MESH ELEMENTS
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,CM%ID_M,BASIS_M,MESH_ELEMENTS_M,ERR,ERROR,*999)
      DO k=1,CM%E_T
        CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(k,MESH_ELEMENTS_M, &
         CM%M(k,1:CM%EN_M),ERR,ERROR,*999)
      END DO
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH,CM%ID_M,ERR,ERROR,*999)


  !DEFINE MESH FOR VELOCITY ELEMENTS
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,CM%ID_V,BASIS_V,MESH_ELEMENTS_V,ERR,ERROR,*999)
      DO k=1,CM%E_T
        CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(k,MESH_ELEMENTS_V, &
         CM%V(k,1:CM%EN_V),ERR,ERROR,*999)
      END DO
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH,CM%ID_V,ERR,ERROR,*999)



!DEFINE MESH FOR PRESSURE ELEMENTS
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,CM%ID_P,BASIS_P,MESH_ELEMENTS_P,ERR,ERROR,*999)
      DO k=1,CM%E_T
        CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(k,MESH_ELEMENTS_P, &
         CM%P(k,1:CM%EN_P),ERR,ERROR,*999)
      END DO
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH,CM%ID_P,ERR,ERROR,*999)

    CALL MESH_CREATE_FINISH(REGION,MESH,ERR,ERROR,*999)
  
  ! END CAREFUL MESH


  !Create a decomposition for mesh
  ! SAME FOR STOKES

  NULLIFY(DECOMPOSITION)
  CALL DECOMPOSITION_CREATE_START(1,MESH,DECOMPOSITION,ERR,ERROR,*999)

      !Set the decomposition to be a general decomposition with the specified number of domains
      CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)
      CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)

  !Finish the decomposition creation
  CALL DECOMPOSITION_CREATE_FINISH(MESH,DECOMPOSITION,ERR,ERROR,*999)



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! START CAREFUL GEOMETRIC FIELD
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! -> how to define different mesh components if nodes have to be in one array???
  ! 	e.g. cubic for V and quadratuc for P

  !Start to create a default (geometric) field on the region
  NULLIFY(GEOMETRIC_FIELD)
  CALL FIELD_CREATE_START(1,REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)

      !Set the decomposition to use
      CALL FIELD_MESH_DECOMPOSITION_SET(GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*999)
      !Set the domain to be used by the field components

    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1,1,ERR,ERROR,*999)
    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,2,1,ERR,ERROR,*999)
    IF(CM%D==3) THEN
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,3,1,ERR,ERROR,*999)
    ENDIF

  !Finish creating the field
  CALL FIELD_CREATE_FINISH(REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)


  ! END CAREFUL GEOMETRIC FIELD

! !   !Set the geometric field values
!
!	auch hier noch CM%N_T anpassen fuer die verschiedenen Faelle
!
!  

! SAME FOR STOKE SO FAR! but be careful with the array given by input

  DO k=1,CM%N_M
    DO j=1,CM%D
      CALL FIELD_PARAMETER_SET_UPDATE_NODE(GEOMETRIC_FIELD,FIELD_VALUES_SET_TYPE,1,k,j, &
           & FIELD_U_VARIABLE_TYPE,CM%N(k,j),ERR,ERROR,*999)
    END DO
  END DO
  
  CALL FIELD_PARAMETER_SET_UPDATE_START(GEOMETRIC_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(GEOMETRIC_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
  

  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD)) GEOMETRIC_FIELD=>REGION%FIELDS%FIELDS(1)%PTR
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   HIER GEHTS LOS UND DA SIMMA DABEI !!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!	HIER WIRD GLAUB ICH RICHTIG INTERESSANT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!->->-> EQUATION SET

  !Create the equations_set
  NULLIFY(EQUATIONS_SET)

  !!! hier passiert anscheinend noch nix
  CALL EQUATIONS_SET_CREATE_START(1,REGION,GEOMETRIC_FIELD,EQUATIONS_SET,ERR,ERROR,*999)
  !Set the equations set to be a Stokes Flow problem

  !!! SUBROUTINE TO BE CHANGED: SPECIFICATION CHANGES  
  CALL EQUATIONS_SET_SPECIFICATION_SET(EQUATIONS_SET,EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_STOKES_FLUID_TYPE, &
    & EQUATIONS_SET_STANDARD_STOKES_SUBTYPE,ERR,ERROR,*999)
  !Finish creating the equations set

  !!! SUBROUTINE TO BE CHANGED: SELECT SET CHANGES
  CALL EQUATIONS_SET_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
! ! ! ! ! ! ! 
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !->->-> DEPENDENT FIELD: HIER MUESSSEN ALLE ABHANEGIGKEITEN GEKLAERT WERDEN: HAUPTSACHE CASE SELECTS
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !-> define with 4 components: p,u,v,w !!!!!!!
! ! ! ! ! ! ! 
! ! ! ! ! ! !   !Create the equations set dependent field variables
! ! ! ! ! ! !   CALL EQUATIONS_SET_DEPENDENT_CREATE_START(EQUATIONS_SET,ERR,ERROR,*999)
! ! ! ! ! ! !   !Finish the equations set dependent field variables
! ! ! ! ! ! !   CALL EQUATIONS_SET_DEPENDENT_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
! ! ! ! ! ! ! 
! ! ! ! ! ! !   !Temporary until fix fixed conditions set
! ! ! ! ! ! !   !Find the first and last dof numbers and ranks
! ! ! ! ! ! !   DEPENDENT_DOF_MAPPING=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
! ! ! ! ! ! ! 
! ! ! ! ! ! ! 
! ! ! ! ! ! !   !Create the problem fixed conditions
! ! ! ! ! ! !   CALL EQUATIONS_SET_FIXED_CONDITIONS_CREATE_START(EQUATIONS_SET,ERR,ERROR,*999)
! ! ! ! ! ! !   !Set bc's
! ! ! ! ! ! !   CALL EQUATIONS_SET_FIXED_CONDITIONS_SET_DOF(EQUATIONS_SET,1,EQUATIONS_SET_FIXED_BOUNDARY_CONDITION, &
! ! ! ! ! ! !     & 42.0_DP,ERR,ERROR,*999)
! ! ! ! ! ! !   !Finish the problem fixed conditions
! ! ! ! ! ! !   CALL EQUATIONS_SET_FIXED_CONDITIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! ! ! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! ! ! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! HIER GEHTS ANS EINGEMACHTE -> ANPASSEN DER MATRIZEN UND DEFINITION LAGRANGE.F90
! ! ! ! ! ! ! !  
! ! ! ! ! ! ! 
! ! ! ! ! ! !   !Create the equations set equations
! ! ! ! ! ! !   NULLIFY(EQUATIONS)
! ! ! ! ! ! !   CALL EQUATIONS_SET_EQUATIONS_CREATE_START(EQUATIONS_SET,ERR,ERROR,*999)
! ! ! ! ! ! !   !Get the equations
! ! ! ! ! ! !   CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
! ! ! ! ! ! !   !Set the equations matrices sparsity type
! ! ! ! ! ! !   CALL EQUATIONS_SPARSITY_TYPE_SET(EQUATIONS,EQUATIONS_SPARSE_MATRICES,ERR,ERROR,*999)
! ! ! ! ! ! !   !CALLEQUATIONS_SPARSITY_TYPE_SET(EQUATIONS,EQUATIONS_FULL_MATRICES,ERR,ERROR,*999)
! ! ! ! ! ! !   !Set the equations set output
! ! ! ! ! ! !   !CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_TIMING_OUTPUT,ERR,ERROR,*999)
! ! ! ! ! ! !   !CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_MATRIX_OUTPUT,ERR,ERROR,*999)
! ! ! ! ! ! !   CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_ELEMENT_MATRIX_OUTPUT,ERR,ERROR,*999)
! ! ! ! ! ! !   CALL EQUATIONS_SET_EQUATIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999) 
! ! ! ! ! ! ! 
! ! ! ! ! ! !   !Create the problem
! ! ! ! ! ! !   NULLIFY(PROBLEM)
! ! ! ! ! ! !   CALL PROBLEM_CREATE_START(1,PROBLEM,ERR,ERROR,*999)
! ! ! ! ! ! !   !Set the problem to be a standard Stokes problem
! ! ! ! ! ! !   CALL PROBLEM_SPECIFICATION_SET(PROBLEM,PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_STOKES_FLUID_TYPE, &
! ! ! ! ! ! !     & PROBLEM_STANDARD_STOKES_SUBTYPE,ERR,ERROR,*999)
! ! ! ! ! ! !   !Finish creating the problem
! ! ! ! ! ! !   CALL PROBLEM_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)
! ! ! ! ! ! ! 
! ! ! ! ! ! !   !Create the problem control loop
! ! ! ! ! ! !   CALL PROBLEM_CONTROL_LOOP_CREATE_START(PROBLEM,ERR,ERROR,*999)
! ! ! ! ! ! !   !Finish creating the problem control loop
! ! ! ! ! ! !   CALL PROBLEM_CONTROL_LOOP_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)
! ! ! ! ! ! ! 
! ! ! ! ! ! !   !Start the creation of the problem solvers
! ! ! ! ! ! !   NULLIFY(SOLVER)
! ! ! ! ! ! !   CALL PROBLEM_SOLVERS_CREATE_START(PROBLEM,ERR,ERROR,*999)
! ! ! ! ! ! !   CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,1,SOLVER,ERR,ERROR,*999)
! ! ! ! ! ! !   !CALL SOLVER_LINEAR_TYPE_SET(SOLVER,SOLVER_LINEAR_DIRECT_SOLVE_TYPE,ERR,ERROR,*999)
! ! ! ! ! ! !   !CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_NO_OUTPUT,ERR,ERROR,*999)
! ! ! ! ! ! !   !CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_TIMING_OUTPUT,ERR,ERROR,*999)
! ! ! ! ! ! !   !CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_SOLVER_OUTPUT,ERR,ERROR,*999)
! ! ! ! ! ! !   CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_MATRIX_OUTPUT,ERR,ERROR,*999)
! ! ! ! ! ! !   !Finish the creation of the problem solver
! ! ! ! ! ! !   CALL PROBLEM_SOLVERS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)
! ! ! ! ! ! ! 
! ! ! ! ! ! !   !Create the problem solver equations
! ! ! ! ! ! !   NULLIFY(SOLVER)
! ! ! ! ! ! !   NULLIFY(SOLVER_EQUATIONS)
! ! ! ! ! ! !   CALL PROBLEM_SOLVER_EQUATIONS_CREATE_START(PROBLEM,ERR,ERROR,*999)
! ! ! ! ! ! !   CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,1,SOLVER,ERR,ERROR,*999)
! ! ! ! ! ! !   CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
! ! ! ! ! ! !   CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
! ! ! ! ! ! !   !CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER,SOLVER_FULL_MATRICES,ERR,ERROR,*999)
! ! ! ! ! ! !   !Add in the equations set
! ! ! ! ! ! !   CALL PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD(PROBLEM,CONTROL_LOOP_NODE,1,EQUATIONS_SET,EQUATIONS_SET_INDEX,ERR,ERROR,*999)
! ! ! ! ! ! !   !Finish the problem solver equations
! ! ! ! ! ! !   CALL PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)
! ! ! ! ! ! ! 
! ! ! ! ! ! !   !Solve the problem
! ! ! ! ! ! !   CALL PROBLEM_SOLVE(PROBLEM,ERR,ERROR,*999)
! ! ! ! ! ! ! 
   EXPORT_FIELD=.TRUE.
   METHOD="FORTRAN"
   IF(EXPORT_FIELD) THEN
     FILE="StokesFlowExample"
     CALL FIELD_IO_NODES_EXPORT(REGION%FIELDS, FILE, METHOD, ERR,ERROR,*999)  
     CALL FIELD_IO_ELEMENTS_EXPORT(REGION%FIELDS, FILE, METHOD, ERR,ERROR,*999)
   ENDIF
   
   !Calculate the stop times and write out the elapsed user and system times
   CALL CPU_TIMER(USER_CPU,STOP_USER_TIME,ERR,ERROR,*999)
   CALL CPU_TIMER(SYSTEM_CPU,STOP_SYSTEM_TIME,ERR,ERROR,*999)
 
   CALL WRITE_STRING_TWO_VALUE(GENERAL_OUTPUT_TYPE,"User time = ",STOP_USER_TIME(1)-START_USER_TIME(1),", System time = ", &
     & STOP_SYSTEM_TIME(1)-START_SYSTEM_TIME(1),ERR,ERROR,*999)
   
!   this causes issues
!   CALL CMISS_FINALISE(ERR,ERROR,*999)
 
   WRITE(*,'(A)') "Program successfully completed."
  
   STOP
999 CALL CMISS_WRITE_ERROR(ERR,ERROR)
   STOP


END PROGRAM StokesFlow
