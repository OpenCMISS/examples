!> \file
!> $Id: InputCMHeartMesh.f90 20 2009-04-08 20:22:52Z cpb $
!> \author Sebastian Krittian
!> \brief This is an example program which imports CMHeart input files.
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

!> \example InputCMHeartMesh/src/InputCMHeartMeshExample.f90
!! Example program which imports CMHeart input files.
!<

!> Main program

PROGRAM MeshFromCMHeart

! OpenCMISS Modules

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
 
! cmHeart input module
  USE IMPORT_CMHEART

IMPLICIT NONE
  

   !Program types
  TYPE(BASIS_TYPE), POINTER :: BASIS
  TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
  TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
  TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
  TYPE(MESH_TYPE), POINTER :: MESH
  TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
  TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
  TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
  TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD,DEPENDENT_FIELD
  TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
  TYPE(REGION_TYPE), POINTER :: REGION,WORLD_REGION
  TYPE(SOLVER_TYPE), POINTER :: SOLVER
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
  TYPE(VARYING_STRING) :: FILE,METHOD

   !Program variables
  INTEGER(INTG) :: NUMBER_OF_DOMAINS
  INTEGER(INTG) :: MPI_IERROR
  INTEGER(INTG) :: EQUATIONS_SET_INDEX
  LOGICAL :: EXPORT_FIELD,IMPORT_FIELD,LAPLACE
  REAL(SP) :: START_USER_TIME(1),STOP_USER_TIME(1),START_SYSTEM_TIME(1),STOP_SYSTEM_TIME(1)
  INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES
  INTEGER(INTG) :: MY_COMPUTATIONAL_NODE_NUMBER
  TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS,MESH_ELEMENTS2
  TYPE(NODES_TYPE), POINTER :: NODES
  INTEGER(INTG) :: ERR
  TYPE(VARYING_STRING) :: ERROR
  INTEGER(INTG) :: DIAG_LEVEL_LIST(5)
  CHARACTER(LEN=MAXSTRLEN) :: DIAG_ROUTINE_LIST(1),TIMING_ROUTINE_LIST(1)

   !User types
  TYPE(EXPORT_CONTAINER):: CM

   !User variables 
  INTEGER:: DECOMPOSITION_USER_NUMBER
  INTEGER:: COORDINATE_USER_NUMBER 
  INTEGER:: REGION_USER_NUMBER
  INTEGER:: k,l,m,n


#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program starts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parametersDEPENDENT_FIELD%DECOMPOSITION
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Import cmHeart Information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   CALL READ_CMHEART_EXE
   CALL RECV_CMHEART_EXE(CM)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Intialise cmiss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   NULLIFY(WORLD_REGION)
   CALL CMISS_INITIALISE(WORLD_REGION,ERR,ERROR,*999)
  
!   !Set all diganostic levels on for testing
!   DIAG_LEVEL_LIST(1)=1
!   DIAG_LEVEL_LIST(2)=2
!   DIAG_LEVEL_LIST(3)=3
!   DIAG_LEVEL_LIST(4)=4
!   DIAG_LEVEL_LIST(5)=5
!   !DIAG_ROUTINE_LIST(1)=""
!   !CALL DIAGNOSTICS_SET_ON(ALL_DIAG_TYPE,DIAG_LEVEL_LIST,"MoreComplexMeshExample",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
!   CALL DIAGNOSTICS_SET_ON(ALL_DIAG_TYPE,DIAG_LEVEL_LIST,"InputCMHeartMeshExample",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
 
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Start the creation of a new RC coordinate system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(COORDINATE_SYSTEM)
  COORDINATE_USER_NUMBER=1

  CALL COORDINATE_SYSTEM_CREATE_START(COORDINATE_USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
    !Set the coordinate system to be DIMENSIONS
    CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,CM%D,ERR,ERROR,*999)
  CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR,ERROR,*999)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Start the creation of a region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(REGION)
  REGION_USER_NUMBER=1

  CALL REGION_CREATE_START(REGION_USER_NUMBER,WORLD_REGION,REGION,ERR,ERROR,*999)
    !Set the regions coordinate system to the RC coordinate system that we have created
    CALL REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*999)
  CALL REGION_CREATE_FINISH(REGION,ERR,ERROR,*999)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Start the creation of a basis for spatial field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(BASIS)


  CALL BASIS_CREATE_START(CM%ID_M,BASIS,ERR,ERROR,*999)  
      CALL BASIS_TYPE_SET(BASIS,CM%IT_T,ERR,ERROR,*999)  
      CALL BASIS_NUMBER_OF_XI_SET(BASIS,CM%D,ERR,ERROR,*999)

      IF (CM%D==2) THEN
	CALL BASIS_INTERPOLATION_XI_SET(BASIS,(/CM%IT_M,CM%IT_M/),ERR,ERROR,*999)
      ELSE IF (CM%D==3) THEN
	CALL BASIS_INTERPOLATION_XI_SET(BASIS,(/CM%IT_M,CM%IT_M,CM%IT_M/),ERR,ERROR,*999)
      ELSE
	GOTO 999
      END IF

  CALL BASIS_CREATE_FINISH(BASIS,ERR,ERROR,*999)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Create a mesh with one mesh component
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  NULLIFY(NODES)
  NULLIFY(MESH)
  NULLIFY(MESH_ELEMENTS)

  CALL NODES_CREATE_START(REGION,CM%N_T,NODES,ERR,ERROR,*999)
  CALL NODES_CREATE_FINISH(NODES,ERR,ERROR,*999)

  CALL MESH_CREATE_START(CM%ID_M,REGION,CM%D,MESH,ERR,ERROR,*999)
  CALL MESH_NUMBER_OF_ELEMENTS_SET(MESH,CM%E_T,ERR,ERROR,*999)

    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,1,BASIS,MESH_ELEMENTS,ERR,ERROR,*999)
      DO k=1,CM%E_M
	CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(k,MESH_ELEMENTS, &
	  CM%M(k,1:CM%EN_M),ERR,ERROR,*999)
      END DO
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH,1,ERR,ERROR,*999)
  
  CALL MESH_CREATE_FINISH(MESH,ERR,ERROR,*999)
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Create a decomposition for mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(DECOMPOSITION)
  DECOMPOSITION_USER_NUMBER=1

  CALL DECOMPOSITION_CREATE_START(DECOMPOSITION_USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*999)
      CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)
      CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)
  CALL DECOMPOSITION_CREATE_FINISH(MESH,DECOMPOSITION,ERR,ERROR,*999)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Define geometric field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(GEOMETRIC_FIELD)
  REGION_USER_NUMBER=1

  CALL FIELD_CREATE_START(REGION_USER_NUMBER,REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)
      CALL FIELD_MESH_DECOMPOSITION_SET(GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*999)
      !Set the domain to be used by the field components
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1,1,ERR,ERROR,*999)
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,2,1,ERR,ERROR,*999)
    
      IF(CM%D==3) THEN
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,3,1,ERR,ERROR,*999)
      ENDIF

  CALL FIELD_CREATE_FINISH(GEOMETRIC_FIELD,ERR,ERROR,*999)

  !Set the geometric field values
  DO k=1,CM%N_M
    DO j=1,CM%D
      CALL FIELD_PARAMETER_SET_UPDATE_NODE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,k,j, &
        & CM%N(k,j),ERR,ERROR,*999)
    END DO
!      WRITE(*,*)'k',k,'CM%N(k,1:CM%D)',CM%N(k,1:CM%D)
  END DO

  CALL FIELD_PARAMETER_SET_UPDATE_START(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

  LAPLACE=.FALSE.
  IF(LAPLACE) THEN
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Create equations set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  NULLIFY(EQUATIONS_SET)

  CALL EQUATIONS_SET_CREATE_START(1,REGION,GEOMETRIC_FIELD,EQUATIONS_SET,ERR,ERROR,*999)
  !Set the equations set to be a standard Laplace problem
    CALL EQUATIONS_SET_SPECIFICATION_SET(EQUATIONS_SET,EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_LAPLACE_EQUATION_TYPE, &
    & EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE,ERR,ERROR,*999)
  CALL EQUATIONS_SET_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)

  !Create the equations set dependent field variables
  CALL EQUATIONS_SET_DEPENDENT_CREATE_START(EQUATIONS_SET,2,DEPENDENT_FIELD,ERR,ERROR,*999)
  CALL EQUATIONS_SET_DEPENDENT_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Define equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(EQUATIONS)

  !Create the equations set equations
  CALL EQUATIONS_SET_EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
    !Set the equations matrices sparsity type
    CALL EQUATIONS_SPARSITY_TYPE_SET(EQUATIONS,EQUATIONS_SPARSE_MATRICES,ERR,ERROR,*999)
     CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_ELEMENT_MATRIX_OUTPUT,ERR,ERROR,*999)
  CALL EQUATIONS_SET_EQUATIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Define boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create the equations set boundary conditions
  NULLIFY(BOUNDARY_CONDITIONS)
  CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
  !Set bc's
  CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD_U_VARIABLE_TYPE,1,BOUNDARY_CONDITION_FIXED, &
      & 42.0_DP,ERR,ERROR,*999)
  CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Define problem and solver settings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(PROBLEM)
  CALL PROBLEM_CREATE_START(1,PROBLEM,ERR,ERROR,*999)
  !Set the problem to be a standard Laplace problem
    CALL PROBLEM_SPECIFICATION_SET(PROBLEM,PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_LAPLACE_EQUATION_TYPE, &
    & PROBLEM_STANDARD_LAPLACE_SUBTYPE,ERR,ERROR,*999)
  CALL PROBLEM_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !Create the problem control loop
  CALL PROBLEM_CONTROL_LOOP_CREATE_START(PROBLEM,ERR,ERROR,*999)
  CALL PROBLEM_CONTROL_LOOP_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !Start the creation of the problem solvers
  NULLIFY(SOLVER)
  CALL PROBLEM_SOLVERS_CREATE_START(PROBLEM,ERR,ERROR,*999)
    CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,1,SOLVER,ERR,ERROR,*999)
    CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_MATRIX_OUTPUT,ERR,ERROR,*999)
  CALL PROBLEM_SOLVERS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !Create the problem solver equations
  NULLIFY(SOLVER)
  NULLIFY(SOLVER_EQUATIONS)
  CALL PROBLEM_SOLVER_EQUATIONS_CREATE_START(PROBLEM,ERR,ERROR,*999)
    CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,1,SOLVER,ERR,ERROR,*999)
    CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
    CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
    CALL SOLVER_EQUATIONS_EQUATIONS_SET_ADD(SOLVER_EQUATIONS,EQUATIONS_SET,EQUATIONS_SET_INDEX,ERR,ERROR,*999)
  !Finish the problem solver equations
  CALL PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Solve the problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL PROBLEM_SOLVE(PROBLEM,ERR,ERROR,*999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!AFTERBURNER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ELSE
  WRITE(*,*)'Export mesh only'
END IF

   EXPORT_FIELD=.TRUE.
   METHOD="FORTRAN"
   IF(EXPORT_FIELD) THEN
     FILE="InputCMHeartMeshExample"
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


END PROGRAM MeshFromCMHeart
