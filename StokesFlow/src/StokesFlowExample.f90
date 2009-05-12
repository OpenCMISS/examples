
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

! OpenCMISS Modules

   USE BASE_ROUTINES				! needed
   USE BASIS_ROUTINES				! needed
   USE CMISS					! needed
   USE CMISS_MPI				! needed with MPI
   USE COMP_ENVIRONMENT				! needed with MPI
   USE CONSTANTS				! basic definitions: pi, e, I, ...
   USE CONTROL_LOOP_ROUTINES			! needed
   USE COORDINATE_ROUTINES			! needed
   USE DOMAIN_MAPPINGS				! needed
   USE EQUATIONS_ROUTINES			! needed
   USE EQUATIONS_SET_CONSTANTS			! needed
   USE EQUATIONS_SET_ROUTINES			! needed
   USE FIELD_ROUTINES				! needed
   USE FIELD_IO_ROUTINES			! not needed
   USE INPUT_OUTPUT				! needed
   USE ISO_VARYING_STRING			! needed
   USE KINDS					! needed
   USE MESH_ROUTINES				! needed
   USE MPI					! needed with MPI
   USE PROBLEM_CONSTANTS			! needed
   USE PROBLEM_ROUTINES				! needed
   USE REGION_ROUTINES				! needed
   USE SOLVER_ROUTINES				! needed
   USE TIMER					! needed
   USE TYPES					! needed
!!!!!
#ifdef WIN32
   USE IFQWIN
#endif

! cmHeart input module
  USE IMPORT_CMHEART

IMPLICIT NONE

   !Program types
  TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOF_MAPPING
  TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
  TYPE(MESH_TYPE), POINTER :: MESH
  TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
  TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
  TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
  TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD, DEPENDENT_FIELD
  TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
  TYPE(REGION_TYPE), POINTER :: REGION
  TYPE(SOLVER_TYPE), POINTER :: SOLVER
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
  TYPE(BASIS_TYPE), POINTER :: BASIS_M,BASIS_V,BASIS_P
  TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS_M,MESH_ELEMENTS_P,MESH_ELEMENTS_V
  TYPE(NODES_TYPE), POINTER :: NODES

   !Program variables
  INTEGER(INTG) :: NUMBER_OF_DOMAINS
  INTEGER(INTG) :: MPI_IERROR
  INTEGER(INTG) :: EQUATIONS_SET_INDEX
  LOGICAL :: EXPORT_FIELD,IMPORT_FIELD
  TYPE(VARYING_STRING) :: FILE,METHOD
  REAL(SP) :: START_USER_TIME(1),STOP_USER_TIME(1),START_SYSTEM_TIME(1),STOP_SYSTEM_TIME(1)
  INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES
  INTEGER(INTG) :: MY_COMPUTATIONAL_NODE_NUMBER
  INTEGER(INTG) :: ERR
  TYPE(VARYING_STRING) :: ERROR
  INTEGER(INTG) :: DIAG_LEVEL_LIST(5)
  CHARACTER(LEN=MAXSTRLEN) :: DIAG_ROUTINE_LIST(1),TIMING_ROUTINE_LIST(1)

   !User types
  TYPE(EXPORT_CONTAINER):: CM

   !User variables
  INTEGER:: DECOMPOSITION_USER_NUMBER
  INTEGER:: GEOMETRIC_FIELD_USER_NUMBER
  INTEGER:: DEPENDENT_FIELD_USER_NUMBER
  INTEGER:: DEPENDENT_FIELD_NUMBER_OF_VARIABLES
  INTEGER:: DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
  INTEGER:: REGION_USER_NUMBER
  INTEGER:: COORDINATE_USER_NUMBER
  INTEGER:: MESH_NUMBER_OF_COMPONENTS
  INTEGER:: k,l,m,n
  INTEGER:: X_DIRECTION,Y_DIRECTION,Z_DIRECTION

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
  !Set the window parameters
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Start the creation of a new RC coordinate system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(COORDINATE_SYSTEM)
  COORDINATE_USER_NUMBER=1
  CALL COORDINATE_SYSTEM_CREATE_START(COORDINATE_USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
      !Set the coordinate system dimension
      CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,CM%D,ERR,ERROR,*999)
  CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR,ERROR,*999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Start the creation of a region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(REGION)
  REGION_USER_NUMBER=1
  CALL REGION_CREATE_START(REGION_USER_NUMBER,REGION,ERR,ERROR,*999)
      !Set the regions coordinate system
      CALL REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*999)
  CALL REGION_CREATE_FINISH(REGION,ERR,ERROR,*999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Start the creation of a basis for spatial, velocity and pressure field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(BASIS_M)
  NULLIFY(BASIS_V)
  NULLIFY(BASIS_P)

  !Mesh: Set type to Lagrange/Simplex and define interpolation order
  CALL BASIS_CREATE_START(CM%ID_M,BASIS_M,ERR,ERROR,*999)
      CALL BASIS_TYPE_SET(BASIS_M,CM%IT_T,ERR,ERROR,*999)
      CALL BASIS_NUMBER_OF_XI_SET(BASIS_M,CM%D,ERR,ERROR,*999)

      IF (CM%D==2) THEN
	CALL BASIS_INTERPOLATION_XI_SET(BASIS_M,(/CM%IT_M,CM%IT_M/),ERR,ERROR,*999)
      ELSE IF (CM%D==3) THEN
	CALL BASIS_INTERPOLATION_XI_SET(BASIS_M,(/CM%IT_M,CM%IT_M,CM%IT_M/),ERR,ERROR,*999)
      ELSE
	GOTO 999
      END IF

  CALL BASIS_CREATE_FINISH(BASIS_M,ERR,ERROR,*999)

  !Velocity: Set type to Lagrange/Simplex and define interpolation order
  CALL BASIS_CREATE_START(CM%ID_V,BASIS_V,ERR,ERROR,*999)
      CALL BASIS_TYPE_SET(BASIS_V,CM%IT_T,ERR,ERROR,*999)
      CALL BASIS_NUMBER_OF_XI_SET(BASIS_V,CM%D,ERR,ERROR,*999)

      IF (CM%D==2) THEN
	CALL BASIS_INTERPOLATION_XI_SET(BASIS_V,(/CM%IT_V,CM%IT_V/),ERR,ERROR,*999)
      ELSE IF (CM%D==3) THEN
	CALL BASIS_INTERPOLATION_XI_SET(BASIS_V,(/CM%IT_V,CM%IT_V,CM%IT_V/),ERR,ERROR,*999)
      ELSE
	GOTO 999
      END IF

  CALL BASIS_CREATE_FINISH(BASIS_V,ERR,ERROR,*999)

  !Pressure: Set type to Lagrange/Simplex and define interpolation order
  CALL BASIS_CREATE_START(CM%ID_P,BASIS_P,ERR,ERROR,*999)
      CALL BASIS_TYPE_SET(BASIS_P,CM%IT_T,ERR,ERROR,*999)
      CALL BASIS_NUMBER_OF_XI_SET(BASIS_P,CM%D,ERR,ERROR,*999)

      IF (CM%D==2) THEN
	CALL BASIS_INTERPOLATION_XI_SET(BASIS_P,(/CM%IT_P,CM%IT_P/),ERR,ERROR,*999)
      ELSE IF (CM%D==3) THEN
	CALL BASIS_INTERPOLATION_XI_SET(BASIS_P,(/CM%IT_P,CM%IT_P,CM%IT_P/),ERR,ERROR,*999)
      ELSE
	GOTO 999
      END IF

  CALL BASIS_CREATE_FINISH(BASIS_P,ERR,ERROR,*999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Create a mesh with three mesh components for different field interpolations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(NODES)
  NULLIFY(MESH_ELEMENTS_M)
  NULLIFY(MESH_ELEMENTS_P)
  NULLIFY(MESH_ELEMENTS_V)

  MESH_NUMBER_OF_COMPONENTS=3

  ! Define nodes from user input
  CALL NODES_CREATE_START(CM%N_T,REGION,NODES,ERR,ERROR,*999)
  CALL NODES_CREATE_FINISH(REGION,ERR,ERROR,*999)



  ! Define elements from user input
  CALL MESH_CREATE_START(1,REGION,CM%D,MESH,ERR,ERROR,*999)
      CALL MESH_NUMBER_OF_ELEMENTS_SET(MESH,CM%E_T,ERR,ERROR,*999)
      CALL MESH_NUMBER_OF_COMPONENTS_SET(MESH,MESH_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)


      !Mesh:
      CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,CM%ID_M,BASIS_M,MESH_ELEMENTS_M,ERR,ERROR,*999)
	  DO k=1,CM%E_T
	    CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(k,MESH_ELEMENTS_M, &
	    CM%M(k,1:CM%EN_M),ERR,ERROR,*999)


!	WRITE(*,*)'k',k,'CM%M(k,1:CM%EN_M)',CM%M(k,1:CM%EN_M)
!	WRITE(*,*)'k',k,'CM%M(k,1:CM%EN_M)',CM%M(k,28)

	  END DO
      CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH,CM%ID_M,ERR,ERROR,*999)



      !Velocity:
      CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,CM%ID_V,BASIS_V,MESH_ELEMENTS_V,ERR,ERROR,*999)
	  DO k=1,CM%E_T
	    CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(k,MESH_ELEMENTS_V, &
	    CM%V(k,1:CM%EN_V),ERR,ERROR,*999)

!	WRITE(*,*)'k',k,'CM%V(k,1:CM%EN_V)',CM%V(k,1:CM%EN_V)

	  END DO
      CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH,CM%ID_V,ERR,ERROR,*999)



      !Pressure:
      CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,CM%ID_P,BASIS_P,MESH_ELEMENTS_P,ERR,ERROR,*999)
	  DO k=1,CM%E_T
	    CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(k,MESH_ELEMENTS_P, &
	    CM%P(k,1:CM%EN_P),ERR,ERROR,*999)

!	WRITE(*,*)'k',k,'CM%P(k,1:CM%EN_P)',CM%P(k,1:CM%EN_P)


	  END DO
      CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH,CM%ID_P,ERR,ERROR,*999)

  
  CALL MESH_CREATE_FINISH(REGION,MESH,ERR,ERROR,*999)
  





WRITE(*,*)'5'

  !CALL MESH_CREATE_FINISH(REGION,MESH,ERR,ERROR,*999)

WRITE(*,*)'6'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Create a decomposition for mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(DECOMPOSITION)

  DECOMPOSITION_USER_NUMBER=1

  CALL DECOMPOSITION_CREATE_START(DECOMPOSITION_USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*999)

      !Set the decomposition to be a general decomposition with the specified number of domains
      CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)
      CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)

  CALL DECOMPOSITION_CREATE_FINISH(MESH,DECOMPOSITION,ERR,ERROR,*999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Define geometric field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(GEOMETRIC_FIELD)

  X_DIRECTION=1
  Y_DIRECTION=2
  Z_DIRECTION=3
  GEOMETRIC_FIELD_USER_NUMBER=1

  !Create geometric field
  CALL FIELD_CREATE_START(GEOMETRIC_FIELD_USER_NUMBER,REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)
      CALL FIELD_TYPE_SET(GEOMETRIC_FIELD_USER_NUMBER,REGION,FIELD_GEOMETRIC_TYPE,ERR,ERROR,*999)
      CALL FIELD_MESH_DECOMPOSITION_SET(GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*999)

      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,X_DIRECTION,CM%ID_M,ERR,ERROR,*999)
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,Y_DIRECTION,CM%ID_M,ERR,ERROR,*999)
      IF(CM%D==3) THEN
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,Z_DIRECTION,CM%ID_M,ERR,ERROR,*999)
      ENDIF

  CALL FIELD_CREATE_FINISH(REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)

  !Set geometric field parameters and do update
  DO k=1,CM%N_M	!number of mesh nodes
    DO j=1,CM%D !dimensions
      CALL FIELD_PARAMETER_SET_UPDATE_NODE(GEOMETRIC_FIELD,FIELD_VALUES_SET_TYPE,CM%ID_M,k,j, &
      FIELD_U_VARIABLE_TYPE,CM%N(k,j),ERR,ERROR,*999)
    END DO
  END DO


  CALL FIELD_PARAMETER_SET_UPDATE_START(GEOMETRIC_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(GEOMETRIC_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD)) GEOMETRIC_FIELD=>REGION%FIELDS%FIELDS(1)%PTR



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Create equations set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NULLIFY(EQUATIONS_SET)

  !Set the equations set to be a Stokes Flow problem
  CALL EQUATIONS_SET_CREATE_START(1,REGION,GEOMETRIC_FIELD,EQUATIONS_SET,ERR,ERROR,*999)
    CALL EQUATIONS_SET_SPECIFICATION_SET(EQUATIONS_SET,EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_STOKES_FLUID_TYPE, &
    & EQUATIONS_SET_STANDARD_STOKES_SUBTYPE,ERR,ERROR,*999)
  CALL EQUATIONS_SET_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)

  !Create the equations set dependent field variables
  CALL EQUATIONS_SET_DEPENDENT_CREATE_START(EQUATIONS_SET,ERR,ERROR,*999)
  CALL EQUATIONS_SET_DEPENDENT_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
  DEPENDENT_DOF_MAPPING=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Define boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create the problem fixed conditions
  CALL EQUATIONS_SET_FIXED_CONDITIONS_CREATE_START(EQUATIONS_SET,ERR,ERROR,*999)
    !Set bc's
    CALL EQUATIONS_SET_FIXED_CONDITIONS_SET_DOF(EQUATIONS_SET,1,EQUATIONS_SET_FIXED_BOUNDARY_CONDITION, &
    & 42.0_DP,ERR,ERROR,*999)
  CALL EQUATIONS_SET_FIXED_CONDITIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Define equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  NULLIFY(EQUATIONS)

  CALL EQUATIONS_SET_EQUATIONS_CREATE_START(EQUATIONS_SET,ERR,ERROR,*999)
    !Get the equations
    CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
    !Set the equations matrices sparsity type
    CALL EQUATIONS_SPARSITY_TYPE_SET(EQUATIONS,EQUATIONS_SPARSE_MATRICES,ERR,ERROR,*999)
    CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_ELEMENT_MATRIX_OUTPUT,ERR,ERROR,*999)
  CALL EQUATIONS_SET_EQUATIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Define problem and solver settings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  NULLIFY(PROBLEM)
  !Set the problem to be a standard Stokes problem
  CALL PROBLEM_CREATE_START(1,PROBLEM,ERR,ERROR,*999)
    CALL PROBLEM_SPECIFICATION_SET(PROBLEM,PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_STOKES_FLUID_TYPE, &
    & PROBLEM_STANDARD_STOKES_SUBTYPE,ERR,ERROR,*999)
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

  NULLIFY(SOLVER)
  NULLIFY(SOLVER_EQUATIONS)
  !Create the problem solver equations
  CALL PROBLEM_SOLVER_EQUATIONS_CREATE_START(PROBLEM,ERR,ERROR,*999)
    CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,1,SOLVER,ERR,ERROR,*999)
    CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
    CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
    !Add in the equations set
    CALL PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD(PROBLEM,CONTROL_LOOP_NODE,1,EQUATIONS_SET,EQUATIONS_SET_INDEX,ERR,ERROR,*999)
  CALL PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! !Solve the problem
! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! 
! ! ! WRITE(*,*)'42'
! ! !   CALL PROBLEM_SOLVE(PROBLEM,ERR,ERROR,*999)
! ! ! WRITE(*,*)'79'
! ! ! 
! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! !AFTERBURNER
! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! 
! ! !    EXPORT_FIELD=.FALSE.
! ! !    METHOD="FORTRAN"
! ! !    IF(EXPORT_FIELD) THEN
! ! !      FILE="InputCMHeartMeshExample"
! ! !      CALL FIELD_IO_NODES_EXPORT(REGION%FIELDS, FILE, METHOD, ERR,ERROR,*999)
! ! !      CALL FIELD_IO_ELEMENTS_EXPORT(REGION%FIELDS, FILE, METHOD, ERR,ERROR,*999)
! ! !    ENDIF
! ! ! 
! ! ! 
! ! !    !Calculate the stop times and write out the elapsed user and system times
! ! !    CALL CPU_TIMER(USER_CPU,STOP_USER_TIME,ERR,ERROR,*999)
! ! !    CALL CPU_TIMER(SYSTEM_CPU,STOP_SYSTEM_TIME,ERR,ERROR,*999)
! ! ! 
! ! !    CALL WRITE_STRING_TWO_VALUE(GENERAL_OUTPUT_TYPE,"User time = ",STOP_USER_TIME(1)-START_USER_TIME(1),", System time = ", &
! ! !      & STOP_SYSTEM_TIME(1)-START_SYSTEM_TIME(1),ERR,ERROR,*999)
! ! ! 
! ! ! !   this causes issues
! ! ! !   CALL CMISS_FINALISE(ERR,ERROR,*999)

   WRITE(*,'(A)') "Program successfully completed."

   STOP
999 CALL CMISS_WRITE_ERROR(ERR,ERROR)
   STOP


END PROGRAM StokesFlow
