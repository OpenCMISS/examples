!> \file
!> \author David Nickerson <nickerso@users.sourceforge.net>
!> \brief This example illustrates the definition of a field on the created mesh and some simple manipulations of that field. In this example we directly access the field data rather than using the field API.
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

!> Doxygen comments ignored when generating example pages

!> \example simple-field-manipulation-direct-access/src/simple-field-manipulation-direct-accessExample.f90
!! Example illustrating the creation of a regular mesh and then defining and manipulating a field over that mesh.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/simple-field-manipulation-direct-access/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/simple-field-manipulation-direct-access/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM SimpleFieldManipulationExample
  
#ifndef NOMPIMOD
  USE MPI
#endif

  USE OpenCMISS
  USE OpenCMISS_Iron

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Test program parameters

  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=3.0_CMISSRP

  !Program variables
  
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS

  INTEGER(CMISSIntg), PARAMETER :: field_geometry_id=1
  INTEGER(CMISSIntg), PARAMETER :: field_geometry_number_of_variables=1
  INTEGER(CMISSIntg) :: field_geometry_number_of_components

  INTEGER(CMISSIntg), PARAMETER :: field_general_id=field_geometry_id+1
  INTEGER(CMISSIntg), PARAMETER :: field_general_number_of_variables=1
  INTEGER(CMISSIntg), PARAMETER :: field_general_number_of_components=3
  
  INTEGER(CMISSIntg) :: NUMBER_COMPUTATIONAL_NODES
  INTEGER(CMISSIntg) :: MY_COMPUTATIONAL_NODE_NUMBER
  INTEGER(CMISSIntg) :: MPI_IERROR

  TYPE(BASIS_TYPE), POINTER :: BASIS
  TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
  TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
  TYPE(MESH_TYPE), POINTER :: MESH
  TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
  TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD,GENERAL_FIELD
  TYPE(REGION_TYPE), POINTER :: REGION,WORLD_REGION

  TYPE(DOMAIN_TYPE), POINTER :: GENERAL_FIELD_DOMAIN
  TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
  REAL(CMISSRP) :: field_value

  LOGICAL :: EXPORT_FIELD
  TYPE(VARYING_STRING) :: FILE,METHOD

  REAL(SP) :: START_USER_TIME(1),STOP_USER_TIME(1),START_SYSTEM_TIME(1),STOP_SYSTEM_TIME(1)

  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: ERR
  TYPE(VARYING_STRING) :: ERROR
  
  INTEGER(CMISSIntg) :: DIAG_LEVEL_LIST(5)
  CHARACTER(LEN=MAXSTRLEN) :: DIAG_ROUTINE_LIST(1)

  INTEGER(CMISSIntg) :: nc,node_idx,total_number_of_nodes,global_node_number

  !Intialise cmiss
  NULLIFY(WORLD_REGION)
  CALL cmfe_Initialise(WORLD_REGION,ERR)
  
  !Set all diganostic levels on for testing
  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5
  !CALL DIAGNOSTICS_SET_ON(ALL_DIAG_TYPE,DIAG_LEVEL_LIST,"",DIAG_ROUTINE_LIST,ERR)
  
  !Calculate the start times
  CALL CPU_TIMER(USER_CPU,START_USER_TIME,ERR)
  CALL CPU_TIMER(SYSTEM_CPU,START_SYSTEM_TIME,ERR)

  !Get the number of computational nodes
  NUMBER_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
  IF(ERR/=0) GOTO 999
  !Get my computational node number
  MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
  IF(ERR/=0) GOTO 999
  
  NUMBER_GLOBAL_X_ELEMENTS=2
  NUMBER_GLOBAL_Y_ELEMENTS=2
  NUMBER_GLOBAL_Z_ELEMENTS=0
  NUMBER_OF_DOMAINS=2
  
  !Read in the number of elements in the X & Y directions, and the number of partitions on the master node (number 0)
  IF(MY_COMPUTATIONAL_NODE_NUMBER==0) THEN
    WRITE(*,'("Enter the number of elements in the X direction (2):")')
    READ(*,*) NUMBER_GLOBAL_X_ELEMENTS
    WRITE(*,'("Enter the number of elements in the Y direction (2):")')
    READ(*,*) NUMBER_GLOBAL_Y_ELEMENTS
    WRITE(*,'("Enter the number of elements in the Z direction (0):")')
    READ(*,*) NUMBER_GLOBAL_Z_ELEMENTS
    WRITE(*,'("Enter the number of domains (2):")')
    READ(*,*) NUMBER_OF_DOMAINS
  ENDIF
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR)
  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"COMPUTATIONAL ENVIRONMENT:",ERR)
  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Total number of computaional nodes = ",NUMBER_COMPUTATIONAL_NODES, &
    & ERR)
  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  My computational node number = ",MY_COMPUTATIONAL_NODE_NUMBER,ERR)

  !Start the creation of a new RC coordinate system
  NULLIFY(COORDINATE_SYSTEM)
  CALL COORDINATE_SYSTEM_CREATE_START(1,COORDINATE_SYSTEM,ERR)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,2,ERR)
  ELSE
    !Set the coordinate system to be 3D
    CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,3,ERR)
  ENDIF
  !Finish the creation of the coordinate system
  CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR)

  !Start the creation of the region
  NULLIFY(REGION)
  CALL REGION_CREATE_START(1,WORLD_REGION,REGION,ERR)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR)
  !Finish the creation of the region
  CALL REGION_CREATE_FINISH(REGION,ERR)

  !Start the creation of a basis (default is trilinear lagrange)
  NULLIFY(BASIS)
  CALL BASIS_CREATE_START(1,BASIS,ERR)  
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
     !Set the basis to be a bilinear Lagrange basis
     CALL BASIS_NUMBER_OF_XI_SET(BASIS,2,ERR)
  ELSE
     !Set the basis to be a trilinear Lagrange basis
     CALL BASIS_NUMBER_OF_XI_SET(BASIS,3,ERR)
  ENDIF
  !Finish the creation of the basis
  CALL BASIS_CREATE_FINISH(BASIS,ERR)

  !Start the creation of a generated mesh in the region
  NULLIFY(GENERATED_MESH)
  CALL GENERATED_MESH_CREATE_START(1,REGION,GENERATED_MESH,ERR)
  !Set up a regular mesh
  CALL GENERATED_MESH_TYPE_SET(GENERATED_MESH,1, &
       & ERR)
  CALL GENERATED_MESH_BASIS_SET(GENERATED_MESH,BASIS,ERR)

  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
     CALL GENERATED_MESH_EXTENT_SET(GENERATED_MESH,(/WIDTH,HEIGHT/),ERR)
     CALL GENERATED_MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/), &
          & ERR)
  ELSE
     CALL GENERATED_MESH_EXTENT_SET(GENERATED_MESH,(/WIDTH,HEIGHT,LENGTH/),ERR)
     CALL GENERATED_MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH,(/NUMBER_GLOBAL_X_ELEMENTS, &
          & NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS/), ERR)
  ENDIF

  !Finish the creation of a generated mesh in the region
  CALL GENERATED_MESH_CREATE_FINISH(GENERATED_MESH,1,MESH,ERR) 

  !Create a decomposition
  NULLIFY(DECOMPOSITION)
  CALL DECOMPOSITION_CREATE_START(1,MESH,DECOMPOSITION,ERR)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR)
  CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR)
  CALL DECOMPOSITION_CREATE_FINISH(DECOMPOSITION,ERR)

  !Start to create a default (geometric) field on the region
  NULLIFY(GEOMETRIC_FIELD)
  CALL FIELD_CREATE_START(field_geometry_id,REGION,GEOMETRIC_FIELD,ERR)
  !Set the decomposition to use
  CALL FIELD_MESH_DECOMPOSITION_SET(GEOMETRIC_FIELD,DECOMPOSITION,ERR)
  !Set the domain to be used by the field components
  !NB these are needed now as the default mesh component number is 1
  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1,1,ERR)
  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,2,1,ERR)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
     CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,3,1,ERR)
  ENDIF
  !Finish creating the field
  CALL FIELD_CREATE_FINISH(GEOMETRIC_FIELD,ERR)

  !Update the geometric field parameters
  CALL GENERATED_MESH_GEOMETRIC_PARAMETERS_CALCULATE(GEOMETRIC_FIELD,GENERATED_MESH,ERR)

  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD)) GEOMETRIC_FIELD=>REGION%FIELDS%FIELDS(field_geometry_id)%PTR
  
  FILE="SimpleFieldManipulationExample-geometry"
  METHOD="FORTRAN"
  
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL FIELD_IO_NODES_EXPORT(REGION%FIELDS, FILE, METHOD, ERR)  
    CALL FIELD_IO_ELEMENTS_EXPORT(REGION%FIELDS, FILE, METHOD, ERR)
  ENDIF
  
  ! define the general field
  NULLIFY(GENERAL_FIELD)
  CALL FIELD_CREATE_START(field_general_id,REGION,GENERAL_FIELD,ERR)
  CALL FIELD_TYPE_SET(GENERAL_FIELD,FIELD_GENERAL_TYPE,ERR)  
  CALL FIELD_MESH_DECOMPOSITION_SET(GENERAL_FIELD,DECOMPOSITION,ERR)
  CALL FIELD_GEOMETRIC_FIELD_SET(GENERAL_FIELD,GEOMETRIC_FIELD,ERR)
  CALL FIELD_NUMBER_OF_VARIABLES_SET(GENERAL_FIELD,field_general_number_of_variables,ERR)
  CALL FIELD_NUMBER_OF_COMPONENTS_SET(GENERAL_FIELD,FIELD_U_VARIABLE_TYPE,field_general_number_of_components,ERR)
  DO nc=1,field_general_number_of_components
     CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GENERAL_FIELD,FIELD_U_VARIABLE_TYPE,nc,1,ERR)
  END DO
  !Finish creating the field
  CALL FIELD_CREATE_FINISH(GENERAL_FIELD,ERR)

  FILE="SimpleFieldManipulationExample-general"
  IF(EXPORT_FIELD) THEN
    CALL FIELD_IO_NODES_EXPORT(REGION%FIELDS, FILE, METHOD, ERR)  
    CALL FIELD_IO_ELEMENTS_EXPORT(REGION%FIELDS, FILE, METHOD, ERR)
  ENDIF

  ! update the general field values for each node
  ! get the field's domain for mesh component 1 (there is only 1?)
  GENERAL_FIELD_DOMAIN => GENERAL_FIELD%DECOMPOSITION%DOMAIN(1)%PTR
  DOMAIN_NODES => GENERAL_FIELD_DOMAIN%TOPOLOGY%NODES
  DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
     global_node_number = DOMAIN_NODES%NODES(node_idx)%GLOBAL_NUMBER
     WRITE(*,*) 'Computational node number: ',MY_COMPUTATIONAL_NODE_NUMBER, &
          & '; global node number: ',global_node_number
     field_value = global_node_number
     DO nc=1,field_general_number_of_components
        field_value = field_value + (1.0_CMISSRP * (MY_COMPUTATIONAL_NODE_NUMBER+1))
        CALL FIELD_PARAMETER_SET_UPDATE_NODE(GENERAL_FIELD,FIELD_U_VARIABLE_TYPE, &
             & FIELD_VALUES_SET_TYPE,1,node_idx,nc,field_value, &
             & ERR)
     END DO
  END DO
  
  ! synchronise all ranks
  CALL FIELD_PARAMETER_SET_UPDATE_START(GENERAL_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR)
  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(GENERAL_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR)
  
  FILE="SimpleFieldManipulationExample-general-update1"
  IF(EXPORT_FIELD) THEN
    CALL FIELD_IO_NODES_EXPORT(REGION%FIELDS, FILE, METHOD, ERR)  
    CALL FIELD_IO_ELEMENTS_EXPORT(REGION%FIELDS, FILE, METHOD, ERR)
  ENDIF

  !Output timing summary
  !Calculate the stop times and write out the elapsed user and system times
  CALL CPU_TIMER(USER_CPU,STOP_USER_TIME,ERR)
  CALL CPU_TIMER(SYSTEM_CPU,STOP_SYSTEM_TIME,ERR)

  CALL WRITE_STRING_TWO_VALUE(GENERAL_OUTPUT_TYPE,"User time = ",STOP_USER_TIME(1)-START_USER_TIME(1),", System time = ", &
    & STOP_SYSTEM_TIME(1)-START_SYSTEM_TIME(1),ERR)
  
  ! Finalise cmiss
  CALL CMFE_FINALISE(ERR)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
999 CALL CMFE_WRITE_ERROR(ERR,ERROR)
  STOP

END PROGRAM SimpleFieldManipulationExample

