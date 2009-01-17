!> \file
!> $Id: cellml.f90 $
!> \author David Nickerson <nickerso@users.sourceforge.net>
!> \brief This module is a openCMISS(cm) buffer module to openCMISS(cellml).
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

!> Doxygen comments get ignored when writing this out in the documentation

! An example application demonstrating the complete openCMISS(cellml) API.
! This example is designed to both provide a guide to the use of CellML models 
! in openCMISS and to test (guide?) the openCMISS(cellml) API implementation.
! In this example we aim to define a regular mesh, define a field ....
PROGRAM CELLMLEXAMPLE

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CMISS
  USE CMISS_MPI
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
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
  USE PROBLEM_CONSTANTS
  USE PROBLEM_ROUTINES
  USE REGION_ROUTINES
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  !Program parameters
  !==================

  ! The mesh size (mm) that will be used.
  REAL(DP), PARAMETER :: HEIGHT=1.0_DP
  REAL(DP), PARAMETER :: WIDTH=2.0_DP
  REAL(DP), PARAMETER :: LENGTH=3.0_DP

  !Program types
  !=============
  
  !Program variables
  !=================

  !Generic CMISS variables
  !=======================
  
  ! The error flag and string
  INTEGER(INTG) :: ERR
  TYPE(VARYING_STRING) :: ERROR

  ! Diagnostic ????
  INTEGER(INTG) :: DIAG_LEVEL_LIST(5)
  CHARACTER(LEN=MAXSTRLEN) :: DIAG_ROUTINE_LIST(1),TIMING_ROUTINE_LIST(1)

  ! Intialise CMISS
  CALL CMISS_INITIALISE(ERR,ERROR,*999)

  ! Finalise CMISS
  CALL CMISS_FINALISE(ERR,ERROR,*999)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
999 CALL CMISS_WRITE_ERROR(ERR,ERROR)
  STOP

END PROGRAM CELLMLEXAMPLE
