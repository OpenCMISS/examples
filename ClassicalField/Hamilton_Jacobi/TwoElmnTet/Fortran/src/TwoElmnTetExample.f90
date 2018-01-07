!> \file
!> \author Ali PASHAEI
!> \brief This is an example program to solve a Hamilton-Jacobi (Eikonal form) equation using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: 0.1
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
!> The Initial Developer of the Original Code is  University of Pompeu Fabra,
!> Barcelona, Spain, University of Auckland, Auckland, New Zealand and 
!> University of Oxford, Oxford, United Kingdom. Portions created by the 
!> University of Auckland and University of Oxford are Copyright (C) 2007 
!> by the University of Pompeu Fabra, the University of Auckland and
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

!> \example ClassicalField/Hamilton_Jacobi/TwoElmnTet/src/TwoElmnTetExample.f90
!! Example program to solve a fast electrophysiology equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/TwoElementTriLinear/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/TwoElementTriLinear/build-gnu'>Linux GNU Build</a>
!<

!> Main program

!...................................................................................
PROGRAM TwoElmnTet
!...................................................................................

! OpenCMISS Modules

  USE OpenCMISS
  USE OpenCMISS_Iron

#ifdef WIN32
  USE IFQWIN
#endif

!...................................................................................
  IMPLICIT NONE
!...................................................................................

  !Program variables

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2

  !local variables

  INTEGER(CMISSIntg) :: I,J,K,ELEMENT_NUMBER,MIN_TRIAL_NODE,TRIAL_STATUS,CHANGED_STATUS,SEED_POINT,TOTAL_NUMBER_OF_CONNECTIVITY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_PER_ELEMENT,TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_SNAPSHOTS
  REAL(CMISSRP) :: MIN_VALUE,MAX_VALUE,TIME_STEP
  REAL(CMISSRP) :: SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)
  CHARACTER (LEN=100) :: INPUT_FILE_NAME
  CHARACTER (LEN=100) :: OUTPUT_FILE_NAME
  CHARACTER (LEN=300) :: OUTPUT_FILE_FIELD_TITLE
  CHARACTER (LEN=10)  :: INPUT_FILE_FORMAT
  CHARACTER (LEN=10)  :: OUTPUT_FILE_FORMAT
  CHARACTER (LEN=12)  :: MATERIAL_BEHAVIOUR
  CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_SPEED_FUNCTION
  CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_SEED_VALUE
  CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_CONDUCTIVITY
  CHARACTER (LEN=100) :: STRING

  INTEGER(CMISSIntg) :: Err

  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_RegionType) :: WorldRegion
  
  REAL(CMISSRP), ALLOCATABLE, DIMENSION(:,:):: NODE_LIST
  REAL(CMISSRP), ALLOCATABLE, DIMENSION(:,:):: SPEED_FUNCTION_TABLE
  REAL(CMISSRP), ALLOCATABLE, DIMENSION(:,:):: CONDUCTIVITY_TENSOR
  REAL(CMISSRP), ALLOCATABLE, DIMENSION(:,:):: SPEED_FUNCTION_TABLE_ON_CONNECTIVITY
  REAL(CMISSRP), ALLOCATABLE, DIMENSION(:,:):: CONDUCTIVITY_TENSOR_ON_CONNECTIVITY
  REAL(CMISSRP), ALLOCATABLE, DIMENSION(:)  :: SEED_VALUE
  REAL(CMISSRP), ALLOCATABLE, DIMENSION(:)  :: SNAPSHOT_SEED_VALUE
  REAL(CMISSRP), ALLOCATABLE, DIMENSION(:)  :: RELATIVE_ERROR
  REAL(CMISSRP), ALLOCATABLE, DIMENSION(:)  :: SEED_VALUE_EXACT_SOLUTION
  REAL(CMISSRP), ALLOCATABLE, DIMENSION(:)  :: TRACE_NODES
  
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:,:):: CONNECTIVITY_LIST
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:,:):: ELEMENT_LIST
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:)  :: CONNECTIVITY_NUMBER
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:)  :: TRACE_NODE

  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:)  :: COLUMN_INDEX
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:)  :: TEMP_ARRAY
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:)  :: RAW_INDEX

  CHARACTER (LEN=10), ALLOCATABLE, DIMENSION(:)  :: STATUS_MASK

#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

! INPUT_FILE_FORMAT :   "TABC"   -> tabc format
!			"TETGEN" -> tetgen format
  INPUT_FILE_FORMAT  = "TETGEN"

! INPUT_FILE_NAME :   in order to calculate TOTAL_NUMBER_OF_NODES
  INPUT_FILE_NAME = "src/TwoElmnTet"

! TOTAL_NUMBER_OF_NODES : to set input parameters
  CALL CMFE_NUMBER_OF_INPUT_NODES(INPUT_FILE_NAME,INPUT_FILE_FORMAT,TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,&
  & TOTAL_NUMBER_OF_CONNECTIVITY,Err)
!  Allocate parameters
  ALLOCATE(NODE_LIST(TOTAL_NUMBER_OF_NODES,3),STAT=ERR)
  ALLOCATE(SPEED_FUNCTION_TABLE(TOTAL_NUMBER_OF_NODES,3),STAT=ERR)
  ALLOCATE(CONDUCTIVITY_TENSOR(TOTAL_NUMBER_OF_NODES,9),STAT=ERR)
  ALLOCATE(CONNECTIVITY_LIST(TOTAL_NUMBER_OF_NODES,50),STAT=ERR)
  ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
  
  ALLOCATE(CONNECTIVITY_NUMBER(TOTAL_NUMBER_OF_NODES),STAT=ERR)
  ALLOCATE(SEED_VALUE(TOTAL_NUMBER_OF_NODES),STAT=ERR)
  ALLOCATE(SNAPSHOT_SEED_VALUE(TOTAL_NUMBER_OF_NODES),STAT=ERR)
  ALLOCATE(RELATIVE_ERROR(TOTAL_NUMBER_OF_NODES),STAT=ERR)
  ALLOCATE(SEED_VALUE_EXACT_SOLUTION(TOTAL_NUMBER_OF_NODES),STAT=ERR)
  ALLOCATE(STATUS_MASK(TOTAL_NUMBER_OF_NODES),STAT=ERR)
  ALLOCATE(TRACE_NODE(TOTAL_NUMBER_OF_NODES),STAT=ERR)
  ALLOCATE(TRACE_NODES(TOTAL_NUMBER_OF_NODES),STAT=ERR)

  ALLOCATE(RAW_INDEX(TOTAL_NUMBER_OF_NODES+1),STAT=ERR)
  ALLOCATE(COLUMN_INDEX(TOTAL_NUMBER_OF_CONNECTIVITY),STAT=ERR)
  ALLOCATE(SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(TOTAL_NUMBER_OF_CONNECTIVITY,3),STAT=ERR)
  ALLOCATE(CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(TOTAL_NUMBER_OF_CONNECTIVITY,9),STAT=ERR)  
  
! MATERIAL_BEHAVIOUR : "ISOTROPIC"   -> isotropic material properties
!		               "ANISOTROPIC" -> anisotropic material properties
  MATERIAL_BEHAVIOUR = "ANISOTROPIC"

! if material is ANISOTROPIC
! INPUT_TYPE_FOR_CONDUCTIVITY : "TENSOR" -> three EigenVectors
!		                        "VECTOR" -> first EigenVector
  INPUT_TYPE_FOR_CONDUCTIVITY = "VECTOR" 

! INPUT_TYPE_FOR_SEED_VALUE : "FILE"  -> refers to file extension ".ESTM" with the same "INPUT_FILE_NAME" name
!		                      "LIST"  -> refers to NODAL LIST below
  INPUT_TYPE_FOR_SEED_VALUE = "LIST"

    ! SEED POINT 2
    SEED_POINT = 2
    STATUS_MASK(SEED_POINT)="SEED POINT"
    SEED_VALUE(SEED_POINT)=0.0_CMISSRP

! INPUT_TYPE_FOR_SPEED_FUNCTION : "FILE"   -> refers to file extension ".COND" with the same "INPUT_FILE_NAME" name
!		                          "FIXED"  -> refers to constant values below
  INPUT_TYPE_FOR_SPEED_FUNCTION = "FIXED" 

    ! enter the SPEED FUNCTIONS along EIGEN-VECTORS
    ! EIGEN-VECTOR 1
    SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1) = 10.0_CMISSRP

    ! EIGEN-VECTOR 2
    SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2) = 1.0_CMISSRP

    ! EIGEN-VECTOR 3
    SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3) = 1.0_CMISSRP

  CALL CMFE_PRE_PROCESS_INFORMATION(MATERIAL_BEHAVIOUR,INPUT_FILE_NAME,INPUT_FILE_FORMAT,TOTAL_NUMBER_OF_NODES,&
&INPUT_TYPE_FOR_SEED_VALUE,INPUT_TYPE_FOR_SPEED_FUNCTION,SPEED_FUNCTION_ALONG_EIGEN_VECTOR,INPUT_TYPE_FOR_CONDUCTIVITY,&
&STATUS_MASK,NODE_LIST,CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,SEED_VALUE,CONNECTIVITY_NUMBER,&
&SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,TOTAL_NUMBER_OF_CONNECTIVITY,&
&CONNECTIVITY_LIST,ELEMENT_LIST,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_NODES_PER_ELEMENT,Err)

  CALL CMFE_SOLVE_PROBLEM_FMM_CONNECTIVITY(TOTAL_NUMBER_OF_NODES,NODE_LIST,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,&
                       &SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,TOTAL_NUMBER_OF_CONNECTIVITY,&
                       &SEED_VALUE,STATUS_MASK)

! OUTPUT_FILE_FORMAT :  "TABC" -> tabc format
!			"VTKTET"  -> VTK tetrahedral format
  OUTPUT_FILE_FORMAT = "VTKTET"
  OUTPUT_FILE_NAME = "src/Output"
  OUTPUT_FILE_FIELD_TITLE = "depolarizatin_time_comp"

  CALL CMFE_POST_PROCESS_DATA(MATERIAL_BEHAVIOUR,OUTPUT_FILE_NAME,OUTPUT_FILE_FORMAT,TOTAL_NUMBER_OF_NODES,NODE_LIST,&
&CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,SEED_VALUE,CONNECTIVITY_NUMBER,OUTPUT_FILE_FIELD_TITLE,&
&CONNECTIVITY_LIST,ELEMENT_LIST,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_NODES_PER_ELEMENT,Err)
  
  STOP

!...................................................................................
END PROGRAM TwoElmnTet
!...................................................................................
