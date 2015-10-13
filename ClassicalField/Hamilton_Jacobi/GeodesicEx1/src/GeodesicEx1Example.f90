!> \file
!> \author Ali PASHAEI
!> \brief This is an example program to find the geodesic using OpenCMISS calls.
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

!> \example ClassicalField/Hamilton_Jacobi/GeodesicEx1/src/GeodesicEx1Example.f90
!! Example program to solve a Laplace equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/Laplace/Laplace/history.html
!!
!<

PROGRAM GeodesicEx1
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

  INTEGER(CMISSIntg) :: I,J,SEED_POINT,TOTAL_NUMBER_OF_CONNECTIVITY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_PER_ELEMENT,TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: NODE_1,NODE_2
  REAL(CMISSRP) :: SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)
  CHARACTER (LEN=300) :: INPUT_FILE_NAME
  CHARACTER (LEN=300) :: OUTPUT_FILE_NAME
  CHARACTER (LEN=300) :: OUTPUT_FILE_FIELD_TITLE
  CHARACTER (LEN=10)  :: INPUT_FILE_FORMAT
  CHARACTER (LEN=10)  :: OUTPUT_FILE_FORMAT
  CHARACTER (LEN=12)  :: MATERIAL_BEHAVIOUR
  CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_SPEED_FUNCTION
  CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_SEED_VALUE
  CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_CONDUCTIVITY

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
  INPUT_FILE_NAME = "src/GeodesicEx1"

! Geodesic points :   Node_1 and Node_2 are Node Number of two ends of the geodesic path in the domain.
  NODE_1=451
  NODE_2=1

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
  
  DO I=1,TOTAL_NUMBER_OF_NODES
    TRACE_NODES(I)=0.0
  ENDDO

  MATERIAL_BEHAVIOUR = "ISOTROPIC"
  INPUT_TYPE_FOR_CONDUCTIVITY = "VECTOR" 
  INPUT_TYPE_FOR_SEED_VALUE = "LIST"

  SEED_POINT = NODE_2
  STATUS_MASK(SEED_POINT)="SEED POINT"
  SEED_VALUE(SEED_POINT)=0.0_CMISSRP
  TRACE_NODE(SEED_POINT)=0

  INPUT_TYPE_FOR_SPEED_FUNCTION = "FIXED" 
  SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1) = 1.0_CMISSRP


  CALL CMFE_PRE_PROCESS_INFORMATION(MATERIAL_BEHAVIOUR,INPUT_FILE_NAME,INPUT_FILE_FORMAT,TOTAL_NUMBER_OF_NODES,&
&INPUT_TYPE_FOR_SEED_VALUE,INPUT_TYPE_FOR_SPEED_FUNCTION,SPEED_FUNCTION_ALONG_EIGEN_VECTOR,INPUT_TYPE_FOR_CONDUCTIVITY,&
&STATUS_MASK,NODE_LIST,CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,SEED_VALUE,CONNECTIVITY_NUMBER,&
&SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,&
&TOTAL_NUMBER_OF_CONNECTIVITY,CONNECTIVITY_LIST,ELEMENT_LIST,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_NODES_PER_ELEMENT,Err)

  CALL CMFE_SOLVE_PROBLEM_GEODESIC_CONNECTIVITY(TOTAL_NUMBER_OF_NODES,NODE_LIST,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,&
                       &SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,TOTAL_NUMBER_OF_CONNECTIVITY,&
                       &SEED_VALUE,STATUS_MASK,TRACE_NODE)

  OUTPUT_FILE_FORMAT = "VTKTET"
  OUTPUT_FILE_NAME = "src/Output"
  OUTPUT_FILE_FIELD_TITLE = "trace"
      
  J=NODE_1
  TRACE_NODES(J)=1.0
  DO WHILE (J .NE. NODE_2)
    J=TRACE_NODE(J)
    TRACE_NODES(J)=1.0
  ENDDO
  
  CALL CMFE_POST_PROCESS_DATA(MATERIAL_BEHAVIOUR,OUTPUT_FILE_NAME,OUTPUT_FILE_FORMAT,TOTAL_NUMBER_OF_NODES,NODE_LIST,&
&CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,TRACE_NODES,CONNECTIVITY_NUMBER,OUTPUT_FILE_FIELD_TITLE,&
&CONNECTIVITY_LIST,ELEMENT_LIST,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_NODES_PER_ELEMENT,Err)
  
  STOP

!...................................................................................
END PROGRAM GeodesicEx1
!...................................................................................
