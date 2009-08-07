MODULE DARCY_PARAMS_CONVSTUDY

  USE DARCY_PARAMETERS

  IMPLICIT NONE

  INTEGER:: PDEGR
  INTEGER:: BC_NUMBER_OF_WALL_NODES, NUMBER_OF_BCS
  DOUBLE PRECISION:: PERM, VIS
  DOUBLE PRECISION:: X1, X2, Y1, Y2, Z1, Z2, GEOM_TOL
  DOUBLE PRECISION:: max_node_spacing
  CHARACTER*60 IN_CHAR
!   CHARACTER*90 PATH
  CHARACTER*90 NIMZ_dummy


  CONTAINS

  SUBROUTINE READ_DARCY_PARAMS

  IMPLICIT NONE

  WRITE(*,*)'Reading Darcy parameters.'

  OPEN(UNIT=42, FILE='./input/CMHEART.inp',STATUS='old')

  READ(42,*) NIMZ_dummy

  NIMZ_dummy = TRIM(NIMZ_dummy)

  OPEN(UNIT=37,FILE=NIMZ_dummy,STATUS='old',action='read') ! Read base file for initial parameters


      DO WHILE (0 < 1)
        READ(37,*,END=50) IN_CHAR
        IF (INDEX(IN_CHAR,'TESTCASE:') == 1)      READ(37,*) TESTCASE

        IF (INDEX(IN_CHAR,'STAB:') == 1)          READ(37,*) STAB
        IF (INDEX(IN_CHAR,'ANALYTIC:') == 1)      READ(37,*) ANALYTIC
        IF (INDEX(IN_CHAR,'DEBUG:') == 1)         READ(37,*) DEBUG

        IF (INDEX(IN_CHAR,'LENGTH:') == 1)        READ(37,*) LENGTH
        IF (INDEX(IN_CHAR,'GEOM_TOL:') == 1)      READ(37,*) GEOM_TOL
        IF (INDEX(IN_CHAR,'X1:') == 1)            READ(37,*) X1
        IF (INDEX(IN_CHAR,'X2:') == 1)            READ(37,*) X2
        IF (INDEX(IN_CHAR,'Y1:') == 1)            READ(37,*) Y1
        IF (INDEX(IN_CHAR,'Y2:') == 1)            READ(37,*) Y2
        IF (INDEX(IN_CHAR,'Z1:') == 1)            READ(37,*) Z1
        IF (INDEX(IN_CHAR,'Z2:') == 1)            READ(37,*) Z2
        IF (INDEX(IN_CHAR,'PERM:') == 1)          READ(37,*) PERM
        IF (INDEX(IN_CHAR,'VIS:') == 1)           READ(37,*) VIS

!         IF (INDEX(IN_CHAR,'PATH:') == 1) READ(37,*) PATH

        IF (INDEX(IN_CHAR,'PDEGR:') == 1) READ(37,*) PDEGR !polynomial degree

        IF (INDEX(IN_CHAR,'BC_NUMBER_OF_WALL_NODES:') == 1) READ(37,*) BC_NUMBER_OF_WALL_NODES
        IF (INDEX(IN_CHAR,'NUMBER_OF_BCS:') == 1) READ(37,*) NUMBER_OF_BCS
      END DO

50 CLOSE(37)

  IF( DEBUG ) THEN
    write(*,*)'Read from NIMZ_dummy = ',NIMZ_dummy
    write(*,*)'Press ENTER to continue.'
    read(*,*)
  END IF

  IF( ABS(VIS) > 1.0E-14 ) THEN
    PERM_OVER_VIS = PERM / VIS
  ELSE
    WRITE(*,*)'Darcy_parameters: VIS cannot be machine zero.'
    STOP
  END IF

  max_node_spacing = 2.0

  END SUBROUTINE READ_DARCY_PARAMS

END MODULE DARCY_PARAMS_CONVSTUDY
