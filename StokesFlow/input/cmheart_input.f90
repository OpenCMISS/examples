MODULE IMPORT_CMHEART

  IMPLICIT NONE

  !1=M, 2=V, 3=P !

  INCLUDE 'cheart_opencmiss.h'


  TYPE (ARRAY_PROBLEM_BASE) BASE_INFO
  TYPE (ARRAY_MESH) MESH_INFO(3)
  INTEGER I,J, FLD, DIMEN, OPENCMISS_INTERPOLATION(3),a,b
  INTEGER NumberOfNodesPerElement(3), NumberOfNodesDefined(3), NumberOfElementsDefined(3), TotalNumberOfNodes
  DOUBLE PRECISION, ALLOCATABLE::OPENCMISS_NODE_COORD(:,:)
  INTEGER, ALLOCATABLE::OPENCMISS_ELEM_M(:,:),OPENCMISS_ELEM_V(:,:),OPENCMISS_ELEM_P(:,:)
  CHARACTER*60 IN_CHAR
  CHARACTER*90 NIMZ
  CHARACTER*30 NAMz


  CONTAINS

  SUBROUTINE READ_CMHEART_EXE

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*)'Run from bin directory if input needed.'
  WRITE(*,*)
  WRITE(*,*)'REMEMBER: M,V,P need to be defined as required by cmHeart!'
  WRITE(*,*)
  WRITE(*,*)'Press ENTER to start'
  READ(*,*)

  OPEN(UNIT=42, FILE='./input/CMHEART.inp',STATUS='old')
!  OPEN(UNIT=42, FILE='./CMHEART.inp',STATUS='old')


  CALL READ_AUX
  CALL READ_NODES
  CALL READ_ELEMENTS

  ALLOCATE(OPENCMISS_ELEM_M(NumberOfElementsDefined(1),NumberOfNodesPerElement(1)))
  ALLOCATE(OPENCMISS_ELEM_V(NumberOfElementsDefined(2),NumberOfNodesPerElement(2)))
  ALLOCATE(OPENCMISS_ELEM_P(NumberOfElementsDefined(3),NumberOfNodesPerElement(3)))


  CALL MAKE_UNIQUE

  CALL ORDER_NUMBERING(OPENCMISS_ELEM_M,MESH_INFO(1)%T,NumberOfElementsDefined(1),NumberOfNodesPerElement(1),1)
  CALL ORDER_NUMBERING(OPENCMISS_ELEM_V,MESH_INFO(2)%T,NumberOfElementsDefined(2),NumberOfNodesPerElement(2),2) 
  CALL ORDER_NUMBERING(OPENCMISS_ELEM_P,MESH_INFO(3)%T,NumberOfElementsDefined(3),NumberOfNodesPerElement(3),3)  
  
 

  ! CALL PRINT_ON_SCREEN  

  WRITE(*,*)'export finished successfully...'
  WRITE(*,*)

  END SUBROUTINE READ_CMHEART_EXE

! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------


SUBROUTINE READ_AUX
      IMPLICIT NONE
  READ(42,*) NIMZ

  NIMZ = TRIM(NIMZ); BASE_INFO%n_B = 0; BASE_INFO%HEXA = 0; BASE_INFO%DM = 3

    OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read') ! Read base file for initial parameters

      DO WHILE (0 < 1)
        READ(1,*,END=50) IN_CHAR
        IF (INDEX(IN_CHAR,'no_fields!') == 1)         READ(1,*) BASE_INFO%n_B
        IF (INDEX(IN_CHAR,'no_gauss!') == 1)          READ(1,*) BASE_INFO%n_pts
        IF (INDEX(IN_CHAR,'volume!') == 1)            READ(1,*) BASE_INFO%VL
        IF (INDEX(IN_CHAR,'no_gauss_f!') == 1)        READ(1,*) BASE_INFO%n_pts_f
        IF (INDEX(IN_CHAR,'no_ele_faces!') == 1)      READ(1,*) BASE_INFO%FACES
        IF (INDEX(IN_CHAR,'no_ele_nodes_f!') == 1)    READ(1,*) BASE_INFO%FNODES
        IF (INDEX(IN_CHAR,'hexa_basis!') == 1)        BASE_INFO%HEXA = 1
        IF (INDEX(IN_CHAR,'domain_dimension!') == 1)  READ(1,*) BASE_INFO%DM
        IF (INDEX(IN_CHAR,'TRI_BASIS!') == 1)  BASE_INFO%TRI_BASIS = 1
        IF (INDEX(IN_CHAR,'TET_BASIS!') == 1)  BASE_INFO%TET_BASIS = 1
        IF (INDEX(IN_CHAR,'QUAD_BASIS!') == 1) BASE_INFO%QUAD_BASIS = 1
        IF (INDEX(IN_CHAR,'HEX_BASIS!') == 1)  BASE_INFO%HEX_BASIS = 1
      END DO

50 CLOSE(1)

    DIMEN=BASE_INFO%DM


! ! ! !         WRITE(*,*)'BASE_INFO%n_B',BASE_INFO%n_B
! ! ! !         WRITE(*,*)'BASE_INFO%n_pts',BASE_INFO%n_pts
! ! ! !         WRITE(*,*)'BASE_INFO%VL',BASE_INFO%VL
! ! ! !         WRITE(*,*)'BASE_INFO%n_pts_f',BASE_INFO%n_pts_f
! ! ! !         WRITE(*,*)'BASE_INFO%FACES',BASE_INFO%FACES
! ! ! !         WRITE(*,*)'BASE_INFO%FNODES',BASE_INFO%FNODES
! ! ! !         WRITE(*,*)'BASE_INFO%HEXA',BASE_INFO%HEXA
! ! ! !         WRITE(*,*)'BASE_INFO%DM',BASE_INFO%DM
! ! ! !         WRITE(*,*)'BASE_INFO%TRI_BASIS',BASE_INFO%TRI_BASIS
! ! ! !         WRITE(*,*)'BASE_INFO%TET_BASIS',BASE_INFO%TET_BASIS
! ! ! !         WRITE(*,*)'BASE_INFO%QUAD_BASIS',BASE_INFO%QUAD_BASIS
! ! ! !         WRITE(*,*)'BASE_INFO%HEX_BASIS',BASE_INFO%HEX_BASIS

        ALLOCATE(BASE_INFO%B(BASE_INFO%n_B))


    OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read')	! Read base file for initial parameters

      DO WHILE (0 < 1)
        READ(1,*,END=52) IN_CHAR
        IF (INDEX(IN_CHAR,'no_basis_M!') == 1) READ(1,*) BASE_INFO%B(1)%n
        IF (INDEX(IN_CHAR,'no_basis_V!') == 1) READ(1,*) BASE_INFO%B(2)%n
        IF (INDEX(IN_CHAR,'no_basis_P!') == 1) READ(1,*) BASE_INFO%B(3)%n

        IF (INDEX(IN_CHAR,'dim_field_M!') == 1) READ(1,*) BASE_INFO%B(1)%DM
        IF (INDEX(IN_CHAR,'dim_field_V!') == 1) READ(1,*) BASE_INFO%B(2)%DM
        IF (INDEX(IN_CHAR,'dim_field_P!') == 1) READ(1,*) BASE_INFO%B(3)%DM
      END DO
52 CLOSE(1)


IF (BASE_INFO%QUAD_BASIS/= 1.AND.BASE_INFO%HEX_BASIS /= 1.AND.BASE_INFO%TRI_BASIS/= 1.AND.BASE_INFO%TET_BASIS /= 1)THEN
  
    WRITE(*,*)'Cubic Hermite not implemented yet'
ELSE

DO I=1,3

        IF (BASE_INFO%TRI_BASIS== 1.OR.BASE_INFO%TET_BASIS == 1)THEN
            IF(DIMEN==2) THEN
              IF(BASE_INFO%B(I)%n==3) THEN
                    NumberOfNodesPerElement(I)=3
                    OPENCMISS_INTERPOLATION(I)=7
              ELSE IF(BASE_INFO%B(I)%n==6) THEN
                    NumberOfNodesPerElement(I)=6
                    OPENCMISS_INTERPOLATION(I)=8
              ELSE IF(BASE_INFO%B(I)%n==10) THEN
                    NumberOfNodesPerElement(I)=10
                    OPENCMISS_INTERPOLATION(I)=9
              ELSE
                  STOP
              END IF
           ELSE IF(DIMEN==3) THEN
              IF(BASE_INFO%B(I)%n==4) THEN
                    NumberOfNodesPerElement(I)=4
                    OPENCMISS_INTERPOLATION(I)=7
              ELSE IF(BASE_INFO%B(I)%n==10) THEN
                    NumberOfNodesPerElement(I)=10
                    OPENCMISS_INTERPOLATION(I)=8
              ELSE IF(BASE_INFO%B(I)%n==20) THEN
                    NumberOfNodesPerElement(I)=20
                    OPENCMISS_INTERPOLATION(I)=9
              ELSE
                  STOP
              END IF
              ELSE 
                STOP
           END IF
        END IF
 

        IF (BASE_INFO%QUAD_BASIS== 1.OR.BASE_INFO%HEX_BASIS == 1)THEN
            IF(BASE_INFO%B(I)%n==2) THEN
            !2D/3D LINEAR LAGRANGE
                  OPENCMISS_INTERPOLATION(I)=1
                  IF(DIMEN==2) THEN
                    NumberOfNodesPerElement(I)=4
                  ELSE IF(DIMEN==3) THEN
                    NumberOfNodesPerElement(I)=8
                  ELSE 
                    STOP
                  END IF



             ELSE IF(BASE_INFO%B(I)%n==3) THEN
            !2D/3D QUADRATIC LAGRANGE
                  OPENCMISS_INTERPOLATION(I)=2
                  IF(DIMEN==2) THEN
                    NumberOfNodesPerElement(I)=9
                  ELSE IF(DIMEN==3) THEN
                    NumberOfNodesPerElement(I)=27
                  ELSE 
                    STOP
                  END IF



            ELSE IF(BASE_INFO%B(I)%n==4) THEN
            !2D/3D CUBIC LAGRANGE
                  OPENCMISS_INTERPOLATION(I)=3
                  IF(DIMEN==2) THEN
                    NumberOfNodesPerElement(I)=16
                  ELSE IF(DIMEN==3) THEN
                    NumberOfNodesPerElement(I)=64
                  ELSE 
                    STOP
                  END IF

            ELSE
                  STOP
            END IF
           END IF

END DO



END IF

END SUBROUTINE READ_AUX

! ----------------------------------------------------------------------------------

SUBROUTINE ORDER_NUMBERING(NEW,OLD,n,m,I)

      IMPLICIT NONE
  INTEGER::n,m,I
  INTEGER::NEW(n,m),OLD(n,m)


!  NEW=OLD

      DO J=1,n
      IF (BASE_INFO%QUAD_BASIS == 1) THEN
        IF(BASE_INFO%B(I)%n==2) THEN
        !2D HEX LINEAR
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        ELSE IF(BASE_INFO%B(I)%n==3) THEN
        !2D HEX QUADR
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,5)
        NEW(J,3)=OLD(J,2)
        NEW(J,4)=OLD(J,6)
        NEW(J,5)=OLD(J,7)
        NEW(J,6)=OLD(J,8)
        NEW(J,7)=OLD(J,3)
        NEW(J,8)=OLD(J,9)
        NEW(J,9)=OLD(J,4)
        ELSE IF(BASE_INFO%B(I)%n==4) THEN
        !2D HEX CUB
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,5)
        NEW(J,3)=OLD(J,6)
        NEW(J,4)=OLD(J,2)
        NEW(J,5)=OLD(J,7)
        NEW(J,6)=OLD(J,8)
        NEW(J,7)=OLD(J,9)
        NEW(J,8)=OLD(J,10)
        NEW(J,9)=OLD(J,11)
        NEW(J,10)=OLD(J,12)
        NEW(J,11)=OLD(J,13)
        NEW(J,12)=OLD(J,14)
        NEW(J,13)=OLD(J,3)
        NEW(J,14)=OLD(J,15)
        NEW(J,15)=OLD(J,16)
        NEW(J,16)=OLD(J,4)
        ELSE
           STOP
        END IF

      ELSE IF (BASE_INFO%HEX_BASIS == 1) THEN


        IF(BASE_INFO%B(I)%n==2) THEN
        !3D HEX LINEAR
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        NEW(J,5)=OLD(J,5)
        NEW(J,6)=OLD(J,6)
        NEW(J,7)=OLD(J,7)
        NEW(J,8)=OLD(J,8)



        ELSE IF(BASE_INFO%B(I)%n==3) THEN
        !3D HEX QUADR
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,9)
        NEW(J,3)=OLD(J,2)
        NEW(J,4)=OLD(J,10)
        NEW(J,5)=OLD(J,11)
        NEW(J,6)=OLD(J,12)
        NEW(J,7)=OLD(J,3)
        NEW(J,8)=OLD(J,13)
        NEW(J,9)=OLD(J,4)
        NEW(J,10)=OLD(J,14)
        NEW(J,11)=OLD(J,15)
        NEW(J,12)=OLD(J,16)
        NEW(J,13)=OLD(J,17)
        NEW(J,14)=OLD(J,18)
        NEW(J,15)=OLD(J,19)
        NEW(J,16)=OLD(J,20)
        NEW(J,17)=OLD(J,21)
        NEW(J,18)=OLD(J,22)
        NEW(J,19)=OLD(J,5)
        NEW(J,20)=OLD(J,23)
        NEW(J,21)=OLD(J,6)
        NEW(J,22)=OLD(J,24)
        NEW(J,23)=OLD(J,25)
        NEW(J,24)=OLD(J,26)
        NEW(J,25)=OLD(J,7)
        NEW(J,26)=OLD(J,27)
        NEW(J,27)=OLD(J,8)


        ELSE IF(BASE_INFO%B(I)%n==4) THEN
    !3D HEX CUB
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,9)
        NEW(J,3)=OLD(J,10)
        NEW(J,4)=OLD(J,2)
        NEW(J,5)=OLD(J,11)
        NEW(J,6)=OLD(J,12)
        NEW(J,7)=OLD(J,13)
        NEW(J,8)=OLD(J,14)
        NEW(J,9)=OLD(J,15)
        NEW(J,10)=OLD(J,16)
        NEW(J,11)=OLD(J,17)
        NEW(J,12)=OLD(J,18)
        NEW(J,13)=OLD(J,3)
        NEW(J,14)=OLD(J,19)
        NEW(J,15)=OLD(J,20)
        NEW(J,16)=OLD(J,4)
        NEW(J,17)=OLD(J,21)
        NEW(J,18)=OLD(J,22)
        NEW(J,19)=OLD(J,23)
        NEW(J,20)=OLD(J,24)
        NEW(J,21)=OLD(J,25)
        NEW(J,22)=OLD(J,26)
        NEW(J,23)=OLD(J,27)
        NEW(J,24)=OLD(J,28)
        NEW(J,25)=OLD(J,29)
        NEW(J,26)=OLD(J,30)
        NEW(J,27)=OLD(J,31)
        NEW(J,28)=OLD(J,32)
        NEW(J,29)=OLD(J,33)
        NEW(J,30)=OLD(J,34)
        NEW(J,31)=OLD(J,35)
        NEW(J,32)=OLD(J,36)
        NEW(J,33)=OLD(J,37)
        NEW(J,34)=OLD(J,38)
        NEW(J,35)=OLD(J,39)
        NEW(J,36)=OLD(J,40)
        NEW(J,37)=OLD(J,41)
        NEW(J,38)=OLD(J,42)
        NEW(J,39)=OLD(J,43)
        NEW(J,40)=OLD(J,44)
        NEW(J,41)=OLD(J,45)
        NEW(J,42)=OLD(J,46)
        NEW(J,43)=OLD(J,47)
        NEW(J,44)=OLD(J,48)
        NEW(J,45)=OLD(J,49)
        NEW(J,46)=OLD(J,50)
        NEW(J,47)=OLD(J,51)
        NEW(J,48)=OLD(J,52)
        NEW(J,49)=OLD(J,5)
        NEW(J,50)=OLD(J,53)
        NEW(J,51)=OLD(J,54)
        NEW(J,52)=OLD(J,6)
        NEW(J,53)=OLD(J,55)
        NEW(J,54)=OLD(J,56)
        NEW(J,55)=OLD(J,57)
        NEW(J,56)=OLD(J,58)
        NEW(J,57)=OLD(J,59)
        NEW(J,58)=OLD(J,60)
        NEW(J,59)=OLD(J,61)
        NEW(J,60)=OLD(J,62)
        NEW(J,61)=OLD(J,7)
        NEW(J,62)=OLD(J,63)
        NEW(J,63)=OLD(J,64)
        NEW(J,64)=OLD(J,8)        

        ELSE
           STOP
        END IF

      ELSE IF (BASE_INFO%TRI_BASIS == 1) THEN
        IF(BASE_INFO%B(I)%n==3) THEN
        !2D TET LINEAR
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        ELSE IF(BASE_INFO%B(I)%n==6) THEN
        !2D TET QUAD
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        NEW(J,5)=OLD(J,6)
        NEW(J,6)=OLD(J,5)
        ELSE IF(BASE_INFO%B(I)%n==10) THEN
        !2D TET CUB
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        NEW(J,5)=OLD(J,5)
        NEW(J,6)=OLD(J,8)
        NEW(J,7)=OLD(J,9)
        NEW(J,8)=OLD(J,7)
        NEW(J,9)=OLD(J,6)
        NEW(J,10)=OLD(J,10)
        ELSE
          STOP
        END IF

      ELSE IF (BASE_INFO%TET_BASIS == 1) THEN
        IF(BASE_INFO%B(I)%n==4) THEN
        !3D TET LINEAR
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        ELSE IF(BASE_INFO%B(I)%n==10) THEN
        !3D TET QUAD
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        NEW(J,5)=OLD(J,5)
        NEW(J,6)=OLD(J,6)
        NEW(J,7)=OLD(J,7)
        NEW(J,8)=OLD(J,8)
        NEW(J,9)=OLD(J,10)
        NEW(J,10)=OLD(J,9)

        ELSE IF(BASE_INFO%B(I)%n==20) THEN
        !3D TET CUB
        NEW(J,1)=OLD(J,1)
        NEW(J,2)=OLD(J,2)
        NEW(J,3)=OLD(J,3)
        NEW(J,4)=OLD(J,4)
        NEW(J,5)=OLD(J,5)
        NEW(J,6)=OLD(J,6)
        NEW(J,7)=OLD(J,7)
        NEW(J,8)=OLD(J,8)
        NEW(J,9)=OLD(J,9)
        NEW(J,10)=OLD(J,10)
        NEW(J,11)=OLD(J,11)
        NEW(J,12)=OLD(J,12)
        NEW(J,13)=OLD(J,15)
        NEW(J,14)=OLD(J,16)
        NEW(J,15)=OLD(J,13)
        NEW(J,16)=OLD(J,14)
        NEW(J,17)=OLD(J,17)
        NEW(J,18)=OLD(J,18)
        NEW(J,19)=OLD(J,19)
        NEW(J,20)=OLD(J,20)
        ELSE
          STOP
        END IF
        
       ELSE
        STOP
       END IF

    END DO



	

END SUBROUTINE ORDER_NUMBERING

! ----------------------------------------------------------------------------------

SUBROUTINE MAKE_UNIQUE
      IMPLICIT NONE

! NOW, THE NODE NUMBERING NEEDS TO BE CHANGED FOR ALL QUADRATIC AND CUBIC ELEMENTS


! NOW, SAME NODES NEED MAKE_UNIQUE/IDENTICAL DEFINITION

! M is considered as reference and check V

    MESH_INFO(1)%T=MESH_INFO(1)%T
    MESH_INFO(2)%T=MESH_INFO(2)%T+NumberOfNodesDefined(1)
    MESH_INFO(3)%T=MESH_INFO(3)%T+NumberOfNodesDefined(1)+NumberOfNodesDefined(2)

    IF(NumberOfNodesDefined(1)==NumberOfNodesDefined(2)) THEN
    ! copy all node numbers from 2 -> 1
    MESH_INFO(2)%T(:,:)= MESH_INFO(1)%T(:,:)
    ELSE
      IF (BASE_INFO%TRI_BASIS== 1) THEN
        MESH_INFO(2)%T(:,1:3)=MESH_INFO(1)%T(:,1:3)
      ELSE IF (BASE_INFO%TET_BASIS == 1)THEN
        MESH_INFO(2)%T(:,1:4)=MESH_INFO(1)%T(:,1:4)
      ELSE IF (BASE_INFO%QUAD_BASIS == 1)THEN
        MESH_INFO(2)%T(:,1:4)=MESH_INFO(1)%T(:,1:4)
      ELSE IF (BASE_INFO%HEX_BASIS == 1)THEN
        MESH_INFO(2)%T(:,1:8)=MESH_INFO(1)%T(:,1:8)
      ELSE
        STOP
      END IF
    END IF

    IF(NumberOfNodesDefined(1)==NumberOfNodesDefined(3)) THEN
    ! copy all node numbers from 2 -> 1
    MESH_INFO(3)%T(:,:)=MESH_INFO(1)%T(:,:)
    ELSE IF(NumberOfNodesDefined(2)==NumberOfNodesDefined(3)) THEN
    MESH_INFO(3)%T(:,:)=MESH_INFO(2)%T(:,:)
    ELSE
      IF (BASE_INFO%TRI_BASIS== 1) THEN
        MESH_INFO(3)%T(:,1:3)=MESH_INFO(1)%T(:,1:3)
      ELSE IF (BASE_INFO%TET_BASIS == 1)THEN
        MESH_INFO(3)%T(:,1:4)=MESH_INFO(1)%T(:,1:4)
      ELSE IF (BASE_INFO%QUAD_BASIS == 1)THEN
        MESH_INFO(3)%T(:,1:4)=MESH_INFO(1)%T(:,1:4)
      ELSE IF (BASE_INFO%HEX_BASIS == 1)THEN
        MESH_INFO(3)%T(:,1:8)=MESH_INFO(1)%T(:,1:8)
      ELSE
        STOP
      END IF
    END IF

END SUBROUTINE MAKE_UNIQUE

! ----------------------------------------------------------------------------------

SUBROUTINE READ_NODES

      IMPLICIT NONE

      INTEGER:: a,b

      READ(42,*) NAMz

      OPEN(UNIT = 1, FILE=NAMz,STATUS='old')
      READ(1,*) NumberOfNodesDefined(1:3)

      TotalNumberOfNodes=NumberOfNodesDefined(1)+NumberOfNodesDefined(2)+NumberOfNodesDefined(3)


! ALLOCATE AND READ MESH NODE INFORMATION

  WRITE(*,*)'Reading Nodes...'

  DO I=1,3
      MESH_INFO(I)%Lx=NumberOfNodesDefined(I)
      ALLOCATE(MESH_INFO(I)%X(MESH_INFO(I)%Lx,3))
        DO J = 1,MESH_INFO(I)%Lx
          READ(1,*,END=35) MESH_INFO(I)%X(J,1:3)
        END DO
  END DO




      CLOSE(1)

    ALLOCATE(OPENCMISS_NODE_COORD(TotalNumberOfNodes,3))
    a=1
    b=0


    DO I=1,3
    a=b+1
    b=b+NumberOfNodesDefined(I)
    OPENCMISS_NODE_COORD(a:b,1:3)=MESH_INFO(I)%X(1:NumberOfNodesDefined(I),1:3)
    END DO





        RETURN

 35     PRINT *, 'FAILS'
        STOP



END SUBROUTINE READ_NODES

! ----------------------------------------------------------------------------------

SUBROUTINE READ_ELEMENTS

      IMPLICIT NONE

      READ(42,*) NAMz
      CLOSE(42)

      OPEN(UNIT = 1, FILE=NAMz,STATUS='old')
      READ(1,*) NumberOfElementsDefined(1:3)

! ALLOCATE AND READ MESH ELEMENT INFORMATION

      TotalNumberOfNodes=NumberOfNodesDefined(1)+NumberOfNodesDefined(2)+NumberOfNodesDefined(3)


  WRITE(*,*)'Reading Elements...'
  WRITE(*,*)



  DO I=1,3
      MESH_INFO(I)%Lt=NumberOfElementsDefined(I)
      ALLOCATE(MESH_INFO(I)%T(MESH_INFO(I)%Lt,NumberOfNodesPerElement(I)))
        DO J = 1,MESH_INFO(I)%Lt
          READ(1,*,END=30) MESH_INFO(I)%T(J,1:NumberOfNodesPerElement(I))
        END DO
  END DO
      CLOSE(1)

        RETURN

 30     PRINT *, 'FAILS'
        STOP


END SUBROUTINE READ_ELEMENTS

! ----------------------------------------------------------------------------------

SUBROUTINE PRINT_ON_SCREEN

      IMPLICIT NONE

  ! where are the node coordinates stored -> 1 MATRIX



     DO I = 1,TotalNumberOfNodes
               WRITE(*,'("Node ",(I0,4x),1000( F5.3,2x ))')I,OPENCMISS_NODE_COORD(I,1:3)
     END DO


  ! where are the element nodes stored -> 3 MATRICES


    WRITE(*,*)
    WRITE(*,*)
  
    DO I = 1,NumberOfElementsDefined(1)
    WRITE(*,'("M-Elements: ", (I0,3x), (1000(I0, 1x)) )')I, &
      & OPENCMISS_ELEM_M(I,1:NumberOfNodesPerElement(1))
    END DO

    WRITE(*,*)

    DO I = 1,NumberOfElementsDefined(2)
    WRITE(*,'("V-Elements: ", (I0,3x), (1000(I0, 1x)) )')I, &
      & OPENCMISS_ELEM_V(I,1:NumberOfNodesPerElement(2))
    END DO

    WRITE(*,*)

    DO I = 1,NumberOfElementsDefined(3)
    WRITE(*,'("P-Elements: ", (I0,3x), (1000(I0, 1x)) )')I, &
      & OPENCMISS_ELEM_P(I,1:NumberOfNodesPerElement(3))
    END DO


    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*)


END SUBROUTINE PRINT_ON_SCREEN

! ----------------------------------------------------------------------------------

SUBROUTINE RECV_CMHEART_EXE(EXPORT)

  TYPE (EXPORT_CONTAINER):: EXPORT  
  TYPE (EXPORT_CONTAINER):: TMP  

  ALLOCATE(TMP%M(NumberOfElementsDefined(1),NumberOfNodesPerElement(1)))
  ALLOCATE(TMP%V(NumberOfElementsDefined(2),NumberOfNodesPerElement(2)))
  ALLOCATE(TMP%P(NumberOfElementsDefined(3),NumberOfNodesPerElement(3)))
  ALLOCATE(TMP%N(TotalNumberOfNodes,3))


  TMP%M=OPENCMISS_ELEM_M
  TMP%V=OPENCMISS_ELEM_V
  TMP%P=OPENCMISS_ELEM_P

  TMP%N=OPENCMISS_NODE_COORD

  TMP%D=DIMEN	
  TMP%F=BASE_INFO%n_B

  TMP%ID_M=1
  TMP%ID_V=2
  TMP%ID_P=3

  TMP%IT_M=OPENCMISS_INTERPOLATION(1)
  TMP%IT_V=OPENCMISS_INTERPOLATION(2)
  TMP%IT_P=OPENCMISS_INTERPOLATION(3)
  
  IF (BASE_INFO%HEXA==1) THEN
    !LAGRANGIAN BASIS
    TMP%IT_T=1
  ELSE 
    ! SIMPLEX BASIS
    TMP%IT_T=2
  END IF

  TMP%E_M=NumberOfElementsDefined(1)
  TMP%E_V=NumberOfElementsDefined(2)
  TMP%E_P=NumberOfElementsDefined(3)
  TMP%E_T=NumberOfElementsDefined(3)

  TMP%EN_M=NumberOfNodesPerElement(1)
  TMP%EN_V=NumberOfNodesPerElement(2)
  TMP%EN_P=NumberOfNodesPerElement(3)
  TMP%EN_T=TMP%EN_M+TMP%EN_V+TMP%EN_P


  TMP%N_M=NumberOfNodesDefined(1)
  TMP%N_V=NumberOfNodesDefined(2)
  TMP%N_P=NumberOfNodesDefined(3)
  TMP%N_T=TMP%N_M+TMP%N_V+TMP%N_P

  EXPORT=TMP

END SUBROUTINE RECV_CMHEART_EXE

! ----------------------------------------------------------------------------------

END MODULE IMPORT_CMHEART