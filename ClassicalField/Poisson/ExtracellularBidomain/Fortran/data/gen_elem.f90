        PROGRAM GENERATE_ELEMENTS
        IMPLICIT NONE

        INTEGER :: START_VAL,ELEM,X,I,J,K,T,M,A,P,Q,INC_ELEM_VER,INC_START
        INTEGER :: N1,N2,N3,N4,N5,N6,N7,N8,LH,LV,BH,BV,NUM_CUBES,NUM_CUBES_VER,NUM_INTER,NUM_ELEM,NUM_NODES
        INTEGER :: N1_EL1,N2_EL1,N3_EL1,N4_EL1,N5_EL1,N6_EL1,N7_EL1,N8_EL1


        !DEFINITON OF VARIABLES
        !X: NUMBER OF NODES PER LINE(FIBER)
        !LH: NUMBER OF LINES PER BLOCK IN HORIZONTAL DIRECTION
        !LV: NUMBER OF LINES PER BLOCK IN VERTICAL DIRECTION
        !BH: NUMBER OF BLOCKS IN HORIZONTAL DIRECTION
        !BV: NUMBER OF BLOCKS IN VERTICAL DIRECTION
        !NUM_CUBES: NUMBER OF BLOCKS IN THE RECURRING PATTERN
        !ELEM: ELEMENT NUMBER

        !VARIABLES WHICH CAN BE ALTERED BY THE USER
        LH=3
        LV=3
        BH=1
        BV=1
        X=349
        !END ---------------------------------------


        !INITIALIZATION
        START_VAL=0
        N1=0
        N2=0
        N3=0
        N4=0
        N5=0
        N6=0
        N7=0
        N8=0
        I=1
        J=1
        K=1
        Q=1
        P=1
        T=1
        M=1
        A=1
        INC_ELEM_VER=0
        INC_START=0

        NUM_CUBES=BH*BV
        WRITE (*,*) NUM_CUBES


        OPEN(UNIT=2,FILE='DEBUG.TXT')

        OPEN(UNIT=3,FILE='ELEMENTNODES.TXT')


        NUM_ELEM=((LH*BH-1)*(LV*BV-1))*(X-1)
        NUM_NODES=(LH*BH)*(LV*BV)*X
        WRITE (3,*) NUM_ELEM, NUM_NODES
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! PART 1 - RECURRING BLOCKS                                !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111111
        DO !BLOCK WISE
           IF (K.GT.NUM_CUBES) EXIT
           DO !ROW-WISE
              IF (M.GT.(LV-1)) EXIT
              WRITE (*,*) 'M', M
              INC_ELEM_VER=(LV-1)*(M-1)*(X-1)
              INC_START=LH*X*(M-1)
              WRITE (*,*) 'INC', INC_ELEM_VER,INC_START
  
              DO !COLUMN-WISE
                 IF (J.GT.(LH-1)) EXIT
                 ELEM=(LH-1)*(LV-1)*(K-1)*(X-1)+(J-1)*(X-1)+INC_ELEM_VER+1
                 START_VAL=LH*LV*X*(K-1)+(J-1)*X+INC_START+1
                 N1_EL1=START_VAL
                 N3_EL1=N1_EL1+X
                 N5_EL1=N1_EL1+LH*X
                 N7_EL1=N5_EL1+X
                 N6_EL1=N5_EL1+1
                 N8_EL1=N7_EL1+1
                 N2_EL1=N1_EL1+1
                 N4_EL1=N3_EL1+1
      
                 N1=N1_EL1
 
                 WRITE(3,*) ELEM,N1_EL1,N2_EL1,N3_EL1,N4_EL1,N5_EL1,N6_EL1,N7_EL1,N8_EL1
                 WRITE (*,*) 'START', START_VAL
                 ELEM=ELEM+1
                 DO !FIBER-WISE
                    IF (T.GT.(X-2)) EXIT !RUN THIS LOOP FOR (NUMBER OF ELEMENTS IN ONE FIBER DIRECTION - 1) = [(NUM OF NODES IN FIBER DIR - 1) -1]
                    N1=N1+1
                    N3=N1+X
                    N5=N1+LH*X
                    N7=N5+X
                    N6=N5+1
                    N8=N7+1
                    N2=N1+1
                    N4=N3+1
                    T=T+1
        
                    WRITE(3,*) ELEM,N1,N2,N3,N4,N5,N6,N7,N8
                    ELEM=ELEM+1
                    END DO !FIBER-WISE
                J=J+1
                T=1
                END DO !COLUMN-WISE
              M=M+1
              J=1
              WRITE(*,*) 'M', M
              END DO !ROW-WISE
        T=1
        M=1
        K=K+1
        END DO !BLOCK-WISE


        !!!! PART 1 - END

        !ELEM VARIABLE IS ALREADY INITIALIZED FOR THE NEXT PART

        !INITALIZE ALL VARIABLES THAT WILL BE USED FOR PART 2 AGAIN
        K=1
        M=1
        J=1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! PART 2- VERTICAL INTERSECTIONS                           !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !WHOLE SYSTEM
        DO !VERTICAL BLOCK
           IF (Q.GT.BV) EXIT
           DO  !HORIZONTAL BLOCK
             IF (P.GT.(BH-1)) EXIT
             !BLOCK WISE
             DO !ROW-WISE
              IF (M.GT.(LV-1)) EXIT
              INC_ELEM_VER=(LV-1)*(K-1)*(X-1)+(M-1)*(X-1)
              INC_START=(A-1)*(LH*LV)*X+(M-1)*LH*X
              
                 !COLUMN-WISE
                 ELEM=(NUM_CUBES*(LH-1)*(LV-1)*(X-1)+1)+INC_ELEM_VER
                 START_VAL=(LH-1)*X+1+INC_START
                 WRITE(2,*) 'START VAL', START_VAL
                 WRITE(2,*) 'ELEM',ELEM
                 N1_EL1=START_VAL
                 N3_EL1=N1_EL1+(LH*(LV-1)+1)*X
                 N5_EL1=N1_EL1+LH*X
                 N7_EL1=N5_EL1+(LH*(LV-1)+1)*X
                 N6_EL1=N5_EL1+1
                 N8_EL1=N7_EL1+1
                 N2_EL1=N1_EL1+1
                 N4_EL1=N3_EL1+1

                 N1=N1_EL1

                 WRITE(3,*) ELEM,N1_EL1,N2_EL1,N3_EL1,N4_EL1,N5_EL1,N6_EL1,N7_EL1,N8_EL1
                 WRITE (2,*) 'START', START_VAL
                 ELEM=ELEM+1

                    DO !FIBER-WISE
                       IF (T.GT.(X-2)) EXIT !RUN THIS LOOP FOR (NUMBER OF ELEMENTS IN ONE FIBER DIRECTION - 1) = [(NUM OF NODES IN FIBER DIR - 1) -1]
                       N1=N1+1
                       N3=N1+(LH*(LV-1)+1)*X
                       N5=N1+LH*X
                       N7=N5+(LH*(LV-1)+1)*X
                       N6=N5+1
                       N8=N7+1
                       N2=N1+1
                       N4=N3+1
                       T=T+1

                       WRITE(3,*) ELEM,N1,N2,N3,N4,N5,N6,N7,N8
                       ELEM=ELEM+1
                       END DO !FIBER-WISE
                J=J+1
                T=1
                M=M+1
                !END - COLUMN-WISE
                J=1
                WRITE(2,*) 'M', M
                END DO !ROW-WISE
                ! END - BLOCK-WISE
             T=1
             M=1
             K=K+1
             P=P+1
             A=A+1
             END DO !HORIZONTAL
        P=1
        Q=Q+1
        A=A+1
        END DO !VERTICAL
        ! END - WHOLE SYSTEM
        
        Q=1
        P=1
        J=1
        A=1

        !!! PART 2 - END
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! PART 3- HORIZONTAL INTERSECTIONS                           !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! NUM_CUBES_VER=(BH-1)*BV

        !WHOLE SYSTEM
        DO !VERTICAL BLOCK
           IF (Q.GT.(BV-1)) EXIT
           WRITE (2,*) 'Q', Q
           DO  !HORIZONTAL BLOCK
             IF (P.GT.BH) EXIT
             WRITE (2,*) 'P', P
             !BLOCK (CUBE) WISE
             !ROW-WISE
               DO !COLUMN-WISE
               IF (J.GT.(LV-1)) EXIT
                 INC_ELEM_VER=(M-1)*(X-1)
                 INC_START=LH*LV*(A-1)*X+(J-1)*X
                 WRITE (2,*) 'A', A
                 WRITE (2,*) 'INC', INC_ELEM_VER,INC_START

                 ELEM=(NUM_CUBES*(LH-1)*(LV-1)+(BH-1)*BV*(LV-1))*(X-1)+1+INC_ELEM_VER
                 START_VAL=LH*(LV-1)*X+1+INC_START
                 WRITE(2,*) 'START VAL', START_VAL
                 WRITE(2,*) 'ELEM',ELEM
                 N1_EL1=START_VAL
                 N3_EL1=N1_EL1+X
                 N5_EL1=N1_EL1+(LH*(LV+1))*X
                 N7_EL1=N5_EL1+X
                 N6_EL1=N5_EL1+1
                 N8_EL1=N7_EL1+1
                 N2_EL1=N1_EL1+1
                 N4_EL1=N3_EL1+1

                 N1=N1_EL1

                 WRITE(3,*) ELEM,N1_EL1,N2_EL1,N3_EL1,N4_EL1,N5_EL1,N6_EL1,N7_EL1,N8_EL1
                 WRITE (2,*) 'START', START_VAL
                 ELEM=ELEM+1

                    DO !FIBER-WISE
                       IF (T.GT.(X-2)) EXIT !RUN THIS LOOP FOR (NUMBER OF ELEMENTS IN ONE FIBER DIRECTION - 1) = [(NUM OF NODES IN FIBER DIR - 1) -1]
                       N1=N1+1
                       N3=N1+X
                       N5=N1+(LH*(LV+1))*X
                       N7=N5+X
                       N6=N5+1
                       N8=N7+1
                       N2=N1+1
                       N4=N3+1
                       T=T+1

                       WRITE(3,*) ELEM,N1,N2,N3,N4,N5,N6,N7,N8
                       ELEM=ELEM+1
                       END DO !FIBER-WISE
                J=J+1
                T=1
                M=M+1
                END DO !COLUMN-WISE
                !END - ROW-WISE
                !END - BLOCK-WISE
             J=1
             T=1
             K=K+1
             P=P+1
             A=A+1
             END DO !HORIZONTAL

        P=1
        Q=Q+1
        END DO !VERTICAL
        ! END - WHOLE SYSTEM

        Q=1
        P=1
        J=1
        M=1

        !!! PART 3 - END
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! PART 4- INTERSECTIONS                           !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        NUM_INTER=(BH-1)*(BV-1)

        WRITE(2,*) 'BEGIN PART 4'

        DO
         IF (Q.GT.(BV-1)) EXIT
          DO
          IF (P.GT.(BH-1)) EXIT
          !IF (M.GT.NUM_INTER) EXIT
                 INC_ELEM_VER=(J-1)*(X-1)
                 INC_START=((M-1)*LH*LV)*X
                 WRITE (2,*) 'INC', INC_ELEM_VER,INC_START

                 ELEM=((LH*BH-1)*(LV*BV-1)-NUM_INTER)*(X-1)+1+INC_ELEM_VER
                 START_VAL=(LH*LV-1)*X+1+INC_START
                 WRITE(2,*) 'START VAL', START_VAL
                 WRITE(2,*) 'ELEM',ELEM
                 N1_EL1=START_VAL
                 N3_EL1=N1_EL1+(LH*(LV-1)+1)*X
                 N5_EL1=N1_EL1+(LH*LV+LH)*X
                 N7_EL1=N5_EL1+(LH*(LV-1)+1)*X
                 N6_EL1=N5_EL1+1
                 N8_EL1=N7_EL1+1
                 N2_EL1=N1_EL1+1
                 N4_EL1=N3_EL1+1

                 N1=N1_EL1

                 WRITE(3,*) ELEM,N1_EL1,N2_EL1,N3_EL1,N4_EL1,N5_EL1,N6_EL1,N7_EL1,N8_EL1
                 WRITE (2,*) 'START', START_VAL
                 ELEM=ELEM+1

                    DO !FIBER-WISE
                       IF (T.GT.(X-2)) EXIT !RUN THIS LOOP FOR (NUMBER OF ELEMENTS IN ONE FIBER DIRECTION - 1) = [(NUM OF NODES IN FIBER DIR - 1) -1]
                       N1=N1+1
                       N3=N1+(LH*(LV-1)+1)*X
                       N5=N1+(LH*LV+LH)*X
                       N7=N5+(LH*(LV-1)+1)*X
                       N6=N5+1
                       N8=N7+1
                       N2=N1+1
                       N4=N3+1
                       WRITE(3,*) ELEM,N1,N2,N3,N4,N5,N6,N7,N8
                       ELEM=ELEM+1
                       T=T+1
                    END DO !FIBER-WISE
                T=1
                P=P+1
                M=M+1
                A=A+1
              END DO
            P=1
            Q=Q+1
            M=M+1
            J=J+1
            END DO

        !!! PART 4 - END

         END PROGRAM GENERATE_ELEMENTS
