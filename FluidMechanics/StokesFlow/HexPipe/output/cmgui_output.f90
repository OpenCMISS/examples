MODULE EXPORT_CMGUI

  USE BASE_ROUTINES
  USE LISTS
  USE BASIS_ROUTINES
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE COMP_ENVIRONMENT
  USE COORDINATE_ROUTINES
  USE ISO_VARYING_STRING
  USE REGION_ROUTINES
  USE MACHINE_CONSTANTS
  USE KINDS
  USE FIELD_ROUTINES
  USE ISO_VARYING_STRING
  USE ISO_C_BINDING
  USE STRINGS
  USE TYPES
  USE CONSTANTS
  USE MPI
  USE CMISS_MPI
  USE INPUT_OUTPUT

  IMPLICIT NONE

  !1=M, 2=V, 3=P !


  INTEGER, DIMENSION(:), ALLOCATABLE:: NodesPerElement
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE::ElementNodesScales
  INTEGER, DIMENSION(:,:), ALLOCATABLE::ElementNodes
  INTEGER:: NumberOfFields
  INTEGER:: NumberOfDimensions
  INTEGER:: ValueIndex
  INTEGER:: NumberOfVariableComponents
  INTEGER:: NumberOfMeshComponents
  INTEGER:: NumberOfMaterialComponents
  INTEGER:: NumberOfNodesDefined
  INTEGER:: NumberOfFieldComponent(3)
  INTEGER:: NumberOfElements
  INTEGER:: GlobalElementNumber(10)
  INTEGER:: MaxNodesPerElement
  INTEGER:: MaxNodesPerMeshComponent

  INTEGER:: ELEMENT_NUMBER
  DOUBLE PRECISION:: XI_COORDINATES(3)
   DOUBLE PRECISION:: COORDINATES(3), test

  INTEGER:: lagrange_simplex

  INTEGER, DIMENSION(:), ALLOCATABLE:: NodesPerMeshComponent

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeXValue
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeYValue
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeZValue 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeUValue 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeVValue 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeWValue 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodePValue 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: NodeMUValue 


  TYPE(FIELD_TYPE), POINTER :: FIELD
  TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS
  TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT

  DOUBLE PRECISION:: ScaleFactorsPerElementNodes(10,10)

  INTEGER:: TRI_BASIS, TET_BASIS, QUAD_BASIS, HEX_BASIS
  CHARACTER*2 NMs(99),KNOT




  CONTAINS


  !HERE THE TYPES DEFINED ABOVE ARE FILLED WITH THE DATA PROVIDED  
  SUBROUTINE READ_CMGUI_EXE(REGION)

  INTEGER:: I,J,K,L,M,N

  INTEGER :: ERR !<The error code
  TYPE(VARYING_STRING):: ERROR !<The error string

   TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to get the coordinate system for

       KNOT = '0'
       NMs(1) = '1'
       NMs(2) = '2'
       NMs(3) = '3'
       NMs(4) = '4'
       NMs(5) = '5'
       NMs(6) = '6'
       NMs(7) = '7'
       NMs(8) = '8'
       NMs(9) = '9'
       K = 9
       DO I = 1,9
          K = K + 1
          NMs(K) = TRIM(NMs(I))//TRIM(KNOT)
          DO J = 1,9
             K = K + 1
             NMs(K) = TRIM(NMs(I))//TRIM(NMs(J))
          END DO
       END DO




   NumberOfFields=REGION%fields%number_of_fields
   NumberOfDimensions=REGION%coordinate_system%number_of_dimensions

  NumberOfVariableComponents=REGION%equations_sets%equations_sets(1)%ptr%dependent%dependent_field%&
  &variables(1)%number_of_components

   NumberOfMaterialComponents=REGION%equations_sets%equations_sets(1)%ptr%materials%materials_field%&
   &variables(1)%number_of_components

   NumberOfElements=REGION%meshes%meshes(1)%ptr%number_of_elements
   NumberOfMeshComponents=REGION%meshes%meshes(1)%ptr%number_of_components
   ALLOCATE(NodesPerElement(NumberOfMeshComponents))
   ALLOCATE(NodesPerMeshComponent(NumberOfMeshComponents))
   MaxNodesPerElement=0
   DO I=1,NumberOfMeshComponents
      NodesPerElement(I)=REGION%fields%fields(1)%ptr%geometric_field%decomposition%domain(1)&
      &%ptr%topology%elements%elements(1)%basis%number_of_element_parameters
      NodesPerMeshComponent(I)=REGION%meshes%meshes(1)%ptr%topology(I)%ptr%nodes%number_of_nodes
   END DO


   MaxNodesPerElement=NodesPerElement(1)



   MaxNodesPerMeshComponent=NodesPerMeshComponent(1)



   ALLOCATE(NodeXValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeYValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeZValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeUValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeVValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeWValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodePValue(NodesPerMeshComponent(1)))
   ALLOCATE(NodeMUValue(NodesPerMeshComponent(1)))
   ALLOCATE(ElementNodesScales(NumberOfElements,NodesPerElement(1)))
   ALLOCATE(ElementNodes(NumberOfElements,NodesPerElement(1)))


! THIS NEEDS TO BE ADJUSTED NOW!!!!!!

    CALL ENTERS("CMGUI OUTPUT",ERR,ERROR,*999)

   FIELD=>REGION%equations_sets%equations_sets(1)%ptr%dependent%dependent_field
   CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(FIELD,FIELD_U_VARIABLE_TYPE,INTERPOLATION_PARAMETERS&
   &,ERR,ERROR,*999)
   CALL FIELD_INTERPOLATED_POINT_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,ERR,ERROR,*999)


  DO I=1,NumberOfElements
   DO J=1,NodesPerElement(1)
 
      ELEMENT_NUMBER=I
      XI_COORDINATES(1)=(REGION%equations_sets%equations_sets(1)%ptr%equations%interpolation%&
      &geometric_interp_parameters%bases(1)%ptr%node_position_index(J,1)-1.0)/(REGION%equations_sets%&
      &equations_sets(1)%ptr%equations%interpolation%geometric_interp_parameters%bases(1)&
      &%ptr%number_of_nodes_xi(1)-1.0)
      XI_COORDINATES(2)=(REGION%equations_sets%equations_sets(1)%ptr%equations%interpolation%&
      &geometric_interp_parameters%bases(1)%ptr%node_position_index(J,2)-1.0)/(REGION%equations_sets%&
      &equations_sets(1)%ptr%equations%interpolation%geometric_interp_parameters%bases(1)&
      &%ptr%number_of_nodes_xi(2)-1.0)
      XI_COORDINATES(3)=(REGION%equations_sets%equations_sets(1)%ptr%equations%interpolation%&
      &geometric_interp_parameters%bases(1)%ptr%node_position_index(J,3)-1.0)/(REGION%equations_sets%&
      &equations_sets(1)%ptr%equations%interpolation%geometric_interp_parameters%bases(1)&
      &%ptr%number_of_nodes_xi(3)-1.0)

      !K is global node number
      K=REGION%meshes%meshes(1)%ptr%topology(1)%ptr%elements%elements(I)%global_element_nodes(J)


      COORDINATES=(/1,1,1/)

      CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,&
      &INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
 
      CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI_COORDINATES,INTERPOLATED_POINT,ERR,ERROR,*999)

      NodeXValue(K)=REGION%equations_sets%equations_sets(1)%ptr%geometry%geometric_field%variables(1)&
      &%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K)

      NodeYValue(K)=REGION%equations_sets%equations_sets(1)%ptr%geometry%geometric_field%variables(1)&
      &%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+NodesPerMeshComponent(1))

      NodeZValue(K)=REGION%equations_sets%equations_sets(1)%ptr%geometry%geometric_field%variables(1)&
      &%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(K+2*NodesPerMeshComponent(1))


      NodeUValue(K)=INTERPOLATED_POINT%VALUES(1,1)
      NodeVValue(K)=INTERPOLATED_POINT%VALUES(2,1)
      NodeWValue(K)=INTERPOLATED_POINT%VALUES(3,1)
      NodePValue(K)=INTERPOLATED_POINT%VALUES(4,1)




   END DO 
  END DO



   NodeMUValue=REGION%equations_sets%equations_sets(1)%ptr%materials%materials_field%variables(1)%&
   &parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(1)

  lagrange_simplex=REGION%equations_sets%equations_sets(1)%ptr%equations%interpolation%geometric_field%type

  NumberOfFieldComponent(1)=NumberOfDimensions
  NumberOfFieldComponent(2)=NumberOfVariableComponents
  NumberOfFieldComponent(3)=NumberOfMaterialComponents


    DO I=1,NumberOfElements
    DO J=1,NodesPerElement(1)
      ElementNodes(I,J)=REGION%meshes%meshes(1)%ptr%topology(1)%&
      &ptr%elements%elements(I)%global_element_nodes(J)
      ElementNodesScales(I,J)=1.0000000000000000E+00
    END DO
    END DO
    CALL EXITS("CMGUI OUTPUT")
    RETURN
999 CALL ERRORS("CMGUI OUTPUT",ERR,ERROR)


  END SUBROUTINE READ_CMGUI_EXE


! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------

  !HERE THE TYPES DEFINED ABOVE ARE FILLED WITH THE DATA PROVIDED
  SUBROUTINE SEND_CMGUI_EXE
     CALL WRITE_NODE_FILE
     CALL WRITE_ELEMENT_FILE
  END SUBROUTINE SEND_CMGUI_EXE

! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------

  SUBROUTINE WRITE_NODE_FILE

  IMPLICIT NONE

  INTEGER:: I,J,K,L,M,N

       OPEN(UNIT=14, FILE='./output/cmgui.exnode',STATUS='unknown')

! WRITING HEADER INFORMATION

       WRITE(14,*) 'Group name: OpenCMISS'
       WRITE(14,*) '#Fields=',TRIM(NMs(NumberOfFields))
       ValueIndex=1
       WRITE(14,*) ' 1) coordinates,  coordinate, rectangular cartesian, #Components=',TRIM(NMs(NumberOfDimensions))
       DO I=1,NumberOfDimensions
         IF(I==1) THEN
           WRITE(14,*) '   x.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0'
         ELSE IF(I==2) THEN
           WRITE(14,*) '   y.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0'
         ELSE
           WRITE(14,*) '   z.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0'
         END IF
         ValueIndex=ValueIndex+1
       END DO
       WRITE(14,*) ' 2) general,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfVariableComponents))
       DO I=1,NumberOfVariableComponents
         WRITE(14,*)  '   ',TRIM(NMs(I)),'.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0' 
         ValueIndex=ValueIndex+1
       END DO
       WRITE(14,*) ' 3) material,  field,  rectangular cartesian, #Components=',TRIM(NMs(NumberOfMaterialComponents))
       DO I=1,NumberOfMaterialComponents
         WRITE(14,*)  '   ',TRIM(NMs(I)),'.  Value index= ',TRIM(NMs(ValueIndex)),',     #Derivatives= 0' 
         ValueIndex=ValueIndex+1
       END DO

! NOW WRITE NODE INFORMATION
! ! ! 
            DO I = 1,NodesPerMeshComponent(1)
               WRITE(14,*) ' Node: ',I



                  WRITE(14,'("    ", es25.16 )')NodeXValue(I)
                  WRITE(14,'("    ", es25.16 )')NodeYValue(I)
                  WRITE(14,'("    ", es25.16 )')NodeZValue(I)
                  WRITE(14,'("    ", es25.16 )')NodeUValue(I)
                  WRITE(14,'("    ", es25.16 )')NodeVValue(I)
                  WRITE(14,'("    ", es25.16 )')NodeWValue(I)
                  WRITE(14,'("    ", es25.16 )')NodePValue(I)
                  WRITE(14,'("    ", es25.16 )')NodeMUValue(I)

            END DO


       WRITE(14,*) ' '
       CLOSE(14)

  WRITE(*,*)'Writing Nodes...'

  END SUBROUTINE WRITE_NODE_FILE


! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------

  SUBROUTINE WRITE_ELEMENT_FILE

   IMPLICIT NONE

        INTEGER:: I,J,K,L,M,N

        CHARACTER*60 ELEM_TYPE


       OPEN(UNIT=5,FILE='./output/cmgui.exelem',STATUS='unknown')
       WRITE(5,*) 'Group name: OpenCMISS'

       IF(lagrange_simplex==2) THEN
         WRITE(5,*) 'Shape.  Dimension=',TRIM(NMs(NumberOfDimensions)),', simplex(2;3)*simplex*simplex'
         IF(MaxNodesPerElement==3) THEN
           WRITE(5,*) '#Scale factor sets= 1'
           WRITE(5,*) ' l.simplex(2)*l.simplex, #Scale factors= ', NodesPerElement(1)
         ELSE IF(MaxNodesPerElement==4) THEN
           WRITE(5,*) '#Scale factor sets= 1'
           WRITE(5,*) ' l.simplex(2;3)*l.simplex*l.simplex, #Scale factors= ', NodesPerElement(1)
         ELSE
           WRITE(5,*) '#Scale factor sets= 0'
         END IF
       ELSE IF (lagrange_simplex==1) THEN
         WRITE(5,*) 'Shape.  Dimension= ',TRIM(NMs(NumberOfDimensions))
         WRITE(5,*) '#Scale factor sets= 1'
         IF(NumberOfDimensions==2) THEN
             IF(MaxNodesPerElement==4) THEN
                   WRITE(5,*) 'l.Lagrange*l.Lagrange, #Scale factors=',NodesPerElement(1)
             ELSE IF(MaxNodesPerElement==9) THEN
                   WRITE(5,*) 'q.Lagrange*q.Lagrange, #Scale factors=',NodesPerElement(1)
             ELSE IF(MaxNodesPerElement==16) THEN
                   WRITE(5,*) 'c.Lagrange*c.Lagrange, #Scale factors=',NodesPerElement(1)
             END IF
         ELSE
             IF(MaxNodesPerElement==8) THEN
                   WRITE(5,*) 'l.Lagrange*l.Lagrange*l.Lagrange, #Scale factors=',NodesPerElement(1)
             ELSE IF(MaxNodesPerElement==27) THEN
                   WRITE(5,*) 'q.Lagrange*q.Lagrange*q.Lagrange, #Scale factors=',NodesPerElement(1)
             ELSE IF(MaxNodesPerElement==64) THEN
                   WRITE(5,*) 'c.Lagrange*c.Lagrange*c.Lagrange, #Scale factors=',NodesPerElement(1)
             END IF
         END IF
      END IF

       WRITE(5,*) '#Nodes= ',TRIM(NMs(NodesPerElement(1)))
       WRITE(5,*) '#Fields= ',TRIM(Nms(NumberOfFields))


       DO I=1,NumberOfFields


          IF(I==1)THEN
           WRITE(5,*)' 1) coordinates,  coordinate, rectangular cartesian, #Components= ',TRIM(NMs(NumberOfDimensions))
          ELSE IF(I==2) THEN
           WRITE(5,*)' 2) general,  field,  rectangular cartesian, #Components= ',TRIM(NMs(NumberOfVariableComponents))
          ELSE IF(I==3) THEN
           WRITE(5,*)' 3) material,  field,  rectangular cartesian, #Components= ',TRIM(NMs(NumberOfMaterialComponents))
          END IF

      DO J=1,NumberOfFieldComponent(I)

        IF(NumberOfDimensions==2) THEN

             IF(I==1)THEN
              IF(J==1) THEN
                IF(MaxNodesPerElement==4)THEN
                    WRITE(5,*)'   x.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==9) THEN
                    WRITE(5,*)'   x.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==16)  THEN
                    WRITE(5,*)'   x.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                END IF 
               ELSE IF(J==2) THEN
                IF(MaxNodesPerElement==4) THEN
                    WRITE(5,*)'   y.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==9)  THEN
                    WRITE(5,*)'   y.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==16)  THEN
                    WRITE(5,*)'   y.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                END IF
               ELSE IF(J==3) THEN
                IF(MaxNodesPerElement==4) THEN
                    WRITE(5,*)'   z.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==9)  THEN
                    WRITE(5,*)'   z.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==16)  THEN
                    WRITE(5,*)'   z.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                END IF
              END IF
             ELSE
                IF(MaxNodesPerElement==4) THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==9)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==16)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                END IF
             END IF


          ELSE IF(NumberOfDimensions==3) THEN

             IF(I==1)THEN
              IF(J==1) THEN
                IF(MaxNodesPerElement==8) THEN
                    WRITE(5,*)'   x.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==27)  THEN
                    WRITE(5,*)'   x.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==64)  THEN
                    WRITE(5,*)'   x.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
                END IF 
               ELSE IF(J==2) THEN
                IF(MaxNodesPerElement==8) THEN
                    WRITE(5,*)'   y.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==27)  THEN
                    WRITE(5,*)'   y.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==64)  THEN
                    WRITE(5,*)'   y.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
                END IF
               ELSE IF(J==3) THEN
                IF(MaxNodesPerElement==8) THEN
                    WRITE(5,*)'   z.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==27)  THEN
                    WRITE(5,*)'   z.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==64)  THEN
                    WRITE(5,*)'   z.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
                END IF
              END IF
             ELSE
                IF(MaxNodesPerElement==8) THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==27)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==64)  THEN
                WRITE(5,*)'   ',TRIM(NMs(J)),'.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
                END IF
             END IF
          END IF



               WRITE(5,*) '   #Nodes= ',TRIM(NMs(MaxNodesPerElement))



           DO K = 1,MaxNodesPerElement
               WRITE(5,*) '    ',TRIM(NMs(K)),'.  #Values=1'
               WRITE(5,*) '     Value indices:     1'
               WRITE(5,*) '     Scale factor indices:   ',TRIM(NMs(K))
           END DO

          END DO

        END DO

       DO K = 1,NumberOfElements
            WRITE(5,*) 'Element:     ', K,' 0  0'
            WRITE(5,*) '   Nodes:'
            WRITE(5,*) '   ', ElementNodes(K,1:NodesPerElement(1))
            WRITE(5,*) '   Scale factors:'
            WRITE(5,*) '   ',ElementNodesScales(K,1:NodesPerElement(1))
       END DO

       WRITE(5,*) ' '
       CLOSE(5)



  WRITE(*,*)'Writing Elements...'

  END SUBROUTINE WRITE_ELEMENT_FILE


! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------

END MODULE EXPORT_CMGUI
