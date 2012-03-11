!> \file
!> \author Tim Wu
!> \brief This is an example program to solve 3D data points projecting onto 1D cartesian elements using OpenCMISS calls.
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
!> Contributor(s): Code based on the examples by Kumar Mithraratne and Prasad Babarenda Gamage
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

!> Main program
PROGRAM DataProjection1DRectangularCartesian

  USE MPI
  USE OPENCMISS

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Program parameters
  INTEGER(CMISSIntg),PARAMETER :: BasisUserNumber=1  
  INTEGER(CMISSIntg),PARAMETER :: CoordinateSystemDimension=3
  INTEGER(CMISSIntg),PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg),PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg),PARAMETER :: FieldUserNumber=1  
  INTEGER(CMISSIntg),PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg),PARAMETER :: RegionUserNumber=1

  REAL(CMISSDP), PARAMETER :: CoordinateSystemOrigin(3)=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
  !Program types

  !Program variables   
  INTEGER(CMISSIntg) :: MeshComponentNumber=1
  INTEGER(CMISSIntg) :: MeshDimensions=1
  INTEGER(CMISSIntg) :: MeshNumberOfElements
  INTEGER(CMISSIntg) :: MeshNumberOfComponents=1
  INTEGER(CMISSIntg) :: NumberOfDomains=2
  INTEGER(CMISSIntg) :: NumberOfNodes
  INTEGER(CMISSIntg) :: NumberOfXi=1
  INTEGER(CMISSIntg) :: BasisInterpolation(1)=(/CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION/)
  INTEGER(CMISSIntg) :: WorldCoordinateSystemUserNumber
  INTEGER(CMISSIntg) :: WorldRegionUserNumber
  
  INTEGER(CMISSIntg) :: FieldNumberOfVariables=1
  INTEGER(CMISSIntg) :: FieldNumberOfComponents=3
  

  INTEGER(CMISSIntg) :: np,el,xi,ver_idx,der_idx,node_idx,comp_idx
    
  REAL(CMISSDP), DIMENSION(5,3) :: DataPointValues !(number_of_data_points,dimension)
  INTEGER(CMISSIntg), DIMENSION(5,2) :: ElementUserNodes  
  REAL(CMISSDP), DIMENSION(2,6,3) :: FieldValues
  
  !Test variables
  REAL(CMISSDP) :: AbsoluteToleranceSet=1.0E-10_CMISSDP !default is 1.0E-8
  REAL(CMISSDP) :: RelativeToleranceSet=1.0E-6_CMISSDP !default is 1.0E-8
  INTEGER(CMISSIntg) :: MaximumNumberOfIterationsSet=30 !default is 25
  REAL(CMISSDP) :: MaximumIterationUpdateSet=0.4_CMISSDP !default is 0.5
  INTEGER(CMISSIntg) :: NumberOfClosestElementsSet=3 !default is 2/4/8 for 1/2/3 dimensional projection 
  INTEGER(CMISSIntg) :: ProjectionTypeSet=CMISS_DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE !same as default
  REAL(CMISSDP) :: StartingXiSet(1)=[0.4_CMISSDP] !default is 0.5
  REAL(CMISSDP) :: AbsoluteToleranceGet
  REAL(CMISSDP) :: RelativeToleranceGet
  INTEGER(CMISSIntg) :: MaximumNumberOfIterationsGet
  REAL(CMISSDP) :: MaximumIterationUpdateGet
  INTEGER(CMISSIntg) :: NumberOfClosestElementsGet
  INTEGER(CMISSIntg) :: ProjectionTypeGet
  REAL(CMISSDP), ALLOCATABLE :: StartingXiGet(:)
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS and MPI variables
  INTEGER(CMISSIntg) :: Err
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS=1 !<number of elements on x axis
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_Y_ELEMENTS=1 !<number of elements on y axis
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_Z_ELEMENTS=1 !<number of elements on z axis  
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS=1      
  INTEGER(CMISSIntg) :: MPI_IERROR  
  
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
    
  !Intialise data points
  DataPointValues(1,:)=[20.5_CMISSDP,1.8_CMISSDP,0.0_CMISSDP]
  DataPointValues(2,:)=[33.2_CMISSDP,-4.8_CMISSDP,0.0_CMISSDP]  
  DataPointValues(3,:)=[9.6_CMISSDP,10.0_CMISSDP,0.0_CMISSDP]  
  DataPointValues(4,:)=[50.0_CMISSDP,-3.0_CMISSDP,6.0_CMISSDP]  
  DataPointValues(5,:)=[44.0_CMISSDP,10.0_CMISSDP,18.6_CMISSDP]  
  
  ElementUserNodes(1,:)=[1,2]
  ElementUserNodes(2,:)=[2,3]
  ElementUserNodes(3,:)=[3,4]
  ElementUserNodes(4,:)=[4,5]
  ElementUserNodes(5,:)=[5,6]        
  
  FieldValues(1,1,:)=[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP] !no der, node 1
  FieldValues(2,1,:)=[10.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP] !first der, node 1
  
  FieldValues(1,2,:)=[10.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP] !no der, node 2
  FieldValues(2,2,:)=[10.0_CMISSDP,-10.0_CMISSDP,0.0_CMISSDP] !first der, node 2
  
  FieldValues(1,3,:)=[20.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP] !no der, node 3
  FieldValues(2,3,:)=[10.0_CMISSDP,20.0_CMISSDP,0.0_CMISSDP] !first der, node 3
  
  FieldValues(1,4,:)=[30.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP] !no der, node 4
  FieldValues(2,4,:)=[10.0_CMISSDP,10.0_CMISSDP,0.0_CMISSDP] !first der, node 4
  
  FieldValues(1,5,:)=[40.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP] !no der, node 5
  FieldValues(2,5,:)=[10.0_CMISSDP,-15.0_CMISSDP,0.0_CMISSDP] !first der, node 5
  
  FieldValues(1,6,:)=[50.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP] !no der, node 6
  FieldValues(2,6,:)=[10.0_CMISSDP,-5.0_CMISSDP,0.0_CMISSDP] !first der, node 6
  
  !Intialise cmiss
  CALL CMISSInitialise(WorldCoordinateSystemUserNumber,WorldRegionUserNumber,Err)
  !Broadcast the number of Elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !=========================================================================================================================
  !Create RC coordinate system
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,Err)
  CALL CMISSCoordinateSystem_TypeSet(CoordinateSystemUserNumber,CMISS_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystemUserNumber,CoordinateSystemDimension,Err)
  CALL CMISSCoordinateSystem_OriginSet(CoordinateSystemUserNumber,CoordinateSystemOrigin,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystemUserNumber,Err) 

  !=========================================================================================================================
  !Create Region and set CS to newly created 3D RC CS
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegionUserNumber,Err)
  CALL CMISSRegion_CoordinateSystemSet(RegionUserNumber,CoordinateSystemUserNumber,Err)
  CALL CMISSRegion_CreateFinish(RegionUserNumber,Err)
    
  !=========================================================================================================================
  !Create Data Points and set the values
  CALL CMISSDataPoints_CreateStart(RegionUserNumber,SIZE(DataPointValues,1),Err)
  DO np=1,SIZE(DataPointValues,1)
    CALL CMISSDataPoints_ValuesSet(RegionUserNumber,np,DataPointValues(np,:),Err)     
  ENDDO
  CALL CMISSDataPoints_CreateFinish(RegionUserNumber,Err)  
  !=========================================================================================================================
  !Define basis function - 1D cubic hermite
  CALL CMISSBasis_CreateStart(BasisUserNumber,Err)
  CALL CMISSBasis_TypeSet(BasisUserNumber,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(BasisUserNumber,NumberOfXi,Err)
  CALL CMISSBasis_InterpolationXiSet(BasisUserNumber,BasisInterpolation,Err)
  CALL CMISSBasis_CreateFinish(BasisUserNumber,Err)  
  !=========================================================================================================================
  !Create a mesh
  MeshNumberOfElements=SIZE(ElementUserNodes,1)
  CALL CMISSMesh_CreateStart(MeshUserNumber,RegionUserNumber,MeshDimensions,Err)
  CALL CMISSMesh_NumberOfComponentsSet(RegionUserNumber,MeshUserNumber,MeshNumberOfComponents,Err)
  CALL CMISSMesh_NumberOfElementsSet(RegionUserNumber,MeshUserNumber,MeshNumberOfElements,Err)
  !define nodes for the mesh
  NumberOfNodes=SIZE(FieldValues,2)
  CALL CMISSNodes_CreateStart(RegionUserNumber,NumberOfNodes,Err)
  CALL CMISSNodes_CreateFinish(RegionUserNumber,Err)  
  !define elements for the mesh
  CALL CMISSMeshElements_CreateStart(RegionUserNumber,MeshUserNumber,MeshComponentNumber,BasisUserNumber,Err)
  Do el=1,MeshNumberOfElements
    CALL CMISSMeshElements_NodesSet(RegionUserNumber,MeshUserNumber,MeshComponentNumber,el,ElementUserNodes(el,:),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(RegionUserNumber,MeshUserNumber,MeshComponentNumber,Err)
  CALL CMISSMesh_CreateFinish(RegionUserNumber,MeshUserNumber,Err)
  !=========================================================================================================================
  !Create a mesh decomposition 
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,RegionUserNumber,MeshUserNumber,Err)
  CALL CMISSDecomposition_TypeSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,NumberOfDomains,Err)
  CALL CMISSDecomposition_CreateFinish(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,Err)
  
  !=========================================================================================================================
  !Create a field to put the geometry
  CALL CMISSField_CreateStart(FieldUserNumber,RegionUserNumber,Err)
  CALL CMISSField_MeshDecompositionSet(RegionUserNumber,FieldUserNumber,MeshUserNumber,DecompositionUserNumber,Err)
  CALL CMISSField_TypeSet(RegionUserNumber,FieldUserNumber,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(RegionUserNumber,FieldUserNumber,FieldNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(RegionUserNumber,FieldUserNumber,CMISS_FIELD_U_VARIABLE_TYPE,FieldNumberOfComponents,Err)
  DO xi=1,NumberOfXi
    CALL CMISSField_ComponentMeshComponentSet(RegionUserNumber,FieldUserNumber,CMISS_FIELD_U_VARIABLE_TYPE,xi,xi,Err)
  ENDDO !xi    
  CALL CMISSField_CreateFinish(RegionUserNumber,FieldUserNumber,Err)
  !node 1
  ver_idx=1 ! version number
  DO der_idx=1,SIZE(FieldValues,1)
    DO node_idx=1,SIZE(FieldValues,2)
      DO comp_idx=1,SIZE(FieldValues,3)
        CALL CMISSField_ParameterSetUpdateNode(RegionUserNumber,FieldUserNumber,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE, &
          & ver_idx,der_idx,node_idx,comp_idx,FieldValues(der_idx,node_idx,comp_idx),Err)
      ENDDO
    ENDDO
  ENDDO
  
  !=========================================================================================================================
  !Create a data projection
  CALL CMISSDataProjection_CreateStart(RegionUserNumber,FieldUserNumber,RegionUserNumber,Err)
  !=========================================================================================================================
  !Test parameter set functions
  CALL CMISSDataProjection_AbsoluteToleranceSet(RegionUserNumber,AbsoluteToleranceSet,Err) !test
  CALL CMISSDataProjection_MaximumIterationUpdateSet(RegionUserNumber,MaximumIterationUpdateSet,Err) !test
  CALL CMISSDataProjection_MaximumNumberOfIterationsSet(RegionUserNumber,MaximumNumberOfIterationsSet,Err) !test
  CALL CMISSDataProjection_NumberOfClosestElementsSet(RegionUserNumber,NumberOfClosestElementsSet,Err) !test
  CALL CMISSDataProjection_ProjectionTypeSet(RegionUserNumber,ProjectionTypeSet,Err)
  CALL CMISSDataProjection_RelativeToleranceSet(RegionUserNumber,RelativeToleranceSet,Err) !test
  CALL CMISSDataProjection_StartingXiSet(RegionUserNumber,StartingXiSet,Err) !test
  !=========================================================================================================================
  !Finish data projection  
  CALL CMISSDataProjection_CreateFinish(RegionUserNumber,Err)
  !=========================================================================================================================
  !Test parameter get functions
  CALL CMISSDataProjection_AbsoluteToleranceGet(RegionUserNumber,AbsoluteToleranceGet,Err) !test
  CALL CMISSDataProjection_MaximumIterationUpdateGet(RegionUserNumber,MaximumIterationUpdateGet,Err) !test
  CALL CMISSDataProjection_MaximumNumberOfIterationsGet(RegionUserNumber,MaximumNumberOfIterationsGet,Err) !test
  CALL CMISSDataProjection_NumberOfClosestElementsGet(RegionUserNumber,NumberOfClosestElementsGet,Err) !test
  CALL CMISSDataProjection_ProjectionTypeGet(RegionUserNumber,ProjectionTypeGet,Err) !test
  CALL CMISSDataProjection_RelativeToleranceGet(RegionUserNumber,RelativeToleranceGet,Err) !test
  CALL CMISSDataProjection_StartingXiGet(RegionUserNumber,StartingXiGet,Err) !test  
  
  !=========================================================================================================================
  !Start data projection
  CALL CMISSDataProjection_Evaluate(RegionUserNumber,Err)
  
  !=========================================================================================================================
  !Destroy used types
  CALL CMISSDataProjection_Destroy(RegionUserNumber,Err)
  CALL CMISSDataPoints_Destroy(RegionUserNumber,Err)
    
  CALL CMISSRegion_Destroy(RegionUserNumber,Err)
  CALL CMISSCoordinateSystem_Destroy(CoordinateSystemUserNumber,Err)  
  
  !=========================================================================================================================
  !Finishing program
  CALL CMISSFinalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  STOP  
  
END PROGRAM DataProjection1DRectangularCartesian
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
