!> \file
!> \author Tim Wu
!> \brief This is an example program to solve 3D data points embedding in 3D cartesian elements using OpenCMISS calls.
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
PROGRAM DataProjection3DRectangularCartesian
#ifndef NOMPIMOD
  USE MPI
#endif

  USE OpenCMISS
  USE OpenCMISS_Iron

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Program parameters
  INTEGER(CMISSIntg),PARAMETER :: BasisUserNumber=1  
  INTEGER(CMISSIntg),PARAMETER :: CoordinateSystemDimension=3
  INTEGER(CMISSIntg),PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg),PARAMETER :: DataProjectionUserNumber=1
  INTEGER(CMISSIntg),PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg),PARAMETER :: FieldUserNumber=1  
  INTEGER(CMISSIntg),PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg),PARAMETER :: RegionUserNumber=1
  REAL(CMISSRP), PARAMETER :: CoordinateSystemOrigin(3)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP]  

  !Program variables
  INTEGER(CMISSIntg) :: MeshComponentNumber=1
  INTEGER(CMISSIntg) :: NumberOfDataPoints 
  INTEGER(CMISSIntg) :: MeshDimensions=3
  INTEGER(CMISSIntg) :: MeshNumberOfElements
  INTEGER(CMISSIntg) :: MeshNumberOfComponents=1
  INTEGER(CMISSIntg) :: NumberOfDomains=1 !NumberOfDomains=2 for parallel processing, need to set up MPI
  INTEGER(CMISSIntg) :: NumberOfNodes
  INTEGER(CMISSIntg) :: NumberOfXi=3
  INTEGER(CMISSIntg) :: BasisInterpolation(3)=[CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION, &
    & CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION]
  INTEGER(CMISSIntg) :: WorldCoordinateSystemUserNumber
  INTEGER(CMISSIntg) :: WorldRegionUserNumber
  
  INTEGER(CMISSIntg) :: FieldNumberOfVariables=1
  INTEGER(CMISSIntg) :: FieldNumberOfComponents=3 
  
  INTEGER(CMISSIntg) :: data_point_idx,elem_idx,ver_idx,der_idx,node_idx,comp_idx
  REAL(CMISSRP), DIMENSION(5,3) :: DataPointValues !(number_of_data_points,dimension)
  REAL(CMISSRP), DIMENSION(5) :: DataPointProjectionDistance !(number_of_data_points)
  INTEGER(CMISSIntg), DIMENSION(5) :: DataPointProjectionElementNumber !(number_of_data_points)
  INTEGER(CMISSIntg), DIMENSION(5) :: DataPointProjectionExitTag !(number_of_data_points)
  REAL(CMISSRP), DIMENSION(5,3) :: DataPointProjectionXi !(number_of_data_points,MeshDimensions)  
  INTEGER(CMISSIntg), DIMENSION(1,8) :: ElementUserNodes  
  REAL(CMISSRP), DIMENSION(8,8,3) :: FieldValues
        
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
    
  !Intialise 5 data points
  DataPointValues(1,:)=[0.1_CMISSRP,0.8_CMISSRP,1.0_CMISSRP]
  DataPointValues(2,:)=[0.5_CMISSRP,0.5_CMISSRP,0.5_CMISSRP]  
  DataPointValues(3,:)=[0.2_CMISSRP,0.5_CMISSRP,0.5_CMISSRP]  
  DataPointValues(4,:)=[0.9_CMISSRP,0.6_CMISSRP,0.9_CMISSRP]
  DataPointValues(5,:)=[0.3_CMISSRP,0.3_CMISSRP,0.3_CMISSRP]
  NumberOfDataPoints=SIZE(DataPointValues,1)
  
  !Intialise 1 element
  ElementUserNodes(1,:)=[1,2,3,4,5,6,7,8]     
  MeshNumberOfElements=SIZE(ElementUserNodes,1)
  
  !Intialise 8 nodes for the element
  FieldValues(1,1,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !no der, node 1
  FieldValues(2,1,:)=[1.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 der, node 1
  FieldValues(3,1,:)=[0.0_CMISSRP,1.0_CMISSRP,0.0_CMISSRP] !s2 der, node 1  
  FieldValues(4,1,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 der, node 1    
  FieldValues(5,1,:)=[0.0_CMISSRP,0.0_CMISSRP,1.0_CMISSRP] !s3 der, node 1
  FieldValues(6,1,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s3 der, node 1  
  FieldValues(7,1,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s2 s3 der, node 1    
  FieldValues(8,1,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 s3 der, node 1    
  
  FieldValues(1,2,:)=[1.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !no der, node 2
  FieldValues(2,2,:)=[1.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 der, node 2
  FieldValues(3,2,:)=[0.0_CMISSRP,1.0_CMISSRP,0.0_CMISSRP] !s2 der, node 2  
  FieldValues(4,2,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 der, node 2    
  FieldValues(5,2,:)=[0.0_CMISSRP,0.0_CMISSRP,1.0_CMISSRP] !s3 der, node 2
  FieldValues(6,2,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s3 der, node 2  
  FieldValues(7,2,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s2 s3 der, node 2    
  FieldValues(8,2,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 s3 der, node 2     
  
  FieldValues(1,3,:)=[0.0_CMISSRP,1.0_CMISSRP,0.0_CMISSRP] !no der, node 3
  FieldValues(2,3,:)=[1.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 der, node 3
  FieldValues(3,3,:)=[0.0_CMISSRP,1.0_CMISSRP,0.0_CMISSRP] !s2 der, node 3  
  FieldValues(4,3,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 der, node 3    
  FieldValues(5,3,:)=[0.0_CMISSRP,0.0_CMISSRP,1.0_CMISSRP] !s3 der, node 3
  FieldValues(6,3,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s3 der, node 3  
  FieldValues(7,3,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s2 s3 der, node 3    
  FieldValues(8,3,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 s3 der, node 3
  
  FieldValues(1,4,:)=[1.0_CMISSRP,1.0_CMISSRP,0.0_CMISSRP] !no der, node 4
  FieldValues(2,4,:)=[1.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 der, node 4
  FieldValues(3,4,:)=[0.0_CMISSRP,1.0_CMISSRP,0.0_CMISSRP] !s2 der, node 4  
  FieldValues(4,4,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 der, node 4    
  FieldValues(5,4,:)=[0.0_CMISSRP,0.0_CMISSRP,1.0_CMISSRP] !s3 der, node 4
  FieldValues(6,4,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s3 der, node 4  
  FieldValues(7,4,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s2 s3 der, node 4    
  FieldValues(8,4,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 s3 der, node 4
  
  FieldValues(1,5,:)=[0.0_CMISSRP,0.0_CMISSRP,1.0_CMISSRP] !no der, node 5
  FieldValues(2,5,:)=[1.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 der, node 5
  FieldValues(3,5,:)=[0.0_CMISSRP,1.0_CMISSRP,0.0_CMISSRP] !s2 der, node 5  
  FieldValues(4,5,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 der, node 5    
  FieldValues(5,5,:)=[0.0_CMISSRP,0.0_CMISSRP,1.0_CMISSRP] !s3 der, node 5
  FieldValues(6,5,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s3 der, node 5  
  FieldValues(7,5,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s2 s3 der, node 5    
  FieldValues(8,5,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 s3 der, node 5    
  
  FieldValues(1,6,:)=[1.0_CMISSRP,0.0_CMISSRP,1.0_CMISSRP] !no der, node 6
  FieldValues(2,6,:)=[1.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 der, node 6
  FieldValues(3,6,:)=[0.0_CMISSRP,1.0_CMISSRP,0.0_CMISSRP] !s2 der, node 6  
  FieldValues(4,6,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 der, node 6    
  FieldValues(5,6,:)=[0.0_CMISSRP,0.0_CMISSRP,1.0_CMISSRP] !s3 der, node 6
  FieldValues(6,6,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s3 der, node 6  
  FieldValues(7,6,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s2 s3 der, node 6    
  FieldValues(8,6,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 s3 der, node 6     
  
  FieldValues(1,7,:)=[0.0_CMISSRP,1.0_CMISSRP,1.0_CMISSRP] !no der, node 7
  FieldValues(2,7,:)=[1.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 der, node 7
  FieldValues(3,7,:)=[0.0_CMISSRP,1.0_CMISSRP,0.0_CMISSRP] !s2 der, node 7  
  FieldValues(4,7,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 der, node 7    
  FieldValues(5,7,:)=[0.0_CMISSRP,0.0_CMISSRP,1.0_CMISSRP] !s3 der, node 7
  FieldValues(6,7,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s3 der, node 7  
  FieldValues(7,7,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s2 s3 der, node 7    
  FieldValues(8,7,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 s3 der, node 7
  
  FieldValues(1,8,:)=[1.0_CMISSRP,1.0_CMISSRP,1.0_CMISSRP] !no der, node 8
  FieldValues(2,8,:)=[1.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 der, node 8
  FieldValues(3,8,:)=[0.0_CMISSRP,1.0_CMISSRP,0.0_CMISSRP] !s2 der, node 8  
  FieldValues(4,8,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 der, node 8    
  FieldValues(5,8,:)=[0.0_CMISSRP,0.0_CMISSRP,1.0_CMISSRP] !s3 der, node 8
  FieldValues(6,8,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s3 der, node 8  
  FieldValues(7,8,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s2 s3 der, node 8    
  FieldValues(8,8,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !s1 s2 s3 der, node 8
  NumberOfNodes=SIZE(FieldValues,2)
  
  !Intialise cmiss
  CALL cmfe_Initialise(WorldCoordinateSystemUserNumber,WorldRegionUserNumber,Err)
  !Broadcast the number of Elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !=========================================================================================================================
  !Create RC coordinate system
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,Err)
  CALL cmfe_CoordinateSystem_TypeSet(CoordinateSystemUserNumber,CMFE_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystemUserNumber,CoordinateSystemDimension,Err)
  CALL cmfe_CoordinateSystem_OriginSet(CoordinateSystemUserNumber,CoordinateSystemOrigin,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystemUserNumber,Err) 

  !=========================================================================================================================
  !Create Region and set CS to newly created 3D RC CS
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegionUserNumber,Err)
  CALL cmfe_Region_CoordinateSystemSet(RegionUserNumber,CoordinateSystemUserNumber,Err)
  CALL cmfe_Region_CreateFinish(RegionUserNumber,Err)
    
  !=========================================================================================================================
  !Create Data Points and set the values
  CALL cmfe_DataPoints_CreateStart(RegionUserNumber,SIZE(DataPointValues,1),Err)
  DO data_point_idx=1,NumberOfDataPoints
    CALL cmfe_DataPoints_ValuesSet(RegionUserNumber,data_point_idx,DataPointValues(data_point_idx,:),Err)     
  ENDDO
  CALL cmfe_DataPoints_CreateFinish(RegionUserNumber,Err)  
  !=========================================================================================================================
  !Define basis function - 1D cubic hermite
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Err)
  CALL cmfe_Basis_TypeSet(BasisUserNumber,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(BasisUserNumber,NumberOfXi,Err)
  CALL cmfe_Basis_InterpolationXiSet(BasisUserNumber,BasisInterpolation,Err)
  CALL cmfe_Basis_CreateFinish(BasisUserNumber,Err)  
  !=========================================================================================================================
  !Create a mesh
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,RegionUserNumber,MeshDimensions,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(RegionUserNumber,MeshUserNumber,MeshNumberOfComponents,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(RegionUserNumber,MeshUserNumber,MeshNumberOfElements,Err)
  !define nodes for the mesh
  CALL cmfe_Nodes_CreateStart(RegionUserNumber,NumberOfNodes,Err)
  CALL cmfe_Nodes_CreateFinish(RegionUserNumber,Err)  
  !define elements for the mesh
  CALL cmfe_MeshElements_CreateStart(RegionUserNumber,MeshUserNumber,MeshComponentNumber,BasisUserNumber,Err)
  Do elem_idx=1,MeshNumberOfElements
    CALL cmfe_MeshElements_NodesSet(RegionUserNumber,MeshUserNumber,MeshComponentNumber,elem_idx,ElementUserNodes(elem_idx,:),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(RegionUserNumber,MeshUserNumber,MeshComponentNumber,Err)
  CALL cmfe_Mesh_CreateFinish(RegionUserNumber,MeshUserNumber,Err)
  !=========================================================================================================================
  !Create a mesh decomposition 
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,RegionUserNumber,MeshUserNumber,Err)
  CALL cmfe_Decomposition_TypeSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,Err)
  
  !=========================================================================================================================
  !Create a field to put the geometry
  CALL cmfe_Field_CreateStart(FieldUserNumber,RegionUserNumber,Err)
  CALL cmfe_Field_MeshDecompositionSet(RegionUserNumber,FieldUserNumber,MeshUserNumber,DecompositionUserNumber,Err)
  CALL cmfe_Field_TypeSet(RegionUserNumber,FieldUserNumber,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(RegionUserNumber,FieldUserNumber,FieldNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(RegionUserNumber,FieldUserNumber,CMFE_FIELD_U_VARIABLE_TYPE,FieldNumberOfComponents,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(RegionUserNumber,FieldUserNumber,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
!  DO xi=1,NumberOfXi
    
!  ENDDO !xi    
  CALL cmfe_Field_CreateFinish(RegionUserNumber,FieldUserNumber,Err)
  !node 1
  ver_idx=1 ! version number
  DO der_idx=1,SIZE(FieldValues,1)
    DO node_idx=1,SIZE(FieldValues,2)
      DO comp_idx=1,SIZE(FieldValues,3)
        CALL cmfe_Field_ParameterSetUpdateNode(RegionUserNumber,FieldUserNumber,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE, &
          & ver_idx,der_idx,node_idx,comp_idx,FieldValues(der_idx,node_idx,comp_idx),Err)
      ENDDO
    ENDDO
  ENDDO
  
  !=========================================================================================================================
  !Create a data projection
  CALL cmfe_DataProjection_CreateStart(DataProjectionUserNumber,RegionUserNumber,MeshUserNumber,RegionUserNumber,Err)
  CALL cmfe_DataProjection_ProjectionTypeSet(DataProjectionUserNumber,RegionUserNumber, &
    & CMFE_DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE,Err) !Set to element projection for data points embedding. The default is boundary/surface projection.
  CALL cmfe_DataProjection_CreateFinish(DataProjectionUserNumber,RegionUserNumber,Err)
  
  !=========================================================================================================================
  !Start data projection
  CALL cmfe_DataProjection_DataPointsProjectionEvaluate(dataProjectionUserNumber,RegionUserNumber,FieldUserNumber, &
    & RegionUserNumber,Err)

  !Retrieve projection results
  DO data_point_idx=1,NumberOfDataPoints
    CALL cmfe_DataProjection_ResultDistanceGet(RegionUserNumber,DataProjectionUserNumber,data_point_idx, &
      & DataPointProjectionDistance(data_point_idx),Err)
    CALL cmfe_DataProjection_ResultElementNumberGet(RegionUserNumber,DataProjectionUserNumber,data_point_idx, &
      & DataPointProjectionElementNumber(data_point_idx),Err)
    CALL cmfe_DataProjection_ResultExitTagGet(RegionUserNumber,DataProjectionUserNumber,data_point_idx, &
      & DataPointProjectionExitTag(data_point_idx),Err)
    CALL cmfe_DataProjection_ResultXiGet(RegionUserNumber,DataProjectionUserNumber,data_point_idx, &
      & DataPointProjectionXi(data_point_idx,:),Err)
  ENDDO  
  
  !=========================================================================================================================
  !Finishing program
  CALL cmfe_Finalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  STOP  
  
END PROGRAM DataProjection3DRectangularCartesian
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
