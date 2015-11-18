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
  INTEGER(CMISSIntg),PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg),PARAMETER :: DataProjectionUserNumber=1
  INTEGER(CMISSIntg),PARAMETER :: FieldUserNumber=1  
  INTEGER(CMISSIntg),PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg),PARAMETER :: RegionUserNumber=1

  REAL(CMISSRP), PARAMETER :: CoordinateSystemOrigin(3)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP]
  
  !Program types

  !Program variables   
  INTEGER(CMISSIntg) :: MeshComponentNumber=1
  INTEGER(CMISSIntg) :: MeshDimensions=1
  INTEGER(CMISSIntg) :: MeshNumberOfElements
  INTEGER(CMISSIntg) :: MeshNumberOfComponents=1
  INTEGER(CMISSIntg) :: NumberOfDomains=2
  INTEGER(CMISSIntg) :: NumberOfNodes
  INTEGER(CMISSIntg) :: NumberOfXi=1
  INTEGER(CMISSIntg) :: BasisInterpolation(1)=[CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION]
  INTEGER(CMISSIntg) :: WorldCoordinateSystemUserNumber
  INTEGER(CMISSIntg) :: WorldRegionUserNumber
  
  INTEGER(CMISSIntg) :: FieldNumberOfVariables=1
  INTEGER(CMISSIntg) :: FieldNumberOfComponents=3
  

  INTEGER(CMISSIntg) :: np,el,xi,ver_idx,der_idx,node_idx,comp_idx
    
  REAL(CMISSRP), DIMENSION(5,3) :: DataPointValues !(number_of_data_points,dimension)
  INTEGER(CMISSIntg), DIMENSION(5,2) :: ElementUserNodes  
  REAL(CMISSRP), DIMENSION(2,6,3) :: FieldValues
  
  !Test variables
  REAL(CMISSRP) :: AbsoluteToleranceSet=1.0E-10_CMISSRP !default is 1.0E-8
  REAL(CMISSRP) :: RelativeToleranceSet=1.0E-6_CMISSRP !default is 1.0E-8
  INTEGER(CMISSIntg) :: MaximumNumberOfIterationsSet=30 !default is 25
  REAL(CMISSRP) :: MaximumIterationUpdateSet=0.4_CMISSRP !default is 0.5
  INTEGER(CMISSIntg) :: NumberOfClosestElementsSet=3 !default is 2/4/8 for 1/2/3 dimensional projection 
  INTEGER(CMISSIntg) :: ProjectionTypeSet=CMFE_DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE !same as default
  REAL(CMISSRP) :: StartingXiSet(1)=[0.4_CMISSRP] !default is 0.5
  REAL(CMISSRP) :: AbsoluteToleranceGet
  REAL(CMISSRP) :: RelativeToleranceGet
  INTEGER(CMISSIntg) :: MaximumNumberOfIterationsGet
  REAL(CMISSRP) :: MaximumIterationUpdateGet
  INTEGER(CMISSIntg) :: NumberOfClosestElementsGet
  INTEGER(CMISSIntg) :: ProjectionTypeGet
  REAL(CMISSRP), ALLOCATABLE :: StartingXiGet(:)
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
  DataPointValues(1,:)=[20.5_CMISSRP,1.8_CMISSRP,0.0_CMISSRP]
  DataPointValues(2,:)=[33.2_CMISSRP,-4.8_CMISSRP,0.0_CMISSRP]  
  DataPointValues(3,:)=[9.6_CMISSRP,10.0_CMISSRP,0.0_CMISSRP]  
  DataPointValues(4,:)=[50.0_CMISSRP,-3.0_CMISSRP,6.0_CMISSRP]  
  DataPointValues(5,:)=[44.0_CMISSRP,10.0_CMISSRP,18.6_CMISSRP]  
  
  ElementUserNodes(1,:)=[1,2]
  ElementUserNodes(2,:)=[2,3]
  ElementUserNodes(3,:)=[3,4]
  ElementUserNodes(4,:)=[4,5]
  ElementUserNodes(5,:)=[5,6]        
  
  FieldValues(1,1,:)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !no der, node 1
  FieldValues(2,1,:)=[10.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !first der, node 1
  
  FieldValues(1,2,:)=[10.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !no der, node 2
  FieldValues(2,2,:)=[10.0_CMISSRP,-10.0_CMISSRP,0.0_CMISSRP] !first der, node 2
  
  FieldValues(1,3,:)=[20.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !no der, node 3
  FieldValues(2,3,:)=[10.0_CMISSRP,20.0_CMISSRP,0.0_CMISSRP] !first der, node 3
  
  FieldValues(1,4,:)=[30.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !no der, node 4
  FieldValues(2,4,:)=[10.0_CMISSRP,10.0_CMISSRP,0.0_CMISSRP] !first der, node 4
  
  FieldValues(1,5,:)=[40.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !no der, node 5
  FieldValues(2,5,:)=[10.0_CMISSRP,-15.0_CMISSRP,0.0_CMISSRP] !first der, node 5
  
  FieldValues(1,6,:)=[50.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !no der, node 6
  FieldValues(2,6,:)=[10.0_CMISSRP,-5.0_CMISSRP,0.0_CMISSRP] !first der, node 6
  
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
  DO np=1,SIZE(DataPointValues,1)
    CALL cmfe_DataPoints_ValuesSet(RegionUserNumber,np,DataPointValues(np,:),Err)     
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
  MeshNumberOfElements=SIZE(ElementUserNodes,1)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,RegionUserNumber,MeshDimensions,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(RegionUserNumber,MeshUserNumber,MeshNumberOfComponents,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(RegionUserNumber,MeshUserNumber,MeshNumberOfElements,Err)
  !define nodes for the mesh
  NumberOfNodes=SIZE(FieldValues,2)
  CALL cmfe_Nodes_CreateStart(RegionUserNumber,NumberOfNodes,Err)
  CALL cmfe_Nodes_CreateFinish(RegionUserNumber,Err)  
  !define elements for the mesh
  CALL cmfe_MeshElements_CreateStart(RegionUserNumber,MeshUserNumber,MeshComponentNumber,BasisUserNumber,Err)
  Do el=1,MeshNumberOfElements
    CALL cmfe_MeshElements_NodesSet(RegionUserNumber,MeshUserNumber,MeshComponentNumber,el,ElementUserNodes(el,:),Err)
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
  DO xi=1,NumberOfXi
    CALL cmfe_Field_ComponentMeshComponentSet(RegionUserNumber,FieldUserNumber,CMFE_FIELD_U_VARIABLE_TYPE,xi,xi,Err)
  ENDDO !xi    
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
  !=========================================================================================================================
  !Test parameter set functions
  CALL cmfe_DataProjection_AbsoluteToleranceSet(DataProjectionUserNumber,RegionUserNumber,AbsoluteToleranceSet,Err)
  CALL cmfe_DataProjection_MaximumIterationUpdateSet(DataProjectionUserNumber,RegionUserNumber,MaximumIterationUpdateSet,Err)
  CALL cmfe_DataProjection_MaximumNumberOfIterationsSet(DataProjectionUserNumber,RegionUserNumber,MaximumNumberOfIterationsSet,Err)
  CALL cmfe_DataProjection_NumberOfClosestElementsSet(DataProjectionUserNumber,RegionUserNumber,NumberOfClosestElementsSet,Err)
  CALL cmfe_DataProjection_ProjectionTypeSet(DataProjectionUserNumber,RegionUserNumber,ProjectionTypeSet,Err)
  CALL cmfe_DataProjection_RelativeToleranceSet(DataProjectionUserNumber,RegionUserNumber,RelativeToleranceSet,Err)
  CALL cmfe_DataProjection_StartingXiSet(DataProjectionUserNumber,RegionUserNumber,StartingXiSet,Err)
  !=========================================================================================================================
  !Finish data projection  
  CALL cmfe_DataProjection_CreateFinish(DataProjectionUserNumber,RegionUserNumber,Err)
  !=========================================================================================================================
  !Test parameter get functions
  CALL cmfe_DataProjection_AbsoluteToleranceGet(DataProjectionUserNumber,RegionUserNumber,AbsoluteToleranceGet,Err) 
  CALL cmfe_DataProjection_MaximumIterationUpdateGet(DataProjectionUserNumber,RegionUserNumber,MaximumIterationUpdateGet,Err)
  CALL cmfe_DataProjection_MaximumNumberOfIterationsGet(DataProjectionUserNumber,RegionUserNumber,MaximumNumberOfIterationsGet,Err)
  CALL cmfe_DataProjection_NumberOfClosestElementsGet(DataProjectionUserNumber,RegionUserNumber,NumberOfClosestElementsGet,Err) 
  CALL cmfe_DataProjection_ProjectionTypeGet(DataProjectionUserNumber,RegionUserNumber,ProjectionTypeGet,Err)
  CALL cmfe_DataProjection_RelativeToleranceGet(DataProjectionUserNumber,RegionUserNumber,RelativeToleranceGet,Err)
  CALL cmfe_DataProjection_StartingXiGet(DataProjectionUserNumber,RegionUserNumber,StartingXiGet,Err)
  
  !=========================================================================================================================
  !Start data projection
  CALL cmfe_DataProjection_DataPointsProjectionEvaluate(DataProjectionUserNumber,RegionUserNumber,FieldUserNumber, &
    & RegionUserNumber,Err)
  
  !=========================================================================================================================
  !Destroy used types
  CALL cmfe_DataProjection_Destroy(DataProjectionUserNumber,RegionUserNumber,Err)
  CALL cmfe_DataPoints_Destroy(RegionUserNumber,Err)
    
  CALL cmfe_Region_Destroy(RegionUserNumber,Err)
  CALL cmfe_CoordinateSystem_Destroy(CoordinateSystemUserNumber,Err)  
  
  !=========================================================================================================================
  !Finishing program
  CALL cmfe_Finalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  STOP  
  
END PROGRAM DataProjection1DRectangularCartesian
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
