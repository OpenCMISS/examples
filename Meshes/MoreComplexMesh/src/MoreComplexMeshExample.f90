!> \file
!> \author Chris Bradley
!> \brief This is an example program which sets up a field which uses a more complex mesh using OpenCMISS calls.
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

!> \example MoreComplexMesh/src/MoreComplexMeshExample.f90
!! Example program which sets up a field which uses a more complex mesh using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MoreComplexMesh/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MoreComplexMesh/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM MORECOMPLEXMESHEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: Basis1UserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: Basis2UserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
 
  !Program types
  
  !Program variables
  
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS

  LOGICAL :: EXPORT_FIELD=.TRUE.

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis1,Basis2
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSMeshElementsType) :: MeshElements1,MeshElements2
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSRegionType) :: Region,WorldRegion

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: Err
  
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

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
  
  NUMBER_OF_DOMAINS=1

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 2D
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Start the creation of a bilinear-Lagrange basis
  CALL CMISSBasis_Initialise(Basis1,Err)
  CALL CMISSBasis_CreateStart(Basis1UserNumber,Basis1,Err)
  !Set the basis to be a bilinear Lagrange basis
  CALL CMISSBasis_NumberOfXiSet(Basis1,2,Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(Basis1,Err)
  
  !Start the creation of a bicubic-Hermite basis
  CALL CMISSBasis_Initialise(Basis2,Err)
  CALL CMISSBasis_CreateStart(Basis2UserNumber,Basis2,Err)
  !Set the basis to be a bilinear Lagrange basis
  CALL CMISSBasis_NumberOfXiSet(Basis2,2,Err)
  !Set the basis to be a bicubic Hermite basis
  CALL CMISSBasis_InterpolationXiSet(Basis2,[CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION,CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION],Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(Basis2,Err)

  !Create a mesh
  !The mesh will consist of two elements. There will be two mesh components. For the first component the first element is
  !bilinear Lagrange and the second element bicubic Hermite. For the second component the first element is bicubic Hermite and
  !the second element is bilinear Lagrange.
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(Region,6,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)
  
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,2,Mesh,Err)
  CALL CMISSMesh_NumberOfElementsSet(Mesh,2,Err)
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,2,Err)
  
  CALL CMISSMeshElements_Initialise(MeshElements1,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,1,Basis2,MeshElements1,Err)
  CALL CMISSMeshElements_NodesSet(MeshElements1,1,[1,2,4,5],Err)
  CALL CMISSMeshElements_BasisSet(MeshElements1,2,Basis1,Err)
  CALL CMISSMeshElements_NodesSet(MeshElements1,2,[2,3,5,6],Err)
  CALL CMISSMeshElements_CreateFinish(MeshElements1,Err)

  CALL CMISSMeshElements_Initialise(MeshElements2,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,2,Basis1,MeshElements2,Err)
  CALL CMISSMeshElements_NodesSet(MeshElements2,1,[1,2,4,5],Err)
  CALL CMISSMeshElements_BasisSet(MeshElements2,2,Basis2,Err)
  CALL CMISSMeshElements_NodesSet(MeshElements2,2,[2,3,5,6],Err)
  CALL CMISSMeshElements_CreateFinish(MeshElements2,Err)

  CALL CMISSMesh_CreateFinish(Mesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Set the geometric field values
  !The x (1) component uses the second mesh component,
  !and the y (2) component uses the first mesh component.
  !So we don't set x derivatives for nodes 1 and 4
  !and y derivatives for nodes 3 and 6
  !X values
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,1,1,1,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,1,2,1,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,2,2,1,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,3,2,1,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,4,2,1,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,1,3,1,2.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,2,3,1,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,3,3,1,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,4,3,1,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,1,4,1,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,1,5,1,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,2,5,1,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,3,5,1,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,4,5,1,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,1,6,1,2.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,2,6,1,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,3,6,1,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,4,6,1,0.0_CMISSDP,Err)
  
  !Y values
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,1,1,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,2,1,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,3,1,2,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,4,1,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,1,2,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,2,2,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,3,2,2,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,4,2,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,1,3,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,1,4,2,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,2,4,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,3,4,2,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,4,4,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,1,5,2,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,2,5,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,3,5,2,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,4,5,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,1,6,2,1.0_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateStart(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
 
  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"MoreComplexMesh","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"MoreComplexMesh","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF

  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM MORECOMPLEXMESHEXAMPLE
