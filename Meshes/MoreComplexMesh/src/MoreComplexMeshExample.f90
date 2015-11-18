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

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


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

  TYPE(cmfe_BasisType) :: Basis1,Basis2
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_FieldType) :: GeometricField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_MeshElementsType) :: MeshElements1,MeshElements2
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_RegionType) :: Region,WorldRegion

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
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  
  NUMBER_OF_DOMAINS=1

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 2D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Start the creation of a bilinear-Lagrange basis
  CALL cmfe_Basis_Initialise(Basis1,Err)
  CALL cmfe_Basis_CreateStart(Basis1UserNumber,Basis1,Err)
  !Set the basis to be a bilinear Lagrange basis
  CALL cmfe_Basis_NumberOfXiSet(Basis1,2,Err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis1,Err)
  
  !Start the creation of a bicubic-Hermite basis
  CALL cmfe_Basis_Initialise(Basis2,Err)
  CALL cmfe_Basis_CreateStart(Basis2UserNumber,Basis2,Err)
  !Set the basis to be a bilinear Lagrange basis
  CALL cmfe_Basis_NumberOfXiSet(Basis2,2,Err)
  !Set the basis to be a bicubic Hermite basis
  CALL cmfe_Basis_InterpolationXiSet(Basis2,[CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION],Err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis2,Err)

  !Create a mesh
  !The mesh will consist of two elements. There will be two mesh components. For the first component the first element is
  !bilinear Lagrange and the second element bicubic Hermite. For the second component the first element is bicubic Hermite and
  !the second element is bilinear Lagrange.
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,6,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,2,Mesh,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,2,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,2,Err)
  
  CALL cmfe_MeshElements_Initialise(MeshElements1,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,1,Basis2,MeshElements1,Err)
  CALL cmfe_MeshElements_NodesSet(MeshElements1,1,[1,2,4,5],Err)
  CALL cmfe_MeshElements_BasisSet(MeshElements1,2,Basis1,Err)
  CALL cmfe_MeshElements_NodesSet(MeshElements1,2,[2,3,5,6],Err)
  CALL cmfe_MeshElements_CreateFinish(MeshElements1,Err)

  CALL cmfe_MeshElements_Initialise(MeshElements2,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,2,Basis1,MeshElements2,Err)
  CALL cmfe_MeshElements_NodesSet(MeshElements2,1,[1,2,4,5],Err)
  CALL cmfe_MeshElements_BasisSet(MeshElements2,2,Basis2,Err)
  CALL cmfe_MeshElements_NodesSet(MeshElements2,2,[2,3,5,6],Err)
  CALL cmfe_MeshElements_CreateFinish(MeshElements2,Err)

  CALL cmfe_Mesh_CreateFinish(Mesh,Err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)
  
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Set the geometric field values
  !The x (1) component uses the second mesh component,
  !and the y (2) component uses the first mesh component.
  !So we don't set x derivatives for nodes 1 and 4
  !and y derivatives for nodes 3 and 6
  !X values
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,1,1,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,2,1,1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,2,2,1,1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,3,2,1,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,4,2,1,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,3,1,2.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,2,3,1,1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,3,3,1,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,4,3,1,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,4,1,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,5,1,1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,2,5,1,1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,3,5,1,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,4,5,1,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,6,1,2.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,2,6,1,1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,3,6,1,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,4,6,1,0.0_CMISSRP,Err)
  
  !Y values
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,1,2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,2,1,2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,3,1,2,1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,4,1,2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,2,2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,2,2,2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,3,2,2,1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,4,2,2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,3,2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,4,2,1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,2,4,2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,3,4,2,1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,4,4,2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,5,2,1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,2,5,2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,3,5,2,1.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,4,5,2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,6,2,1.0_CMISSRP,Err)

  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
 
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"MoreComplexMesh","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"MoreComplexMesh","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF

  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM MORECOMPLEXMESHEXAMPLE
