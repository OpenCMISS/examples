!> \file
!> $Id$
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
  
  NUMBER_OF_DOMAINS=1

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 2D
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,2,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)

  !Start the creation of a bilinear-Lagrange basis
  CALL CMISSBasisTypeInitialise(Basis1,Err)
  CALL CMISSBasisCreateStart(Basis1UserNumber,Basis1,Err)
  !Set the basis to be a bilinear Lagrange basis
  CALL CMISSBasisNumberOfXiSet(Basis1,2,Err)
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(Basis1,Err)
  
  !Start the creation of a bicubic-Hermite basis
  CALL CMISSBasisTypeInitialise(Basis2,Err)
  CALL CMISSBasisCreateStart(Basis2UserNumber,Basis2,Err)
  !Set the basis to be a bilinear Lagrange basis
  CALL CMISSBasisNumberOfXiSet(Basis2,2,Err)
  !Set the basis to be a bicubic Hermite basis
  CALL CMISSBasisInterpolationXiSet(Basis2,(/CMISSBasisCubicHermiteInterpolation,CMISSBasisCubicHermiteInterpolation/),Err)
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(Basis2,Err)

  !Create a mesh
  !The mesh will consist of two elements. There will be two mesh components. For the first component the first element is
  !bilinear Lagrange and the second element bicubic Hermite. For the second component the first element is bicubic Hermite and
  !the second element is bilinear Lagrange.
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSNodesCreateStart(Region,6,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)
  
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,2,Mesh,Err)
  CALL CMISSMeshNumberOfElementsSet(Mesh,2,Err)
  CALL CMISSMeshNumberOfComponentsSet(Mesh,2,Err)
  
  CALL CMISSMeshElementsTypeInitialise(MeshElements1,Err)
  CALL CMISSMeshElementsCreateStart(Mesh,1,Basis1,MeshElements1,Err)
  CALL CMISSMeshElementsNodesSet(MeshElements1,1,(/1,2,4,5/),Err)
  CALL CMISSMeshElementsBasisSet(MeshElements1,2,Basis2,Err)
  CALL CMISSMeshElementsNodesSet(MeshElements1,2,(/2,3,5,6/),Err)
  CALL CMISSMeshElementsCreateFinish(MeshElements1,Err)

  CALL CMISSMeshElementsTypeInitialise(MeshElements2,Err)
  CALL CMISSMeshElementsCreateStart(Mesh,2,Basis2,MeshElements2,Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,1,(/1,2,4,5/),Err)
  CALL CMISSMeshElementsBasisSet(MeshElements2,2,Basis1,Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,2,(/2,3,5,6/),Err)
  CALL CMISSMeshElementsCreateFinish(MeshElements2,Err)

  CALL CMISSMeshCreateFinish(Mesh,Err)

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,1,Err)
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Set the geometric field values
  !X values
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,2,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,2,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,2,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,3,1,2.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,3,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,3,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,3,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,4,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,5,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,5,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,5,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,5,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,6,1,2.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,6,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,6,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,6,1,0.0_CMISSDP,Err)
  
  !Y values
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,1,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,1,2,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,1,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,2,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,2,2,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,2,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,3,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,4,2,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,4,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,4,2,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,4,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,5,2,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,5,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,5,2,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,5,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,6,2,1.0_CMISSDP,Err)

  CALL CMISSFieldParameterSetUpdateStart(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
 
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"MoreComplexMesh","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"MoreComplexMesh","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
  ENDIF

  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM MORECOMPLEXMESHEXAMPLE
