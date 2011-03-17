!> \file
!> $Id: 1DCubicHermiteExample.f90 1528 2010-09-21 01:32:29Z chrispbradley $
!> \author Chris Bradley
!> \brief This is an example program which sets up a 1D cubic Hermite mesh using OpenCMISS calls.
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

!> Main program
PROGRAM ONEDCUBICHERMITEEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: GeometricField1UserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricField2UserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: GeometricField3UserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: GeometricField4UserNumber=9

  REAL(CMISSDP), PARAMETER :: PI=3.141592653589793238462643383279502884197_CMISSDP
  REAL(CMISSDP), PARAMETER :: Half = 1.0_CMISSDP/SQRT(2.0_CMISSDP)
  REAL(CMISSDP), PARAMETER :: OneThird = COS(PI/3.0_CMISSDP)
  REAL(CMISSDP), PARAMETER :: TwoThird = SIN(PI/3.0_CMISSDP)
 
  !Program types
  
  !Program variables
  
  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSFieldType) :: GeometricField1,GeometricField2,GeometricField3,GeometricField4
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSMeshElementsType) :: MeshElements
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSRegionType) :: Region,WorldRegion

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
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

  !Set error handling mode to trap all errors
  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

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
  !Set the Label
  CALL CMISSRegionLabelSet(Region,"1DCubicHermiteRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)

  !Start the creation of a bilinear-Lagrange basis
  CALL CMISSBasisTypeInitialise(Basis,Err)
  CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
  !Set the basis to be a cubic Hermite basis
  CALL CMISSBasisNumberOfXiSet(Basis,1,Err)
  CALL CMISSBasisInterpolationXiSet(Basis,[CMISSBasisCubicHermiteInterpolation],Err)
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(Basis,Err)
  
  !Create a mesh
  !The mesh will consist of three elements. 
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSNodesCreateStart(Region,4,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)
  
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,1,Mesh,Err)
  CALL CMISSMeshNumberOfElementsSet(Mesh,3,Err)
  CALL CMISSMeshNumberOfComponentsSet(Mesh,1,Err)
  
  CALL CMISSMeshElementsTypeInitialise(MeshElements,Err)
  CALL CMISSMeshElementsCreateStart(Mesh,1,Basis,MeshElements,Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,1,[1,2],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,2,[2,3],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,3,[3,4],Err)
  CALL CMISSMeshElementsCreateFinish(MeshElements,Err)

  CALL CMISSMeshCreateFinish(Mesh,Err)

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)
  
  !Start to create a geometric field using unit scaling on the region
  CALL CMISSFieldTypeInitialise(GeometricField1,Err)
  CALL CMISSFieldCreateStart(GeometricField1UserNumber,Region,GeometricField1,Err)
  !Set the label
  CALL CMISSFieldVariableLabelSet(GeometricField1,CMISSFieldUVariableType,"UnitScaling",Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField1,Decomposition,Err)
  !Set the field Scaling
  CALL CMISSFieldScalingTypeSet(GeometricField1,CMISSFieldUnitScaling,Err)
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField1,Err)

  !Start to create a geometric field using arithmetic mean scaling on the region
  CALL CMISSFieldTypeInitialise(GeometricField2,Err)
  CALL CMISSFieldCreateStart(GeometricField2UserNumber,Region,GeometricField2,Err)
  !Set the label
  CALL CMISSFieldVariableLabelSet(GeometricField2,CMISSFieldUVariableType,"ArithmeticMeanScaling",Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField2,Decomposition,Err)
  !Set the field Scaling
  CALL CMISSFieldScalingTypeSet(GeometricField2,CMISSFieldArithmeticMeanScaling,Err)
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField2,Err)

  !Start to create a geometric field using geometric mean scaling on the region
  CALL CMISSFieldTypeInitialise(GeometricField3,Err)
  CALL CMISSFieldCreateStart(GeometricField3UserNumber,Region,GeometricField3,Err)
  !Set the label
  CALL CMISSFieldVariableLabelSet(GeometricField3,CMISSFieldUVariableType,"GeometricMeanScaling",Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField3,Decomposition,Err)
  !Set the field Scaling
  CALL CMISSFieldScalingTypeSet(GeometricField3,CMISSFieldGeometricMeanScaling,Err)
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField3,Err)

  !Start to create a geometric field using harmonic mean scaling on the region
  CALL CMISSFieldTypeInitialise(GeometricField4,Err)
  CALL CMISSFieldCreateStart(GeometricField4UserNumber,Region,GeometricField4,Err)
  !Set the label
  CALL CMISSFieldVariableLabelSet(GeometricField4,CMISSFieldUVariableType,"HarmonicMeanScaling",Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField4,Decomposition,Err)
  !Set the field Scaling
  CALL CMISSFieldScalingTypeSet(GeometricField4,CMISSFieldHarmonicMeanScaling,Err)
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField4,Err)

  !Set the geometric field values for geometric field 1
  !Node 1
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,1,1,Half,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,1,2,Half,Err)
  !Node 2
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,2,1,Half,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,2,2,-Half,Err)
  !Node 3
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,1,3.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,3,1,OneThird,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,3,2,TwoThird,Err)
  !Node 4
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,1,5.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,4,1,OneThird,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,4,2,-TwoThird,Err)
 
  !Set the geometric field values for geometric field 2
  !Node 1
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,1,1,Half,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,1,2,Half,Err)
  !Node 2
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,2,1,Half,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,2,2,-Half,Err)
  !Node 3
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,1,3.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,3,1,OneThird,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,3,2,TwoThird,Err)
  !Node 4
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,1,5.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,4,1,OneThird,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,4,2,-TwoThird,Err)
 
  !Set the geometric field values for geometric field 3
  !Node 1
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,1,1,Half,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,1,2,Half,Err)
  !Node 2
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,2,1,Half,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,2,2,-Half,Err)
  !Node 3
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,1,3.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,3,1,OneThird,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,3,2,TwoThird,Err)
  !Node 4
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,1,5.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,4,1,OneThird,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,4,2,-TwoThird,Err)

  !Set the geometric field values for geometric field 4
  !Node 1
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,1,1,Half,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,1,2,Half,Err)
  !Node 2
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,2,1,Half,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,2,2,-Half,Err)
  !Node 3
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,1,3.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,3,1,OneThird,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,3,2,TwoThird,Err)
  !Node 4
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,1,5.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,4,1,OneThird,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2,4,2,-TwoThird,Err)

  !Update fields
  CALL CMISSFieldParameterSetUpdateStart(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField1,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateStart(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField2,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateStart(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField3,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateStart(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField4,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
 
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"1DCubicHermite","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"1DCubicHermite","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM ONEDCUBICHERMITEEXAMPLE
