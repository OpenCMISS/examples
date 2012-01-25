!> \file
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
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

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
  !Set the Label
  CALL CMISSRegion_LabelSet(Region,"1DCubicHermiteRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Start the creation of a bilinear-Lagrange basis
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  !Set the basis to be a cubic Hermite basis
  CALL CMISSBasis_NumberOfXiSet(Basis,1,Err)
  CALL CMISSBasis_InterpolationXiSet(Basis,[CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION],Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(Basis,Err)
  
  !Create a mesh
  !The mesh will consist of three elements. 
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(Region,4,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)
  
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,1,Mesh,Err)
  CALL CMISSMesh_NumberOfElementsSet(Mesh,3,Err)
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,1,Err)
  
  CALL CMISSMeshElements_Initialise(MeshElements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,1,Basis,MeshElements,Err)
  CALL CMISSMeshElements_NodesSet(MeshElements,1,[1,2],Err)
  CALL CMISSMeshElements_NodesSet(MeshElements,2,[2,3],Err)
  CALL CMISSMeshElements_NodesSet(MeshElements,3,[3,4],Err)
  CALL CMISSMeshElements_CreateFinish(MeshElements,Err)

  CALL CMISSMesh_CreateFinish(Mesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)
  
  !Start to create a geometric field using unit scaling on the region
  CALL CMISSField_Initialise(GeometricField1,Err)
  CALL CMISSField_CreateStart(GeometricField1UserNumber,Region,GeometricField1,Err)
  !Set the label
  CALL CMISSField_VariableLabelSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,"UnitScaling",Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField1,Decomposition,Err)
  !Set the field Scaling
  CALL CMISSField_ScalingTypeSet(GeometricField1,CMISS_FIELD_UNIT_SCALING,Err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField1,Err)

  !Start to create a geometric field using arithmetic mean scaling on the region
  CALL CMISSField_Initialise(GeometricField2,Err)
  CALL CMISSField_CreateStart(GeometricField2UserNumber,Region,GeometricField2,Err)
  !Set the label
  CALL CMISSField_VariableLabelSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,"ArithmeticMeanScaling",Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField2,Decomposition,Err)
  !Set the field Scaling
  CALL CMISSField_ScalingTypeSet(GeometricField2,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField2,Err)

  !Start to create a geometric field using geometric mean scaling on the region
  CALL CMISSField_Initialise(GeometricField3,Err)
  CALL CMISSField_CreateStart(GeometricField3UserNumber,Region,GeometricField3,Err)
  !Set the label
  CALL CMISSField_VariableLabelSet(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,"GeometricMeanScaling",Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField3,Decomposition,Err)
  !Set the field Scaling
  CALL CMISSField_ScalingTypeSet(GeometricField3,CMISS_FIELD_GEOMETRIC_MEAN_SCALING,Err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField3,Err)

  !Start to create a geometric field using harmonic mean scaling on the region
  CALL CMISSField_Initialise(GeometricField4,Err)
  CALL CMISSField_CreateStart(GeometricField4UserNumber,Region,GeometricField4,Err)
  !Set the label
  CALL CMISSField_VariableLabelSet(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,"HarmonicMeanScaling",Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField4,Decomposition,Err)
  !Set the field Scaling
  CALL CMISSField_ScalingTypeSet(GeometricField4,CMISS_FIELD_HARMONIC_MEAN_SCALING,Err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField4,Err)

  !Set the geometric field values for geometric field 1
  !Node 1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,1,1,Half,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,1,2,Half,Err)
  !Node 2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,1, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,2,1,Half,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,2,2,-Half,Err)
  !Node 3
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,1, &
    & 3.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,3,1,OneThird, &
    & Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,3,2,TwoThird, &
    & Err)
  !Node 4
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,1, &
    & 7.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,4,1,OneThird, &
    & Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,4,2, &
    & -TwoThird, &
    & Err)
 
  !Set the geometric field values for geometric field 2
  !Node 1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,1,1,Half,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,1,2,Half,Err)
  !Node 2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,1, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,2,1,Half,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,2,2,-Half,Err)
  !Node 3
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,1, &
    & 3.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,3,1,OneThird, &
    & Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,3,2,TwoThird, &
    & Err)
  !Node 4
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,1, &
    & 7.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,4,1,OneThird, &
    & Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,4,2, &
    & -TwoThird, &
    & Err)
 
  !Set the geometric field values for geometric field 3
  !Node 1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,1,1,Half,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,1,2,Half,Err)
  !Node 2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,1, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,2,1,Half,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,2,2,-Half,Err)
  !Node 3
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,1, &
    & 3.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,3,1,OneThird, &
    & Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,3,2,TwoThird, &
    & Err)
  !Node 4
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,1, &
    & 7.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,4,1,OneThird, &
    & Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,4,2, &
    & -TwoThird, &
    & Err)

  !Set the geometric field values for geometric field 4
  !Node 1
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,1,1,Half,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,1,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,1,2,Half,Err)
  !Node 2
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,1, &
    & 1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,2,1,Half,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,2,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,2,2,-Half,Err)
  !Node 3
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,1, &
    & 3.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,3,1,OneThird, &
    & Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,3,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,3,2,TwoThird, &
    & Err)
  !Node 4
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,1, &
    & 7.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,4,1,OneThird, &
    & Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,4,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2,4,2, &
    & -TwoThird, &
    & Err)

  !Update fields
  CALL CMISSField_ParameterSetUpdateStart(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateStart(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateStart(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField3,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateStart(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField4,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
 
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"1DCubicHermite","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"1DCubicHermite","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM ONEDCUBICHERMITEEXAMPLE
