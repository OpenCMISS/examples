!> \file
!> \author Xiani Nancy Yan
!> \brief This is an example program which sets up a generated mesh (1/2/3D) with multiple mesh components
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
PROGRAM MULTIPLEMESHCOMPONENTSEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------
  !Change "NumberOfXiCoordinates" to switch between 1/2/3D meshes
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=2
  INTEGER(CMISSIntg), PARAMETER :: FirstBasisInterpolation=CMISSBasisQuadraticLagrangeInterpolation
  INTEGER(CMISSIntg), PARAMETER :: SecondBasisInterpolation=CMISSBasisCubicLagrangeInterpolation
  INTEGER(CMISSIntg), PARAMETER :: NumberGlobalXElements=2
  INTEGER(CMISSIntg), PARAMETER :: NumberGlobalYElements=2
  INTEGER(CMISSIntg), PARAMETER :: NumberGlobalZElements=2
  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: LENGTH=100.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=100.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: HEIGHT=100.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: Basis1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Basis2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=19
  INTEGER(CMISSIntg), PARAMETER :: SecondFieldUserNumber=20

  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=3

  !Program variables
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,NumberOfDomains

  !CMISS variables

  TYPE(CMISSRegionType) :: WorldRegion
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
  TYPE(CMISSRegionType) :: Region
  TYPE(CMISSBasisType) :: BasisTypes(2)
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: SecondField
  TYPE(CMISSFieldsType) :: Fields

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

  !Intialise cmiss
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  WRITE(*,'(A)') "Program starting."

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains=1

  !Start the creation of a new RC coordinate system for the first region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM << == '
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemTypeSet(CoordinateSystem,CMISSCoordinateRectangularCartesianType,Err)
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystemOriginSet(CoordinateSystem,[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP],Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Start the creation of the first region
  PRINT *, ' == >> CREATING REGION << == '
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionLabelSet(Region,"Region",Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !Start the creation of the first Basis type
  PRINT *, ' == >> CREATING BASIS(1) << == '
  CALL CMISSBasisTypeInitialise(BasisTypes(1),Err)
  CALL CMISSBasisCreateStart(Basis1UserNumber,BasisTypes(1),Err)
  CALL CMISSBasisTypeSet(BasisTypes(1),CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(BasisTypes(1),NumberOfXiCoordinates,Err)
  SELECT CASE(NumberOfXiCoordinates)
  CASE(1)
    CALL CMISSBasisInterpolationXiSet(BasisTypes(1),[FirstBasisInterpolation],Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisTypes(1),[CMISSBasisMidQuadratureScheme],Err)
  CASE(2)
    CALL CMISSBasisInterpolationXiSet(BasisTypes(1),[FirstBasisInterpolation, &
      & FirstBasisInterpolation],Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisTypes(1), &
      & [CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme],Err)
  CASE(3)
     CALL CMISSBasisInterpolationXiSet(BasisTypes(1),[FirstBasisInterpolation, &
    & FirstBasisInterpolation,FirstBasisInterpolation],Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisTypes(1), &
    & [CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme, &
    & CMISSBasisMidQuadratureScheme],Err)
  END SELECT
  CALL CMISSBasisCreateFinish(BasisTypes(1),Err)

  !Start the creation of the second Basis type
  PRINT *, ' == >> CREATING BASIS(2) << == '
  CALL CMISSBasisTypeInitialise(BasisTypes(2),Err)
  CALL CMISSBasisCreateStart(Basis2UserNumber,BasisTypes(2),Err)
  CALL CMISSBasisTypeSet(BasisTypes(2),CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(BasisTypes(2),NumberOfXiCoordinates,Err)
  SELECT CASE(NumberOfXiCoordinates)
  CASE(1)
    CALL CMISSBasisInterpolationXiSet(BasisTypes(2),[SecondBasisInterpolation],Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisTypes(2),[CMISSBasisMidQuadratureScheme],Err)
  CASE(2)
    CALL CMISSBasisInterpolationXiSet(BasisTypes(2),[SecondBasisInterpolation, &
      & SecondBasisInterpolation],Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisTypes(2), &
      & [CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme],Err)
  CASE(3)
    CALL CMISSBasisInterpolationXiSet(BasisTypes(2),[SecondBasisInterpolation, &
    & SecondBasisInterpolation,SecondBasisInterpolation],Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisTypes(2), &
    & [CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme, &
    & CMISSBasisMidQuadratureScheme],Err)
  END SELECT
  CALL CMISSBasisCreateFinish(BasisTypes(2),Err)

  !Start the creation of a generated mesh in the first region
  PRINT *, ' == >> CREATING GENERATED MESH << == '
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,BasisTypes,Err)
  CALL CMISSGeneratedMeshOriginSet(GeneratedMesh,[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP],Err)
  SELECT CASE(NumberOfXiCoordinates)
  CASE(1)
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[LENGTH,0.0_CMISSDP,0.0_CMISSDP],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements],Err)
  CASE(2)
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[LENGTH,WIDTH,0.0_CMISSDP],Err)
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements, &
    & NumberGlobalYElements],Err)
  CASE(3)
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[LENGTH,WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements, &
      & NumberGlobalYElements,NumberGlobalZElements],Err)
  END SELECT
  !Finish the creation of a generated mesh in the first region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition for mesh
  PRINT *, ' == >> CREATING MESH DECOMPOSITION << == '
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the first region
  PRINT *, ' == >> CREATING MESH GEOMETRIC FIELD << == '
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)
  CALL CMISSFieldVariableLabelSet(GeometricField,CMISSFieldUVariableType,"Geometry",Err)
  !Set the mesh component to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,2,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,2,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,2,Err)
  !Finish creating the first field
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Create a second extra field to use the other mesh component
  PRINT *, ' == >> CREATING SECOND FIELD << == '
  CALL CMISSFieldTypeInitialise(SecondField,Err)
  CALL CMISSFieldCreateStart(SecondFieldUserNumber,Region,SecondField,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(SecondField,Decomposition,Err)
  CALL CMISSFieldTypeSet(SecondField,CMISSFieldGeneralType,Err)
  CALL CMISSFieldGeometricFieldSet(SecondField,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(SecondField,1,Err)
  CALL CMISSFieldNumberOfComponentsSet(SecondField,CMISSFieldUVariableType,1,Err)
  CALL CMISSFieldVariableLabelSet(SecondField,CMISSFieldUVariableType,"Extra Field",Err)
  !Set the mesh component to be used by the field component.
  CALL CMISSFieldComponentMeshComponentSet(SecondField,CMISSFieldUVariableType,1,1,Err)
  !Finish creating the field
  CALL CMISSFieldCreateFinish(SecondField,Err)
  CALL CMISSFieldComponentValuesInitialise(SecondField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,0.0_CMISSDP,Err)

  !Update the geometric field parameters for the first field
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

  !Output solution
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"MultipleMeshComponents","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"MultipleMeshComponents","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR

END PROGRAM MULTIPLEMESHCOMPONENTSEXAMPLE
