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
  INTEGER(CMISSIntg), PARAMETER :: FirstBasisInterpolation=CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: SecondBasisInterpolation=CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION
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

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  WRITE(*,'(A)') "Program starting."

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains=1

  !Start the creation of a new RC coordinate system for the first region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM << == '
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_TypeSet(CoordinateSystem,CMISS_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystem_OriginSet(CoordinateSystem,[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP],Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the first region
  PRINT *, ' == >> CREATING REGION << == '
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_LabelSet(Region,"Region",Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Start the creation of the first Basis type
  PRINT *, ' == >> CREATING BASIS(1) << == '
  CALL CMISSBasis_Initialise(BasisTypes(1),Err)
  CALL CMISSBasis_CreateStart(Basis1UserNumber,BasisTypes(1),Err)
  CALL CMISSBasis_TypeSet(BasisTypes(1),CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(BasisTypes(1),NumberOfXiCoordinates,Err)
  SELECT CASE(NumberOfXiCoordinates)
  CASE(1)
    CALL CMISSBasis_InterpolationXiSet(BasisTypes(1),[FirstBasisInterpolation],Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisTypes(1),[CMISS_BASIS_MID_QUADRATURE_SCHEME],Err)
  CASE(2)
    CALL CMISSBasis_InterpolationXiSet(BasisTypes(1),[FirstBasisInterpolation, &
      & FirstBasisInterpolation],Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisTypes(1), &
      & [CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME],Err)
  CASE(3)
     CALL CMISSBasis_InterpolationXiSet(BasisTypes(1),[FirstBasisInterpolation, &
    & FirstBasisInterpolation,FirstBasisInterpolation],Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisTypes(1), &
    & [CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME, &
    & CMISS_BASIS_MID_QUADRATURE_SCHEME],Err)
  END SELECT
  CALL CMISSBasis_CreateFinish(BasisTypes(1),Err)

  !Start the creation of the second Basis type
  PRINT *, ' == >> CREATING BASIS(2) << == '
  CALL CMISSBasis_Initialise(BasisTypes(2),Err)
  CALL CMISSBasis_CreateStart(Basis2UserNumber,BasisTypes(2),Err)
  CALL CMISSBasis_TypeSet(BasisTypes(2),CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(BasisTypes(2),NumberOfXiCoordinates,Err)
  SELECT CASE(NumberOfXiCoordinates)
  CASE(1)
    CALL CMISSBasis_InterpolationXiSet(BasisTypes(2),[SecondBasisInterpolation],Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisTypes(2),[CMISS_BASIS_MID_QUADRATURE_SCHEME],Err)
  CASE(2)
    CALL CMISSBasis_InterpolationXiSet(BasisTypes(2),[SecondBasisInterpolation, &
      & SecondBasisInterpolation],Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisTypes(2), &
      & [CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME],Err)
  CASE(3)
    CALL CMISSBasis_InterpolationXiSet(BasisTypes(2),[SecondBasisInterpolation, &
    & SecondBasisInterpolation,SecondBasisInterpolation],Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisTypes(2), &
    & [CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME, &
    & CMISS_BASIS_MID_QUADRATURE_SCHEME],Err)
  END SELECT
  CALL CMISSBasis_CreateFinish(BasisTypes(2),Err)

  !Start the creation of a generated mesh in the first region
  PRINT *, ' == >> CREATING GENERATED MESH << == '
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,BasisTypes,Err)
  CALL CMISSGeneratedMesh_OriginSet(GeneratedMesh,[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP],Err)
  SELECT CASE(NumberOfXiCoordinates)
  CASE(1)
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH,0.0_CMISSDP,0.0_CMISSDP],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements],Err)
  CASE(2)
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH,WIDTH,0.0_CMISSDP],Err)
  CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements, &
    & NumberGlobalYElements],Err)
  CASE(3)
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH,WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements, &
      & NumberGlobalYElements,NumberGlobalZElements],Err)
  END SELECT
  !Finish the creation of a generated mesh in the first region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition for mesh
  PRINT *, ' == >> CREATING MESH DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the first region
  PRINT *, ' == >> CREATING MESH GEOMETRIC FIELD << == '
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)
  CALL CMISSField_VariableLabelSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  !Set the mesh component to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,2,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,2,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,2,Err)
  !Finish creating the first field
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Create a second extra field to use the other mesh component
  PRINT *, ' == >> CREATING SECOND FIELD << == '
  CALL CMISSField_Initialise(SecondField,Err)
  CALL CMISSField_CreateStart(SecondFieldUserNumber,Region,SecondField,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(SecondField,Decomposition,Err)
  CALL CMISSField_TypeSet(SecondField,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_GeometricFieldSet(SecondField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(SecondField,1,Err)
  CALL CMISSField_NumberOfComponentsSet(SecondField,CMISS_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMISSField_VariableLabelSet(SecondField,CMISS_FIELD_U_VARIABLE_TYPE,"Extra Field",Err)
  !Set the mesh component to be used by the field component.
  CALL CMISSField_ComponentMeshComponentSet(SecondField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(SecondField,Err)
  CALL CMISSField_ComponentValuesInitialise(SecondField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)

  !Update the geometric field parameters for the first field
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Output solution
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"MultipleMeshComponents","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"MultipleMeshComponents","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR

END PROGRAM MULTIPLEMESHCOMPONENTSEXAMPLE
