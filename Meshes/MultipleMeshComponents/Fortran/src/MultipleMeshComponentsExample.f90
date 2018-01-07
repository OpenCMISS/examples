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


  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------
  !Change "NumberOfXiCoordinates" to switch between 1/2/3D meshes
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=2
  INTEGER(CMISSIntg), PARAMETER :: FirstBasisInterpolation=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: SecondBasisInterpolation=CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: NumberGlobalXElements=2
  INTEGER(CMISSIntg), PARAMETER :: NumberGlobalYElements=2
  INTEGER(CMISSIntg), PARAMETER :: NumberGlobalZElements=2
  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  !Test program parameters

  REAL(CMISSRP), PARAMETER :: LENGTH=100.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=100.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: HEIGHT=100.0_CMISSRP

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

  TYPE(cmfe_RegionType) :: WorldRegion
  TYPE(cmfe_CoordinateSystemType) :: WorldCoordinateSystem
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
  TYPE(cmfe_RegionType) :: Region
  TYPE(cmfe_BasisType) :: BasisTypes(2)
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_FieldType) :: GeometricField
  TYPE(cmfe_FieldType) :: SecondField
  TYPE(cmfe_FieldsType) :: Fields

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
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  WRITE(*,'(A)') "Program starting."

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains=1

  !Start the creation of a new RC coordinate system for the first region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM << == '
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_TypeSet(CoordinateSystem,CMFE_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL cmfe_CoordinateSystem_OriginSet(CoordinateSystem,[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP],Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the first region
  PRINT *, ' == >> CREATING REGION << == '
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"Region",Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Start the creation of the first Basis type
  PRINT *, ' == >> CREATING BASIS(1) << == '
  CALL cmfe_Basis_Initialise(BasisTypes(1),Err)
  CALL cmfe_Basis_CreateStart(Basis1UserNumber,BasisTypes(1),Err)
  CALL cmfe_Basis_TypeSet(BasisTypes(1),CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(BasisTypes(1),NumberOfXiCoordinates,Err)
  SELECT CASE(NumberOfXiCoordinates)
  CASE(1)
    CALL cmfe_Basis_InterpolationXiSet(BasisTypes(1),[FirstBasisInterpolation],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisTypes(1),[CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
  CASE(2)
    CALL cmfe_Basis_InterpolationXiSet(BasisTypes(1),[FirstBasisInterpolation, &
      & FirstBasisInterpolation],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisTypes(1), &
      & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
  CASE(3)
     CALL cmfe_Basis_InterpolationXiSet(BasisTypes(1),[FirstBasisInterpolation, &
    & FirstBasisInterpolation,FirstBasisInterpolation],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisTypes(1), &
    & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME, &
    & CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
  END SELECT
  CALL cmfe_Basis_CreateFinish(BasisTypes(1),Err)

  !Start the creation of the second Basis type
  PRINT *, ' == >> CREATING BASIS(2) << == '
  CALL cmfe_Basis_Initialise(BasisTypes(2),Err)
  CALL cmfe_Basis_CreateStart(Basis2UserNumber,BasisTypes(2),Err)
  CALL cmfe_Basis_TypeSet(BasisTypes(2),CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(BasisTypes(2),NumberOfXiCoordinates,Err)
  SELECT CASE(NumberOfXiCoordinates)
  CASE(1)
    CALL cmfe_Basis_InterpolationXiSet(BasisTypes(2),[SecondBasisInterpolation],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisTypes(2),[CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
  CASE(2)
    CALL cmfe_Basis_InterpolationXiSet(BasisTypes(2),[SecondBasisInterpolation, &
      & SecondBasisInterpolation],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisTypes(2), &
      & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
  CASE(3)
    CALL cmfe_Basis_InterpolationXiSet(BasisTypes(2),[SecondBasisInterpolation, &
    & SecondBasisInterpolation,SecondBasisInterpolation],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisTypes(2), &
    & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME, &
    & CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
  END SELECT
  CALL cmfe_Basis_CreateFinish(BasisTypes(2),Err)

  !Start the creation of a generated mesh in the first region
  PRINT *, ' == >> CREATING GENERATED MESH << == '
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,BasisTypes,Err)
  CALL cmfe_GeneratedMesh_OriginSet(GeneratedMesh,[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP],Err)
  SELECT CASE(NumberOfXiCoordinates)
  CASE(1)
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH,0.0_CMISSRP,0.0_CMISSRP],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements],Err)
  CASE(2)
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH,WIDTH,0.0_CMISSRP],Err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements, &
    & NumberGlobalYElements],Err)
  CASE(3)
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH,WIDTH,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements, &
      & NumberGlobalYElements,NumberGlobalZElements],Err)
  END SELECT
  !Finish the creation of a generated mesh in the first region
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition for mesh
  PRINT *, ' == >> CREATING MESH DECOMPOSITION << == '
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the first region
  PRINT *, ' == >> CREATING MESH GEOMETRIC FIELD << == '
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  !Set the mesh component to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,2,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,2,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,2,Err)
  !Finish creating the first field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Create a second extra field to use the other mesh component
  PRINT *, ' == >> CREATING SECOND FIELD << == '
  CALL cmfe_Field_Initialise(SecondField,Err)
  CALL cmfe_Field_CreateStart(SecondFieldUserNumber,Region,SecondField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(SecondField,Decomposition,Err)
  CALL cmfe_Field_TypeSet(SecondField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_GeometricFieldSet(SecondField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(SecondField,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(SecondField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_VariableLabelSet(SecondField,CMFE_FIELD_U_VARIABLE_TYPE,"Extra Field",Err)
  !Set the mesh component to be used by the field component.
  CALL cmfe_Field_ComponentMeshComponentSet(SecondField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(SecondField,Err)
  CALL cmfe_Field_ComponentValuesInitialise(SecondField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP,Err)

  !Update the geometric field parameters for the first field
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Output solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"MultipleMeshComponents","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"MultipleMeshComponents","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR

END PROGRAM MULTIPLEMESHCOMPONENTSEXAMPLE
