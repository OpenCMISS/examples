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
  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------
  
  !Test program parameters
  
  REAL(CMISSDP), PARAMETER :: LENGTH=100.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=100.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: HEIGHT=100.0_CMISSDP 

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystem1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Region1UserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: Basis1UserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: Basis2UserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMesh1UserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: Mesh1UserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: Decomposition1UserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: GeometricField1UserNumber=19
    
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3
  
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=3
  
  !Program variables
  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,NumberOfDomains

  INTEGER(CMISSIntg) :: MPI_IERROR

  !CMISS variables

  TYPE(CMISSRegionType) :: WorldRegion
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem1
  TYPE(CMISSRegionType) :: Region1
  TYPE(CMISSBasisType) :: Basis1(2)
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh1
  TYPE(CMISSMeshType) :: Mesh1
  TYPE(CMISSDecompositionType) :: Decomposition1
  TYPE(CMISSFieldType) :: GeometricField1
  TYPE(CMISSFieldsType) :: Fields1
  
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
  
  NumberGlobalXElements=2
  NumberGlobalYElements=2
  NumberGlobalZElements=2
  NumberOfDomains=1

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Start the creation of a new RC coordinate system for the first region
  PRINT *, ' == >> CREATING COORDINATE SYSTEM(1) << == '
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem1,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystem1UserNumber,CoordinateSystem1,Err)
  CALL CMISSCoordinateSystemTypeSet(CoordinateSystem1,CMISSCoordinateRectangularCartesianType,Err)
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem1,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystemOriginSet(CoordinateSystem1,[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP],Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem1,Err)
  
  !Start the creation of the first region
  PRINT *, ' == >> CREATING REGION(1) << == '
  CALL CMISSRegionTypeInitialise(Region1,Err)
  CALL CMISSRegionCreateStart(Region1UserNumber,WorldRegion,Region1,Err)
  CALL CMISSRegionLabelSet(Region1,"Region1",Err)
  CALL CMISSRegionCoordinateSystemSet(Region1,CoordinateSystem1,Err)
  CALL CMISSRegionCreateFinish(Region1,Err)
  
  !Start the creation of a tri-linear-Lagrange basis
  PRINT *, ' == >> CREATING BASIS(1,1) << == '
  CALL CMISSBasisTypeInitialise(Basis1(1),Err)
  CALL CMISSBasisCreateStart(Basis1UserNumber,Basis1(1),Err)
  CALL CMISSBasisTypeSet(Basis1(1),CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(Basis1(1),NumberOfXiCoordinates,Err)
  SELECT CASE(NumberOfXiCoordinates)
  CASE(1)
    CALL CMISSBasisInterpolationXiSet(Basis1(1),[CMISSBasisQuadraticLagrangeInterpolation],Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis1(1),[CMISSBasisMidQuadratureScheme],Err)
  CASE(2)
    CALL CMISSBasisInterpolationXiSet(Basis1(1),[CMISSBasisQuadraticLagrangeInterpolation, &
      & CMISSBasisQuadraticLagrangeInterpolation],Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis1(1), &
      & [CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme],Err)
  CASE(3)
     CALL CMISSBasisInterpolationXiSet(Basis1(1),[CMISSBasisQuadraticLagrangeInterpolation, &
    & CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation],Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis1(1), &
    & [CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme, &
    & CMISSBasisMidQuadratureScheme],Err)
  END SELECT
  CALL CMISSBasisCreateFinish(Basis1(1),Err)

  !Start the creation of a tri-linear-Lagrange basis
  PRINT *, ' == >> CREATING BASIS(1,2) << == '
  CALL CMISSBasisTypeInitialise(Basis1(2),Err)
  CALL CMISSBasisCreateStart(Basis2UserNumber,Basis1(2),Err)
  CALL CMISSBasisTypeSet(Basis1(2),CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(Basis1(2),NumberOfXiCoordinates,Err)
  SELECT CASE(NumberOfXiCoordinates)
  CASE(1)
    CALL CMISSBasisInterpolationXiSet(Basis1(2),[CMISSBasisCubicLagrangeInterpolation],Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis1(2),[CMISSBasisMidQuadratureScheme],Err)
  CASE(2)
    CALL CMISSBasisInterpolationXiSet(Basis1(2),[CMISSBasisCubicLagrangeInterpolation, &
      & CMISSBasisCubicLagrangeInterpolation],Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis1(2), &
      & [CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme],Err)
  CASE(3)
    CALL CMISSBasisInterpolationXiSet(Basis1(2),[CMISSBasisCubicLagrangeInterpolation, &
    & CMISSBasisCubicLagrangeInterpolation,CMISSBasisCubicLagrangeInterpolation],Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis1(2), &
    & [CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme, &
    & CMISSBasisMidQuadratureScheme],Err)
  END SELECT
  CALL CMISSBasisCreateFinish(Basis1(2),Err)
  
  !Start the creation of a generated mesh in the first region
  PRINT *, ' == >> CREATING GENERATED MESH(1) << == '
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh1,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMesh1UserNumber,Region1,GeneratedMesh1,Err)
  !Set up a regular x*y mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh1,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh1,Basis1,Err)   
  CALL CMISSGeneratedMeshOriginSet(GeneratedMesh1,[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP],Err)
  SELECT CASE(NumberOfXiCoordinates)
  CASE(1)
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh1,[LENGTH,0.0_CMISSDP,0.0_CMISSDP],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh1,[NumberGlobalXElements],Err)   
  CASE(2)
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh1,[LENGTH,WIDTH,0.0_CMISSDP],Err)
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh1,[NumberGlobalXElements, &
    & NumberGlobalYElements],Err)   
  CASE(3)
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh1,[LENGTH,WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh1,[NumberGlobalXElements, &
      & NumberGlobalYElements,NumberGlobalZElements],Err)   
  END SELECT
  !Finish the creation of a generated mesh in the first region
  CALL CMISSMeshTypeInitialise(Mesh1,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh1,Mesh1UserNumber,Mesh1,Err)
 
  !Create a decomposition for mesh1
  PRINT *, ' == >> CREATING MESH(1) DECOMPOSITION << == '
  CALL CMISSDecompositionTypeInitialise(Decomposition1,Err)
  CALL CMISSDecompositionCreateStart(Decomposition1UserNumber,Mesh1,Decomposition1,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition1,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition1,NumberOfDomains,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition1,Err)
 
  !Start to create a default (geometric) field on the first region
  PRINT *, ' == >> CREATING MESH(1) GEOMETRIC FIELD << == '
  CALL CMISSFieldTypeInitialise(GeometricField1,Err)
  CALL CMISSFieldCreateStart(GeometricField1UserNumber,Region1,GeometricField1,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField1,Decomposition1,Err)
  CALL CMISSFieldVariableLabelSet(GeometricField1,CMISSFieldUVariableType,"Geometry",Err)
  CALL CMISSFieldTypeSet(GeometricField1,CMISSFieldGeometricType,Err)  
  CALL CMISSFieldNumberOfVariablesSet(GeometricField1,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField1,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)  
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField1,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField1,CMISSFieldUVariableType,2,2,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField1,CMISSFieldUVariableType,3,2,Err)
  !Finish creating the first field
  CALL CMISSFieldCreateFinish(GeometricField1,Err)
  
  !Update the geometric field parameters for the first field
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField1,GeneratedMesh1,Err)
  
  !Output solution  
  CALL CMISSFieldsTypeInitialise(Fields1,Err)
  CALL CMISSFieldsTypeCreate(Region1,Fields1,Err)
  CALL CMISSFieldIONodesExport(Fields1,"MultipleMeshComponents","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields1,"MultipleMeshComponents","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields1,Err)
  

  WRITE(*,'(A)') "Program successfully completed."

  STOP
  
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
  
END PROGRAM MULTIPLEMESHCOMPONENTSEXAMPLE
