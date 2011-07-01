!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a finite elasticity equation using openCMISS calls.
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
!> The Original Code is openCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): Jack Lee
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

!> \example FiniteElasticity/UniAxialExtension/src/UniAxialExtensionExample.f90
!! Example program to solve a finite elasticity equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/UniAxialExtension/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/UniAxialExtension/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM QUADRATICELLIPSOIDCOSTAREADINEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters
  
  REAL(CMISSDP) :: PI=3.14159_CMISSDP
  REAL(CMISSDP) :: LONG_AXIS=2.0_CMISSDP
  REAL(CMISSDP) :: SHORT_AXIS=1.0_CMISSDP
  REAL(CMISSDP) :: WALL_THICKNESS=0.5_CMISSDP
  REAL(CMISSDP) :: CUTOFF_ANGLE=1.5708_CMISSDP
  REAL(CMISSDP) :: FIBRE_SLOPE_ENDO=1.73205_CMISSDP !Slope of fibres in endocardium = 60 degrees
  REAL(CMISSDP) :: FIBRE_SLOPE_EPI=-3.4641_CMISSDP !Slope of fibres in endocardium = -60 degrees 
  REAL(CMISSDP) :: SHEET_SLOPE_BASE_ENDO=1.0_CMISSDP !Slope of sheet at base endocardium 
  !REAL(CMISSDP), PARAMETER :: COSTA_PARAMS (1:7) =  (/ 0.2, 30.0, 12.0, 14.0, 14.0, 10.0, 18.0 /)
  !REAL(CMISSDP), DIMENSION(:),ALLOCATABLE :: COSTA_PARAMS ! =  (/ 0.2, 30.0, 12.0, 14.0, 14.0, 10.0, 18.0 /) ! a bff bfs bfn bss bsn bnn
  REAL(CMISSDP) :: COSTA_PARAMS (1:7)
  REAL(CMISSDP), PARAMETER :: INNER_PRESSURE=2.0_CMISSDP  !Positive is compressive
  REAL(CMISSDP), PARAMETER :: OUTER_PRESSURE=0.0_CMISSDP  !Positive is compressive

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS=4  ! X ==NUMBER_GLOBAL_CIRCUMFERENTIAL_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_Y_ELEMENTS=4  ! Y ==NUMBER_GLOBAL_LONGITUDINAL_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_Z_ELEMENTS=1  ! Z ==NUMBER_GLOBAL_TRANSMURAL_ELEMENTS
  INTEGER(CMISSIntg) :: NumberOfDomains

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticCollapsedBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: LinearCollapsedBasisUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DerivativeUserNumber=1
  
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfElements=1

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponents=7

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariables=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponents=4

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program types


  !Program variables

  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  REAL(CMISSDP) :: FibreFieldAngle(3) 
  REAL(CMISSDP) :: nu,theta,omega,XI3,XI3delta,XI2delta, zero
  INTEGER(CMISSIntg) ::i,j,k,component_idx,node_idx,TOTAL_NUMBER_NODES_XI(3)
  !For grabbing surfaces
  INTEGER(CMISSIntg) :: InnerNormalXi,OuterNormalXi,TopNormalXi
  INTEGER(CMISSIntg), ALLOCATABLE :: InnerSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: OuterSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: TopSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: G(:)
  INTEGER(CMISSIntg) :: NotCornerNode, CornerNode, NumberOfCornerNodes, gn, CorrectNodeNumber,TotalNumberOfNodes
  INTEGER(CMISSIntg) :: NN,NODE,NodeDomain
  REAL(CMISSDP) :: XCoord,YCoord,ZCoord
  LOGICAL :: X_FIXED,Y_FIXED,X_OKAY,Y_OKAY
  INTEGER(CMISSIntg) :: PARAMETERSETNR
  !CMISS variables

  TYPE(CMISSBasisType) :: QuadraticBasis,QuadraticCollapsedBasis,LinearBasis,LinearCollapsedBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,FibreField,MaterialField,DependentField,EquationsSetField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: Err
  CHARACTER(len=20) :: HEIGHT, WIDTH, WALL, CUTOFF,X_ELEM,Y_ELEM,Z_ELEM
  CHARACTER(len=20) :: FIBRE_ENDO, FIBRE_EPI, SHEET_ENDO_SLOPE
  CHARACTER(len=20) :: COSTA_1, COSTA_2, COSTA_3, COSTA_4, COSTA_5, COSTA_6, COSTA_7
  CHARACTER(len=4) :: PARAMETERSET
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
  CALL GETARG(1,PARAMETERSET)
  CALL GETARG(2,X_ELEM)
  CALL GETARG(3,Y_ELEM)
  CALL GETARG(4,Z_ELEM)
  CALL GETARG(5,HEIGHT)
  CALL GETARG(6,WIDTH)
  CALL GETARG(7,WALL)
  !CALL GETARG(8,CUTOFF)
  CALL GETARG(8,FIBRE_ENDO)
  CALL GETARG(9,FIBRE_EPI)
  CALL GETARG(10,SHEET_ENDO_SLOPE)
  CALL GETARG(11,COSTA_1)
  CALL GETARG(12,COSTA_2)
  CALL GETARG(13,COSTA_3)
  CALL GETARG(14,COSTA_4)
  CALL GETARG(15,COSTA_5)
  CALL GETARG(16,COSTA_6)
  CALL GETARG(17,COSTA_7)
  READ(HEIGHT,*)LONG_AXIS
  READ(WIDTH,*)SHORT_AXIS
  READ(WALL,*)WALL_THICKNESS
  !READ(CUTOFF,*)CUTOFF_ANGLE
  READ(FIBRE_ENDO,*)FIBRE_SLOPE_ENDO
  READ(FIBRE_EPI,*)FIBRE_SLOPE_EPI
  READ(SHEET_ENDO_SLOPE,*)SHEET_SLOPE_BASE_ENDO
  READ(COSTA_1,*)COSTA_PARAMS(1)
  READ(COSTA_2,*)COSTA_PARAMS(2)
  READ(COSTA_3,*)COSTA_PARAMS(3)
  READ(COSTA_4,*)COSTA_PARAMS(4)
  READ(COSTA_5,*)COSTA_PARAMS(5)
  READ(COSTA_6,*)COSTA_PARAMS(6)
  READ(COSTA_7,*)COSTA_PARAMS(7)
  !Convert PARAMETERSET to integer and measure the number of digits
  !Intialise cmiss
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  CALL CMISSDiagnosticsSetOn(CMISSFromDiagType,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_FINITE_ELEMENT_CALCULATE"/),Err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
!   CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  READ(X_ELEM,*)NUMBER_GLOBAL_X_ELEMENTS
  READ(Y_ELEM,*)NUMBER_GLOBAL_Y_ELEMENTS
  READ(Z_ELEM,*)NUMBER_GLOBAL_Z_ELEMENTS 
  NumberOfDomains=NumberOfComputationalNodes

  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemTypeSet(CoordinateSystem,CMISSCoordinateRectangularCartesianType,Err)
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystemOriginSet(CoordinateSystem,(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/),Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Create a region and assign the CS to the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange, each with collapsed variant
    !Quadratic Basis
  CALL CMISSBasisTypeInitialise(QuadraticBasis,Err)
  CALL CMISSBasisCreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasisInterpolationXiSet(QuadraticBasis,(/CMISSBasisQuadraticLagrangeInterpolation, &
    & CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err) !Have to do this
  CALL CMISSBasisCreateFinish(QuadraticBasis,Err)

    !Collapsed Quadratic Basis
  CALL CMISSBasisTypeInitialise(QuadraticCollapsedBasis,Err)
  CALL CMISSBasisCreateStart(QuadraticCollapsedBasisUserNumber,QuadraticCollapsedBasis,Err)
  CALL CMISSBasisTypeSet(QuadraticCollapsedBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(QuadraticCollapsedBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(QuadraticCollapsedBasis,(/CMISSBasisQuadraticLagrangeInterpolation, &
       & CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation/),Err)
  CALL CMISSBasisCollapsedXiSet(QuadraticCollapsedBasis,(/CMISSBasisXiCollapsed, &
       & CMISSBasisCollapsedAtXi0,CMISSBasisNotCollapsed/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticCollapsedBasis, &
       & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)  
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(QuadraticCollapsedBasis,.true.,Err) !Have to do this
  CALL CMISSBasisCreateFinish(QuadraticCollapsedBasis,Err)

    !Linear Basis
  CALL CMISSBasisTypeInitialise(LinearBasis,Err)
  CALL CMISSBasisCreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(LinearBasis, &
    & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL CMISSBasisCreateFinish(LinearBasis,Err)

    !Collapsed Linear Basis
  CALL CMISSBasisTypeInitialise(LinearCollapsedBasis,Err)
  CALL CMISSBasisCreateStart(LinearCollapsedBasisUserNumber,LinearCollapsedBasis,Err)
  CALL CMISSBasisTypeSet(LinearCollapsedBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(LinearCollapsedBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(LinearCollapsedBasis,(/CMISSBasisLinearLagrangeInterpolation, &
       & CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation/),Err)
  CALL CMISSBasisCollapsedXiSet(LinearCollapsedBasis,(/CMISSBasisXiCollapsed,CMISSBasisCollapsedAtXi0,CMISSBasisNotCollapsed/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(LinearCollapsedBasis, &
       & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(LinearCollapsedBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL CMISSBasisCreateFinish(LinearCollapsedBasis,Err)

  !Start the creation of a generated ellipsoid mesh
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up an ellipsoid mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshEllipsoidMeshType,Err)
  !Set the quadratic and linear bases
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,[QuadraticBasis,QuadraticCollapsedBasis,LinearBasis,LinearCollapsedBasis],Err)
  !Define the mesh on the region
  CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/LONG_AXIS,SHORT_AXIS,WALL_THICKNESS,CUTOFF_ANGLE/),Err)
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS/),Err)
  
  !Finish the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecompositionCalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry)
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)  
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

  !Create a fibre field and attach it to the geometric field  
  CALL CMISSFieldTypeInitialise(FibreField,Err)
  CALL CMISSFieldCreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSFieldTypeSet(FibreField,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreField,CMISSFieldUVariableType,FieldFibreNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(FibreField,Err)

  !Set Fibre directions (this block is parallel-untested)
  node_idx=0  
  !This is valid only for quadratic basis functions
  TOTAL_NUMBER_NODES_XI(1)=NUMBER_GLOBAL_X_ELEMENTS*2
  TOTAL_NUMBER_NODES_XI(2)=NUMBER_GLOBAL_Y_ELEMENTS*2+1
  TOTAL_NUMBER_NODES_XI(3)=NUMBER_GLOBAL_Z_ELEMENTS*2+1
   
  !Map the correct node number (cn) to geometric node number (gn) G(gn)=cn
  NumberOfCornerNodes=(NUMBER_GLOBAL_Z_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS*NUMBER_GLOBAL_X_ELEMENTS+1)
  TotalNumberOfNodes=(TOTAL_NUMBER_NODES_XI(3))*((TOTAL_NUMBER_NODES_XI(2)-1)*TOTAL_NUMBER_NODES_XI(1)+1)
  
  ALLOCATE(G(TotalNumberOfNodes),STAT=Err)
    
  CornerNode=0
  !Numbering of not corner nodes starts where Corner nodes end
  NotCornerNode=NumberOfCornerNodes
  gn=0

  DO k=1,TOTAL_NUMBER_NODES_XI(3)
     gn=gn+1
     IF (mod(k,2)==0) THEN
        !Not corner node
        NotCornerNode=NotCornerNode+1
        G(gn)=NotCornerNode
     ELSE
        !Corner node
        CornerNode=CornerNode+1
        G(gn)=CornerNode
     ENDIF
     DO j=2,TOTAL_NUMBER_NODES_XI(2)  
        DO i=1,TOTAL_NUMBER_NODES_XI(1)
           gn=gn+1
           IF ((mod(i,2)==0).OR.(mod(j,2)==0).OR.(mod(k,2)==0)) THEN
              !Not corner node
              NotCornerNode=NotCornerNode+1
              G(gn)=NotCornerNode
           ELSE
              !Corner node
              CornerNode=CornerNode+1
              G(gn)=CornerNode
           ENDIF
        ENDDO
     ENDDO
  ENDDO


  XI2delta=(PI-CUTOFF_ANGLE)/(TOTAL_NUMBER_NODES_XI(2)-1)
  XI3=0
  XI3delta=(1.0)/(TOTAL_NUMBER_NODES_XI(3)-1)
  zero=0
  DO k=1, TOTAL_NUMBER_NODES_XI(3)
     !Apex nodes
     j=1
     i=1
     node_idx=node_idx+1
     CorrectNodeNumber=G(node_idx)
     FibreFieldAngle=(/zero,zero,zero/) 
     DO component_idx=1,FieldFibreNumberOfComponents
        CALL CMISSFieldParameterSetUpdateNode(FibreField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,DerivativeUserNumber, &
             & CorrectNodeNumber,component_idx,FibreFieldAngle(component_idx),Err)
     ENDDO
     theta=atan((FIBRE_SLOPE_EPI-FIBRE_SLOPE_ENDO)*XI3+FIBRE_SLOPE_ENDO)
     DO j=2, TOTAL_NUMBER_NODES_XI(2) 
        nu=XI2delta*(j-1)
        omega=PI/2+(-2*nu/(PI-CUTOFF_ANGLE)+1)*atan(SHEET_SLOPE_BASE_ENDO)*(-2*XI3+1)
        !nu=PI-XI2delta*(j-1)
        !omega=PI/2+cos(2*nu)*atan(SHEET_SLOPE_BASE_ENDO)*(-2*XI3+1)
        DO i=1, TOTAL_NUMBER_NODES_XI(1)
           node_idx=node_idx+1
           CorrectNodeNumber=G(node_idx)
           FibreFieldAngle=(/theta,zero,omega/)
           DO component_idx=1,FieldFibreNumberOfComponents
              CALL CMISSFieldParameterSetUpdateNode(FibreField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1, &
                   & DerivativeUserNumber, CorrectNodeNumber,component_idx,FibreFieldAngle(component_idx),Err)
           ENDDO
        ENDDO
     ENDDO
     XI3=XI3+XI3delta
  ENDDO
 
  !Create a material field and attach it to the geometric field  
  CALL CMISSFieldTypeInitialise(MaterialField,Err)
  CALL CMISSFieldCreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL CMISSFieldTypeSet(MaterialField,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialField,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(MaterialField,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialField,CMISSFieldUVariableType,FieldMaterialNumberOfComponents,Err)  
  CALL CMISSFieldComponentInterpolationSet(MaterialField,CMISSFieldUVariableType,1,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(MaterialField,CMISSFieldUVariableType,2,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldCreateFinish(MaterialField,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  !CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2.0_CMISSDP,Err)
  !CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,6.0_CMISSDP,Err)
  !Set Costa material parameters
  DO I=1,7
   CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,I,COSTA_PARAMS(I),Err)
  ENDDO
  !Create the equations_set
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSEquationsSetCreateStart(EquationSetUserNumber,Region,FibreField,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetOrthotropicMaterialCostaSubtype,EquationsSetFieldUserNumber,&
    & EquationsSetField,EquationsSet,Err)
  ! CMISSEquationsSetOrthotropicMaterialCostaSubtype
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the dependent field with 2 variables and 4 components (3 displacement, 1 pressure)
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSFieldCreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL CMISSFieldTypeSet(DependentField,CMISSFieldGeneralType,Err)
  CALL CMISSFieldMeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(DependentField,GeometricField,Err)
  CALL CMISSFieldDependentTypeSet(DependentField,CMISSFieldDependentType,Err)
  CALL CMISSFieldNumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldUVariableType,FieldDependentNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldDelUDelNVariableType,FieldDependentNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,4,LinearMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,4,LinearMeshComponentNumber,Err)
  CALL CMISSFieldScalingTypeSet(DependentField,CMISSFieldUnitScaling,Err)
  CALL CMISSFieldCreateFinish(DependentField,Err)

  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)   

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,-14.0_CMISSDP,Err)

  !Define the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemElasticityClass,CMISSProblemFiniteElasticityType, &
    & CMISSProblemNoSubtype,Err)
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverNewtonJacobianCalculationTypeSet(Solver,CMISSSolverNewtonJacobianFDCalculated,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(Solver,CMISSSolverNewtonJacobianAnalyticCalculated,Err)
  CALL CMISSSolverNewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearDirectSolveType,Err)
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)   
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Grab the list of nodes on inner, outer and top surfaces
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshEllipsoidTopSurfaceType,TopSurfaceNodes,TopNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshEllipsoidInnerSurfaceType,InnerSurfaceNodes,InnerNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshEllipsoidOuterSurfaceType,OuterSurfaceNodes,OuterNormalXi,Err)

  ! ASSIGN BOUNDARY CONDITIONS
  !Fix base of the ellipsoid in z direction
  DO NN=1,SIZE(TopSurfaceNodes,1)
    NODE=TopSurfaceNodes(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NODE,3,ZCoord,Err)
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NODE,3, &
        & CMISSBoundaryConditionFixed,ZCoord,Err)
    ENDIF
  ENDDO

  !Apply inner surface pressure
  !NOTE: Surface pressure goes into pressure_values_set_type of the DELUDELN type
  DO NN=1,SIZE(InnerSurfaceNodes,1)
    NODE=InnerSurfaceNodes(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldDelUDelNVariableType,1,1,NODE, &
        & ABS(InnerNormalXi), &
        & CMISSBoundaryConditionPressure,INNER_PRESSURE,Err)
    ENDIF
  ENDDO

  !Apply outer surface pressure
  DO NN=1,SIZE(OuterSurfaceNodes,1)
    NODE=OuterSurfaceNodes(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldDelUDelNVariableType,1,1,NODE, &
        & ABS(OuterNormalXi), &
        & CMISSBoundaryConditionPressure,OUTER_PRESSURE,Err)
    ENDIF
  ENDDO

  !Fix apex in x and y direction
!!$  DO NN=1, NUMBER_GLOBAL_Z_ELEMENTS*2+1
!!$     NODE=NN+((NN-1)*(NUMBER_GLOBAL_X_ELEMENTS*2)*(NUMBER_GLOBAL_Z_ELEMENTS*2+1))
!!$     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
!!$     IF(NodeDomain==ComputationalNodeNumber) THEN
!!$        CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NODE,1,XCoord,Err)
!!$        CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NODE,2,YCoord,Err)
!!$        CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NODE,1, &
!!$             & CMISSBoundaryConditionFixed,XCoord,Err)
!!$        CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NODE,2, &
!!$             & CMISSBoundaryConditionFixed,YCoord,Err)
!!$     ENDIF
!!$  ENDDO

  !Fix more nodes at the base to stop free body motion
  X_FIXED=.FALSE.
  Y_FIXED=.FALSE.
  DO NN=1,SIZE(TopSurfaceNodes,1)
    NODE=TopSurfaceNodes(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NODE,1,XCoord,Err)
      CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NODE,2,YCoord,Err)
      IF(ABS(XCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NODE,1, &
          & CMISSBoundaryConditionFixed,XCoord,Err)
        WRITE(*,*) "FIXING NODE",NODE,"IN X DIRECTION"
        X_FIXED=.TRUE.
      ENDIF
      IF(ABS(YCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NODE,2, &
          & CMISSBoundaryConditionFixed,YCoord,Err)
        WRITE(*,*) "FIXING NODE",NODE,"IN Y DIRECTION"
        Y_FIXED=.TRUE.
      ENDIF
    ENDIF
  ENDDO
  CALL MPI_REDUCE(X_FIXED,X_OKAY,1,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_REDUCE(Y_FIXED,Y_OKAY,1,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_WORLD,MPI_IERROR)
  IF(ComputationalNodeNumber==0) THEN
    IF(.NOT.(X_OKAY.AND.Y_OKAY)) THEN
      WRITE(*,*) "Free body motion could not be prevented!"
      CALL CMISSFinalise(Err)
      STOP
    ENDIF
  ENDIF

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblemSolve(Problem,Err)

  !Output solution  
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"QuadraticEllipsoidCostaReadIn"// PARAMETERSET,"FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"QuadraticEllipsoidCostaReadIn"// PARAMETERSET,"FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM QUADRATICELLIPSOIDCOSTAREADINEXAMPLE

