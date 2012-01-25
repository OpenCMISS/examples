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

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_FINITE_ELEMENT_CALCULATE"/),Err)

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
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_TypeSet(CoordinateSystem,CMISS_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystem_OriginSet(CoordinateSystem,(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/),Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the CS to the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange, each with collapsed variant
    !Quadratic Basis
  CALL CMISSBasis_Initialise(QuadraticBasis,Err)
  CALL CMISSBasis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasis_InterpolationXiSet(QuadraticBasis,(/CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & (/CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME/),Err)
  CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err) !Have to do this
  CALL CMISSBasis_CreateFinish(QuadraticBasis,Err)

    !Collapsed Quadratic Basis
  CALL CMISSBasis_Initialise(QuadraticCollapsedBasis,Err)
  CALL CMISSBasis_CreateStart(QuadraticCollapsedBasisUserNumber,QuadraticCollapsedBasis,Err)
  CALL CMISSBasis_TypeSet(QuadraticCollapsedBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(QuadraticCollapsedBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(QuadraticCollapsedBasis,(/CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
       & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION/),Err)
  CALL CMISSBasis_CollapsedXiSet(QuadraticCollapsedBasis,(/CMISS_BASIS_XI_COLLAPSED, &
       & CMISS_BASIS_COLLAPSED_AT_XI0,CMISS_BASIS_NOT_COLLAPSED/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticCollapsedBasis, &
       & (/CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME/),Err)  
  CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(QuadraticCollapsedBasis,.true.,Err) !Have to do this
  CALL CMISSBasis_CreateFinish(QuadraticCollapsedBasis,Err)

    !Linear Basis
  CALL CMISSBasis_Initialise(LinearBasis,Err)
  CALL CMISSBasis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearBasis, &
    & (/CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME/),Err)
  CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL CMISSBasis_CreateFinish(LinearBasis,Err)

    !Collapsed Linear Basis
  CALL CMISSBasis_Initialise(LinearCollapsedBasis,Err)
  CALL CMISSBasis_CreateStart(LinearCollapsedBasisUserNumber,LinearCollapsedBasis,Err)
  CALL CMISSBasis_TypeSet(LinearCollapsedBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(LinearCollapsedBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(LinearCollapsedBasis,(/CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
       & CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION/),Err)
  CALL CMISSBasis_CollapsedXiSet(LinearCollapsedBasis,(/CMISS_BASIS_XI_COLLAPSED,CMISS_BASIS_COLLAPSED_AT_XI0, &
    & CMISS_BASIS_NOT_COLLAPSED/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearCollapsedBasis, &
       & (/CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME/),Err)
  CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(LinearCollapsedBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL CMISSBasis_CreateFinish(LinearCollapsedBasis,Err)

  !Start the creation of a generated ellipsoid mesh
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up an ellipsoid mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_ELLIPSOID_MESH_TYPE,Err)
  !Set the quadratic and linear bases
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,[QuadraticBasis,QuadraticCollapsedBasis,LinearBasis,LinearCollapsedBasis],Err)
  !Define the mesh on the region
  CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/LONG_AXIS,SHORT_AXIS,WALL_THICKNESS,CUTOFF_ANGLE/),Err)
  CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS/),Err)
  
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry)
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)  
  CALL CMISSField_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Create a fibre field and attach it to the geometric field  
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(FibreField,Err)

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
        CALL CMISSField_ParameterSetUpdateNode(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
          & DerivativeUserNumber, &
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
              CALL CMISSField_ParameterSetUpdateNode(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
                   & DerivativeUserNumber, CorrectNodeNumber,component_idx,FibreFieldAngle(component_idx),Err)
           ENDDO
        ENDDO
     ENDDO
     XI3=XI3+XI3delta
  ENDDO
 
  !Create a material field and attach it to the geometric field  
  CALL CMISSField_Initialise(MaterialField,Err)
  CALL CMISSField_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL CMISSField_TypeSet(MaterialField,CMISS_FIELD_MATERIAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(MaterialField,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(MaterialField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)  
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_CreateFinish(MaterialField,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  !CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2.0_CMISSDP,Err)
  !CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,6.0_CMISSDP,Err)
  !Set Costa material parameters
  DO I=1,7
   CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,I, &
     & COSTA_PARAMS(I),Err)
  ENDDO
  !Create the equations_set
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSEquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE, &
      & EquationsSetFieldUserNumber,&
    & EquationsSetField,EquationsSet,Err)
  ! CMISS_EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field with 2 variables and 4 components (3 displacement, 1 pressure)
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSField_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL CMISSField_TypeSet(DependentField,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL CMISSField_DependentTypeSet(DependentField,CMISS_FIELD_DEPENDENT_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL CMISSField_ScalingTypeSet(DependentField,CMISS_FIELD_UNIT_SCALING,Err)
  CALL CMISSField_CreateFinish(DependentField,Err)

  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)   

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4, &
    & -14.0_CMISSDP, &
    & Err)

  !Define the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_ELASTICITY_CLASS,CMISS_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMISS_PROBLEM_NO_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)   
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Grab the list of nodes on inner, outer and top surfaces
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_ELLIPSOID_TOP_SURFACE,TopSurfaceNodes,TopNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_ELLIPSOID_INNER_SURFACE,InnerSurfaceNodes,InnerNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_ELLIPSOID_OUTER_SURFACE,OuterSurfaceNodes,OuterNormalXi,Err)

  ! ASSIGN BOUNDARY CONDITIONS
  !Fix base of the ellipsoid in z direction
  DO NN=1,SIZE(TopSurfaceNodes,1)
    NODE=TopSurfaceNodes(NN)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE,3, &
        & ZCoord, &
        & Err)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED,ZCoord,Err)
    ENDIF
  ENDDO

  !Apply inner surface pressure
  !NOTE: Surface pressure goes into pressure_values_set_type of the DELUDELN type
  DO NN=1,SIZE(InnerSurfaceNodes,1)
    NODE=InnerSurfaceNodes(NN)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NODE, &
        & ABS(InnerNormalXi), &
        & CMISS_BOUNDARY_CONDITION_PRESSURE,INNER_PRESSURE,Err)
    ENDIF
  ENDDO

  !Apply outer surface pressure
  DO NN=1,SIZE(OuterSurfaceNodes,1)
    NODE=OuterSurfaceNodes(NN)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NODE, &
        & ABS(OuterNormalXi), &
        & CMISS_BOUNDARY_CONDITION_PRESSURE,OUTER_PRESSURE,Err)
    ENDIF
  ENDDO

  !Fix apex in x and y direction
!!$  DO NN=1, NUMBER_GLOBAL_Z_ELEMENTS*2+1
!!$     NODE=NN+((NN-1)*(NUMBER_GLOBAL_X_ELEMENTS*2)*(NUMBER_GLOBAL_Z_ELEMENTS*2+1))
!!$     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
!!$     IF(NodeDomain==ComputationalNodeNumber) THEN
!!$        CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE,1,XCoord,Err)
!!$        CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE,2,YCoord,Err)
!!$        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,1, &
!!$             & CMISS_BOUNDARY_CONDITION_FIXED,XCoord,Err)
!!$        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,2, &
!!$             & CMISS_BOUNDARY_CONDITION_FIXED,YCoord,Err)
!!$     ENDIF
!!$  ENDDO

  !Fix more nodes at the base to stop free body motion
  X_FIXED=.FALSE.
  Y_FIXED=.FALSE.
  DO NN=1,SIZE(TopSurfaceNodes,1)
    NODE=TopSurfaceNodes(NN)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE,1, &
        & XCoord, &
        & Err)
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE,2, &
        & YCoord, &
        & Err)
      IF(ABS(XCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,1, &
          & CMISS_BOUNDARY_CONDITION_FIXED,XCoord,Err)
        WRITE(*,*) "FIXING NODE",NODE,"IN X DIRECTION"
        X_FIXED=.TRUE.
      ENDIF
      IF(ABS(YCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,2, &
          & CMISS_BOUNDARY_CONDITION_FIXED,YCoord,Err)
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

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Output solution  
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"QuadraticEllipsoidCostaReadIn"// PARAMETERSET,"FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"QuadraticEllipsoidCostaReadIn"// PARAMETERSET,"FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM QUADRATICELLIPSOIDCOSTAREADINEXAMPLE

