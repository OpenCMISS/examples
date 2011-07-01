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
PROGRAM CUBOIDGENERICEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters
  CHARACTER(LEN=100),ALLOCATABLE :: ALL_ARGS(:)
  INTEGER(CMISSIntg) :: Basis1,Basis2
  REAL(CMISSDP) :: X_LENG,Y_LENG,Z_LENG   ! dimension of cuboid
  INTEGER(CMISSIntg)  ::   NumberGlobalXElements
  INTEGER(CMISSIntg)  ::   NumberGlobalYElements
  INTEGER(CMISSIntg)  ::   NumberGlobalZElements
  INTEGER(CMISSIntg)  ::   incrementSteps
  REAL(CMISSDP),ALLOCATABLE :: BC(:,:)
  REAL(CMISSDP), PARAMETER :: C1=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: C2=6.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: CubicBasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: DisplacementMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureMeshComponentNumber=2

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponents=2

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariables=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponents=4

  !INTEGER(CMISSIntg), PARAMETER :: FieldAnalyticUserNumber=1337

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program types


  !Program variables
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  !CMISS variables

  TYPE(CMISSBasisType) :: CubicBasis, QuadraticBasis, LinearBasis, Bases(2)
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,FibreField,MaterialField
  TYPE(CMISSFieldType) :: DependentField,EquationsSetField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSControlLoopType) :: ControlLoop

  !Other variables
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face1Nodes(:),Face2Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face3Nodes(:),Face4Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face5Nodes(:),Face6Nodes(:)
  INTEGER(CMISSIntg),POINTER :: Nodes(:)
  INTEGER(CMISSIntg) :: FaceXi(6)
  INTEGER(CMISSIntg) :: I,VariableType

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

  !Parse command line options: basis_1, basis_2, num_els_x, .. z, x_leng,.. z_leng, BC
  CALL GET_ALL_ARGS(ALL_ARGS)
  CALL GET_ARGUMENT(ALL_ARGS,"basis_displacement",ARG_INT=Basis1)
  CALL GET_ARGUMENT(ALL_ARGS,"basis_pressure",ARG_INT=Basis2)
  CALL GET_ARGUMENT(ALL_ARGS,"NumElsX",ARG_INT=NumberGlobalXElements)
  CALL GET_ARGUMENT(ALL_ARGS,"NumElsY",ARG_INT=NumberGlobalYElements)
  CALL GET_ARGUMENT(ALL_ARGS,"NumElsZ",ARG_INT=NumberGlobalZElements)
  CALL GET_ARGUMENT(ALL_ARGS,"XLeng",ARG_DP=X_LENG)
  CALL GET_ARGUMENT(ALL_ARGS,"YLeng",ARG_DP=Y_LENG)
  CALL GET_ARGUMENT(ALL_ARGS,"ZLeng",ARG_DP=Z_LENG)
  CALL GET_ARGUMENT(ALL_ARGS,"IncrementSteps",ARG_INT=incrementSteps)
  CALL GET_BC(ALL_ARGS,BC)

  !Set all diganostic levels on for testing
  CALL CMISSDiagnosticsSetOn(CMISSFromDiagType,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_FINITE_ELEMENT_CALCULATE"/),Err)
  !CALL CMISSDiagnosticsSetOff(Err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  write(*,*) "NumberOfDomains=",NumberOfComputationalNodes
  NumberOfDomains=NumberOfComputationalNodes

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR) ! NOW A PARAMETER
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR) ! NOW A PARAMETER
  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR) ! NOW A PARAMETER 
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  
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

  !Define basis functions
  CALL CMISSBasisTypeInitialise(LinearBasis,Err)
  CALL CMISSBasisCreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(LinearBasis, &
    & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err)
  CALL CMISSBasisCreateFinish(LinearBasis,Err)

  CALL CMISSBasisTypeInitialise(QuadraticBasis,Err)
  CALL CMISSBasisCreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasisInterpolationXiSet(QuadraticBasis,(/CMISSBasisQuadraticLagrangeInterpolation, &
    & CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err)
  CALL CMISSBasisCreateFinish(QuadraticBasis,Err)

  CALL CMISSBasisTypeInitialise(CubicBasis,Err)
  CALL CMISSBasisCreateStart(CubicBasisUserNumber,CubicBasis,Err)
  CALL CMISSBasisInterpolationXiSet(CubicBasis,(/CMISSBasisCubicLagrangeInterpolation, &
    & CMISSBasisCubicLagrangeInterpolation,CMISSBasisCubicLagrangeInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(CubicBasis, &
    & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(CubicBasis,.true.,Err) !Enable 3D interpolation on faces
  CALL CMISSBasisCreateFinish(CubicBasis,Err)

  SELECT CASE(Basis1)
  CASE(2)
    Bases(1)=QuadraticBasis
  CASE(3)
    Bases(1)=CubicBasis
  CASE DEFAULT
    WRITE(*,*) "Invalid basis type for displacement"
    STOP
  END SELECT

  SELECT CASE(Basis2)
  CASE(1)
    Bases(2)=LinearBasis
  CASE(2)
    Bases(2)=QuadraticBasis
  CASE DEFAULT
    WRITE(*,*) "Invalid basis type for pressure"
    STOP
  END SELECT

  !Start the creation of a generated cylinder mesh
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up an regular (cuboid) mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the quadratic basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Bases,Err)
  !Define the mesh on the region
  CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[X_LENG,Y_LENG,Z_LENG],Err)
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements],Err)
  
  !Finish the creation of generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSRandomSeedsSet(0_CMISSIntg,Err) !To keep the automatic decomposition same each time
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Automatic decomposition
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  !Manual decomposition
!   IF(NumberOfDomains>1) THEN
!     CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionUserDefinedType,Err)
!     !Set all elements but last one to first domain
!     CALL CMISSMeshNumberOfElementsGet(Mesh,NE,Err)
!     do E=1,NE/2
!       CALL CMISSDecompositionElementDomainSet(Decomposition,E,0,Err)
!     enddo
!     do E=NE/2+1,NE
!       CALL CMISSDecompositionElementDomainSet(Decomposition,E,1,Err)
!     enddo
!     CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
!   ENDIF
  CALL CMISSDecompositionCalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry)
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)  
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,DisplacementMeshComponentNumber,Err)
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
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,1,PressureMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,2,PressureMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,3,PressureMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(FibreField,Err)

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
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,C1,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,C2,Err)

  !Create the equations_set
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSEquationsSetCreateStart(EquationSetUserNumber,Region,FibreField,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetMooneyRivlinSubtype,EquationsSetFieldUserNumber,EquationsSetField,&
    & EquationsSet,Err)
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
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,1,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,2,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,3,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,4,PressureMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,1,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,2,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,3,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,4,PressureMeshComponentNumber,Err)
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
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  CALL CMISSControlLoopMaximumIterationsSet(ControlLoop,incrementSteps,Err)  ! this one sets the increment loop counter
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)
  
  !Create the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
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

  !BC Assignment
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Get surfaces - will fix two nodes on bottom face, pressure conditions inside
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularFrontSurface,Face1Nodes,FaceXi(1),Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularBackSurface,Face2Nodes,FaceXi(2),Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularRightSurface,Face3Nodes,FaceXi(3),Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularLeftSurface,Face4Nodes,FaceXi(4),Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularTopSurface,Face5Nodes,FaceXi(5),Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularBottomSurface,Face6Nodes,FaceXi(6),Err)

  DO I=1,SIZE(BC,1)
    IF(INT(BC(I,3))==CMISSBoundaryConditionPressure.OR.INT(BC(I,3))==CMISSBoundaryConditionPressureIncremented) THEN
      VariableType=CMISSFieldDelUDelNVariableType
    ELSE
      VariableType=CMISSFieldUVariableType
    ENDIF
    SELECT CASE (INT(BC(I,1)))
    CASE (1)
      Nodes=>Face1Nodes
    CASE (2)
      Nodes=>Face2Nodes
    CASE (3)
      Nodes=>Face3Nodes
    CASE (4)
      Nodes=>Face4Nodes
    CASE (5)
      Nodes=>Face5Nodes
    CASE (6)
      Nodes=>Face6Nodes
    END SELECT
    CALL SET_BC(Decomposition,DependentField,GeometricField,BoundaryConditions,VariableType,INT(BC(I,3)),Nodes,INT(BC(I,2)), &
      & BC(I,4),ComputationalNodeNumber)
  ENDDO

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblemSolve(Problem,Err)

  !Output solution  
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"outputs/CuboidGeneric","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"outputs/CuboidGeneric","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  CALL CMISSFinalise(Err)
  IF(ALLOCATED(BC)) DEALLOCATE(BC)
  IF(ALLOCATED(Face1Nodes)) DEALLOCATE(Face1Nodes)
  IF(ALLOCATED(Face2Nodes)) DEALLOCATE(Face2Nodes)
  IF(ALLOCATED(Face3Nodes)) DEALLOCATE(Face3Nodes)
  IF(ALLOCATED(Face4Nodes)) DEALLOCATE(Face4Nodes)
  IF(ALLOCATED(Face5Nodes)) DEALLOCATE(Face5Nodes)
  IF(ALLOCATED(Face6Nodes)) DEALLOCATE(Face6Nodes)
  IF(ALLOCATED(ALL_ARGS)) DEALLOCATE(ALL_ARGS)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

  CONTAINS

  !-------------------------------------------------------------------
  !> Converts a string to upper case
  SUBROUTINE UPPER_CASE(LWORD,UWORD)
    CHARACTER(LEN=*),INTENT(IN) :: LWORD
    CHARACTER(LEN=LEN(LWORD)),INTENT(OUT) :: UWORD
    INTEGER(CMISSIntg) ::I,IC,NLEN

    NLEN = LEN_TRIM(LWORD)
    UWORD=LWORD
    DO I=1,NLEN
      IC = ICHAR(LWORD(I:I))
      IF (IC >= 97 .and. IC <= 122) UWORD(I:I) = CHAR(IC-32)
    ENDDO
  END SUBROUTINE UPPER_CASE

  !-------------------------------------------------------------------
  !> Just grabs the whole input file
  SUBROUTINE GET_ALL_ARGS(ALL_ARGS)
    CHARACTER(LEN=100),ALLOCATABLE,INTENT(INOUT) :: ALL_ARGS(:)
    !Local variables
    CHARACTER(LEN=100) :: ARGS(100)
    INTEGER(CMISSIntg) :: I,REASON
    CHARACTER(LEN=100) :: WORD

    ARGS=""
    I=0
    DO
      I=I+1
      IF(I>100) THEN
        WRITE(*,*) "Increase ARGS size"
        STOP
      ENDIF
      READ(5,'(A)',IOSTAT=Reason) WORD
      IF(Reason>0)  THEN
        WRITE(*,*) "Error reading options file"
        STOP
      ELSEIF(Reason<0) THEN
        EXIT
      ELSE
        CALL UPPER_CASE(WORD,ARGS(I))
      ENDIF
    ENDDO

    ALLOCATE(ALL_ARGS(I))
    ALL_ARGS=ARGS(1:I)

  END SUBROUTINE GET_ALL_ARGS

  !-------------------------------------------------------------------
  !> Returns the argument of the requested type as a varying string
  SUBROUTINE GET_ARGUMENT(ALL_ARGS,ARG_TYPE,ARG,ARG_INT,ARG_DP)
    CHARACTER(LEN=*), INTENT(IN) :: ALL_ARGS(:)
    CHARACTER(LEN=*), INTENT(IN) :: ARG_TYPE
    CHARACTER(LEN=256), INTENT(OUT), OPTIONAL :: ARG
    INTEGER(CMISSIntg), INTENT(OUT), OPTIONAL :: ARG_INT
    REAL(CMISSDP), INTENT(OUT), OPTIONAL :: ARG_DP
    !Local variables
    CHARACTER(LEN=LEN_TRIM(ARG_TYPE)) :: ARG_TYPE_UPPER
    CHARACTER(LEN=256) :: WORD_UPPER,ARGOUT
    INTEGER(CMISSIntg) :: NARGS,I,LENG,WORD_LENG

    NARGS=SIZE(ALL_ARGS,1)
    LENG=LEN_TRIM(ARG_TYPE) !STRING TO LOOK FOR
    CALL UPPER_CASE(ARG_TYPE,ARG_TYPE_UPPER)

    DO I=1,NARGS
      WORD_UPPER=ALL_ARGS(I)
      WORD_LENG=LEN_TRIM(WORD_UPPER)
      IF(WORD_UPPER(1:1+LENG)=="-"//TRIM(ARG_TYPE_UPPER)) THEN
        IF(WORD_UPPER(2+LENG:2+LENG)=="=") THEN
          ! USING = AS DELIMITER
          ARGOUT=WORD_UPPER(3+LENG:WORD_LENG)
          EXIT
        ELSE
          ! USING A SPACE AS DELIMITER
          ARGOUT=TRIM(WORD_UPPER(3+LENG:WORD_LENG))
          EXIT
        ENDIF
      ENDIF
    ENDDO

    IF(PRESENT(ARG)) CALL UPPER_CASE(ARGOUT,ARG)
    IF(PRESENT(ARG_INT)) READ(ARGOUT,*) ARG_INT
    IF(PRESENT(ARG_DP)) READ(ARGOUT,*) ARG_DP
  END SUBROUTINE GET_ARGUMENT

  !-------------------------------------------------------------------
  !> Sets a given boundary condition to a group of nodes
  SUBROUTINE SET_BC(Decomposition,DependentField,GeometricField,BoundaryConditions,VariableType,BCType,Nodes,Component,Value, &
    & ComputationalNodeNumber)
    TYPE(CMISSDecompositionType),INTENT(IN) :: Decomposition
    TYPE(CMISSFieldType),INTENT(IN) :: DependentField
    TYPE(CMISSFieldType),INTENT(IN) :: GeometricField    
    TYPE(CMISSBoundaryConditionsType),INTENT(INOUT) :: BoundaryConditions
    INTEGER(CMISSIntg),INTENT(IN) :: VariableType
    INTEGER(CMISSIntg),INTENT(IN) :: BCType
    INTEGER(CMISSIntg),INTENT(IN) :: Nodes(:)
    INTEGER(CMISSIntg),INTENT(IN) :: Component
    REAL(CMISSDP),INTENT(IN) :: Value
    INTEGER(CMISSIntg),INTENT(IN) :: ComputationalNodeNumber
    !Local variables
    INTEGER(CMISSIntg) :: i,node,NodeDomain,BCType2
    REAL(CMISSDP) :: coord

    DO I=1,SIZE(Nodes)
      node=Nodes(I)
      CALL CMISSDecompositionNodeDomainGet(Decomposition,node,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(BCType<0) THEN
          !Displacement type condition - get the undeformed geometric value first
          CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node, &
            & Component,coord,Err)
          coord=coord+Value
          IF(BCType==-777) BCType2=CMISSBoundaryConditionFixed
          IF(BCType==-888) BCType2=CMISSBoundaryConditionFixedIncremented
          CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,VariableType,1,1,node,Component,BCType2,coord,Err)
        ELSE
          CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,VariableType,1,1,node,Component,BCType,Value,Err)
        ENDIF
      ENDIF
    ENDDO
    
  END SUBROUTINE SET_BC

  !-------------------------------------------------------------------
  !> Grabs all specified BC from the command line
  SUBROUTINE GET_BC(ALL_ARGS,BC)
    CHARACTER(LEN=*), INTENT(IN) :: ALL_ARGS(:)
    REAL(CMISSDP),ALLOCATABLE,INTENT(OUT) :: BC(:,:)
    !Local variable
    REAL(CMISSDP) :: BCtemp(100,4)
    INTEGER(CMISSIntg) :: NARG,N,I,J,POS1,POS2
    CHARACTER(LEN=100) :: WORD_UPPER
    CHARACTER(LEN=100) :: WORDS(10)
    
    N=0

    !Find out how many BCs there are
    NARG=SIZE(ALL_ARGS,1)
    DO I=1,NARG
      WORD_UPPER=ALL_ARGS(I)
      IF(WORD_UPPER(1:3)=="-BC") THEN
        N=N+1
        !Split the string into words
        POS1=1
        J=0
        DO
          POS2 = INDEX(TRIM(WORD_UPPER(POS1:)), " ")
          IF(POS2 == 0) THEN
            J=J+1
            WORDS(J) = TRIM(WORD_UPPER(POS1:))
            EXIT
          ENDIF
          J=J+1
          WORDS(J) = WORD_UPPER(POS1:POS1+POS2-2)
          POS1=POS2+POS1
        ENDDO
        ! surface type
        SELECT CASE (TRIM(WORDS(2)))
        CASE ("FRONT")
          BCtemp(N,1)=1
        CASE ("BACK")
          BCtemp(N,1)=2
        CASE ("RIGHT")
          BCtemp(N,1)=3
        CASE ("LEFT")
          BCtemp(N,1)=4
        CASE ("TOP")
          BCtemp(N,1)=5
        CASE ("BOTTOM")
          BCtemp(N,1)=6
        CASE DEFAULT
          WRITE(*,*) "Invalid type of surface in bc specification"
          STOP
        END SELECT
        ! component
        READ(WORDS(3),*) BCtemp(N,2)
        ! bc type
        WORD_UPPER=WORDS(4)
        SELECT CASE (TRIM(WORD_UPPER))
        CASE ("DISPLACEMENT")
          BCtemp(N,3)=-777
        CASE ("DISPLACEMENTINCREMENTED")
          BCtemp(N,3)=-888
        CASE ("FIXED")
          BCtemp(N,3)=CMISSBoundaryConditionFixed
        CASE ("FIXEDINCREMENTED")
          BCtemp(N,3)=CMISSBoundaryConditionFixedIncremented
        CASE ("PRESSURE")
          BCtemp(N,3)=CMISSBoundaryConditionPressure
        CASE ("PRESSUREINCREMENTED")
          BCtemp(N,3)=CMISSBoundaryConditionPressureIncremented
        CASE DEFAULT
          WRITE(*,*) "Invalid type of bc in bc specification"
          STOP
        END SELECT
        ! value
        WORD_UPPER=WORDS(5)
        READ(WORD_UPPER,*) BCtemp(N,4)
      ENDIF
    ENDDO

    ALLOCATE(BC(N,4)) ! face_no,component,bctype,value
    BC=BCtemp(1:N,:)
    
  END SUBROUTINE GET_BC


END PROGRAM CUBOIDGENERICEXAMPLE

