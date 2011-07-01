!> \file
!> \author Adam Reeve
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
!> Contributor(s): Adam Reeve
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

!> \example FiniteElasticity/Cantilever/src/CantileverExample.f90
!! Example program to solve a finite elasticity equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/Cantilever/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/Cantilever/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM CANTILEVEREXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: Height=40.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: Width=60.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: Length=40.0_CMISSDP
  INTEGER(CMISSIntg), PARAMETER :: DisplacementInterpolationType=CMISSBasisCubicHermiteInterpolation
  INTEGER(CMISSIntg), PARAMETER :: PressureInterpolationType=CMISSBasisLinearLagrangeInterpolation
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussXi=4
  REAL(CMISSDP), PARAMETER :: Density=9.0E-4_CMISSDP !in g mm^-3
  REAL(CMISSDP), PARAMETER :: Gravity(3)=[0.0_CMISSDP,0.0_CMISSDP,-9.8_CMISSDP] !in m s^-2

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DisplacementBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldSourceUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx,component_idx,deriv_idx
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg) :: LeftNormalXi
  REAL(CMISSDP) :: Value

  !CMISS variables
  TYPE(CMISSBasisType) :: DisplacementBasis,PressureBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,FibreField,MaterialField,DependentField,SourceField,EquationsSetField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSControlLoopType) :: ControlLoop

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

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberGlobalXElements=2
  NumberGlobalYElements=1
  NumberGlobalZElements=1
  NumberOfDomains=NumberOfComputationalNodes

  !Create a 3D rectangular cartesian coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !Define basis function for displacement
  CALL CMISSBasisTypeInitialise(DisplacementBasis,Err)
  CALL CMISSBasisCreateStart(DisplacementBasisUserNumber,DisplacementBasis,Err)
  SELECT CASE(DisplacementInterpolationType)
  CASE(1,2,3,4)
    CALL CMISSBasisTypeSet(DisplacementBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CASE(7,8,9)
    CALL CMISSBasisTypeSet(DisplacementBasis,CMISSBasisSimplexType,Err)
  END SELECT
  CALL CMISSBasisNumberOfXiSet(DisplacementBasis,3,Err)
  CALL CMISSBasisInterpolationXiSet(DisplacementBasis,[DisplacementInterpolationType,DisplacementInterpolationType, &
      & DisplacementInterpolationType],Err)
  IF(NumberOfGaussXi>0) THEN
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(DisplacementBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL CMISSBasisCreateFinish(DisplacementBasis,Err)

  !Basis for pressure
  CALL CMISSBasisTypeInitialise(PressureBasis,Err)
  CALL CMISSBasisCreateStart(PressureBasisUserNumber,PressureBasis,Err)
  SELECT CASE(PressureInterpolationType)
  CASE(1,2,3,4)
    CALL CMISSBasisTypeSet(PressureBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CASE(7,8,9)
    CALL CMISSBasisTypeSet(PressureBasis,CMISSBasisSimplexType,Err)
  END SELECT
  CALL CMISSBasisNumberOfXiSet(PressureBasis,3,Err)
  CALL CMISSBasisInterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType, &
      & PressureInterpolationType],Err)
  IF(NumberOfGaussXi>0) THEN
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL CMISSBasisCreateFinish(PressureBasis,Err)

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,[DisplacementBasis,PressureBasis],Err)
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[Width,Height],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[Width,Length,Height],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
      & NumberGlobalZElements],Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldVariableLabelSet(GeometricField,CMISSFieldUVariableType,"Geometry",Err)
  CALL CMISSFieldScalingTypeSet(GeometricField,CMISSFieldArithmeticMeanScaling,Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

  !Create a fibre field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(FibreField,Err)
  CALL CMISSFieldCreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSFieldTypeSet(FibreField,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreField,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSFieldVariableLabelSet(FibreField,CMISSFieldUVariableType,"Fibre",Err)
  CALL CMISSFieldCreateFinish(FibreField,Err)

  !Create the equations_set
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationSetUserNumber,Region,FibreField,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetMooneyRivlinSubtype,EquationsSetFieldUserNumber,EquationsSetField, &
    & EquationsSet,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the dependent field
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSFieldVariableLabelSet(DependentField,CMISSFieldUVariableType,"Dependent",Err)
  DO component_idx=1,3
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,component_idx,1,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,component_idx,1,Err)
  ENDDO
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,4,2,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,4,2,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldUVariableType,4,CMISSFieldNodeBasedInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentField,CMISSFieldDelUDelNVariableType,4,CMISSFieldNodeBasedInterpolation,Err)
  CALL CMISSFieldScalingTypeSet(DependentField,CMISSFieldArithmeticMeanScaling,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL CMISSFieldTypeInitialise(MaterialField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL CMISSFieldVariableLabelSet(MaterialField,CMISSFieldUVariableType,"Material",Err)
  CALL CMISSFieldVariableLabelSet(MaterialField,CMISSFieldVVariableType,"Density",Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,6.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldVVariableType,CMISSFieldValuesSetType,1,Density,Err)

  !Create the source field with the gravity vector
  CALL CMISSFieldTypeInitialise(SourceField,Err)
  CALL CMISSEquationsSetSourceCreateStart(EquationsSet,FieldSourceUserNumber,SourceField,Err)
  CALL CMISSEquationsSetSourceCreateFinish(EquationsSet,Err)
  DO component_idx=1,3
    CALL CMISSFieldComponentValuesInitialise(SourceField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
        & component_idx,Gravity(component_idx),Err)
  ENDDO

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsFullMatrices,Err)
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
  CALL CMISSControlLoopMaximumIterationsSet(ControlLoop,2,Err)
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(Solver,CMISSSolverNewtonJacobianFDCalculated,Err)
  CALL CMISSSolverNewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearIterativeSolveType,Err)
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshRegularLeftSurface,LeftSurfaceNodes,LeftNormalXi,Err)

  !Fix x=0 nodes
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      DO component_idx=1,3
        CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NodeNumber, &
          & component_idx,Value,Err)
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NodeNumber, &
          & component_idx, &
          & CMISSBoundaryConditionFixed,Value,Err)
        DO deriv_idx=3,8
          CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,deriv_idx, &
              & NodeNumber,component_idx,Value,Err)
          CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,deriv_idx,NodeNumber, &
            & component_idx, &
            & CMISSBoundaryConditionFixed,Value,Err)
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblemSolve(Problem,Err)

  !Output solution
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"Cantilever","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"Cantilever","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM CANTILEVEREXAMPLE

