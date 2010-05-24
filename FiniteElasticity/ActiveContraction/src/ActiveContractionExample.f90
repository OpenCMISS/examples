!> \file
!> $Id: ActiveContractionExample.f90 20 2007-05-28 20:22:52Z cpb $
!> \author Sander Land
!> \brief This is an example program to solve active contraction based finite elasticity equation using openCMISS calls.
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
!> Contributor(s): Code based on the LargeUniAxialExtensionExample by Chris Bradley
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


!> Main program
PROGRAM ActiveContractionExample

  USE OPENCMISS
  USE MPI

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=1.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldGPUserNumber=5 ! temp/test
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumber=6


  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  REAL(CMISSDP), PARAMETER :: START_TIME = 0.0, END_TIME = 10.0, DT = 1.0  ! ms

  LOGICAL, PARAMETER :: TEST_GAUSS_POINT_FIELD = .FALSE.

  !Program types

  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: D, E, I, y_idx,z_idx,NodeNumber1,Node1Domain

  REAL(CMISSDP) :: TMP

  !CMISS variables
  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,FibreField,MaterialField,DependentField, GPfield, IndependentField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

  TYPE(CMISSControlLoopType) :: ControlLoop

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: Err

  !Intialise cmiss
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  !Set all diganostic levels on for testing
  !CALL CMISSDiagnosticsSetOn(CMISSFromDiagType,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_RESIDUAL_EVALUATE"/),Err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberGlobalXElements=1
  NumberGlobalYElements=1
  NumberGlobalZElements=1  
  NumberOfDomains=NumberOfComputationalNodes

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
  !CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  !CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  !CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  !CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Create a 3D rectangular cartesian coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !Define basis function - tri-linear Lagrange  
  CALL CMISSBasisTypeInitialise(Basis,Err)
  CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err) 
  CALL CMISSBasisCreateFinish(Basis,Err)

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/WIDTH,HEIGHT,LENGTH/),Err)
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NumberGlobalXElements,NumberGlobalYElements, &
      & NumberGlobalZElements/),Err)

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
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

  !Create a fibre field and attach it to the geometric field  
  CALL CMISSFieldTypeInitialise(FibreField,Err)
  CALL CMISSFieldCreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSFieldTypeSet(FibreField,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSFieldCreateFinish(FibreField,Err)

 ! create the gauss point based field, for testing
  IF(TEST_GAUSS_POINT_FIELD) THEN
  WRITE(*,*) '---------<TESTING GAUSS POINT FIELD>---------'
  CALL CMISSFieldTypeInitialise(GPfield,Err)
  CALL CMISSFieldCreateStart(FieldGPUserNumber,Region,GPfield,Err)
  CALL CMISSFieldTypeSet(GPfield,CMISSFieldGeneralType,Err) ! ?
  CALL CMISSFieldMeshDecompositionSet(GPfield,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(GPfield,GeometricField,Err)

  CALL CMISSFieldNumberOfComponentsSet(GPfield,CMISSFieldUVariableType,2,Err)



  ! FOR COMPONENTS
  CALL CMISSFieldComponentInterpolationSet(GPfield,CMISSFieldUVariableType,1,CMISSFieldGaussPointBasedInterpolation,Err) ! GP based
  CALL CMISSFieldComponentInterpolationSet(GPfield,CMISSFieldUVariableType,2,CMISSFieldGaussPointBasedInterpolation,Err) ! GP based

  CALL CMISSFieldCreateFinish(GPfield,Err)

  CALL CMISSFieldComponentValuesInitialise(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,3.14_CMISSDP,Err) ! init!
  CALL CMISSFieldComponentValuesInitialise(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,2.17_CMISSDP,Err) ! init!
  CALL CMISSFieldComponentValuesInitialise(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,4.14_CMISSDP,Err) ! set to const

  ! test gauss point field
  D=0;
  DO E=1,NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements
  DO I=1,8
    CALL CMISSFieldParameterSetGetGaussPoint(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,E,I,1,TMP,Err)
    CALL CMISSFieldParameterSetUpdateGaussPoint(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,E,I,1,TMP+D,Err)
    CALL CMISSFieldParameterSetGetGaussPoint(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,E,I,2,TMP,Err)
    CALL CMISSFieldParameterSetUpdateGaussPoint(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,E,I,2,TMP+D,Err)
    D=D+1
  ENDDO
  ENDDO

  D=0;
  DO E=1,NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements
  DO I=1,8
    CALL CMISSFieldParameterSetGetGaussPoint(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,E,I,1,TMP,Err)
    WRITE(*,*) 'COMPONENT 1 ELEMENT ',E,', GP ', I, ' = ',TMP, ' Should = ', 3.14 + 1 + D
    CALL CMISSFieldParameterSetGetGaussPoint(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,E,I,2,TMP,Err)
    WRITE(*,*) 'COMPONENT 2 ELEMENT ',E,', GP ', I, ' = ',TMP, ' Should = ', 2.17 + D
    D=D+1;
  ENDDO
  ENDDO
  WRITE(*,*) '---------</TESTING GAUSS POINT FIELD>---------'
  ENDIF  ! TEST GP FIELD

  !Create the equations_set
  CALL CMISSEquationsSetCreateStart(EquationSetUserNumber,Region,FibreField,EquationsSet,Err)
  CALL CMISSEquationsSetSpecificationSet(EquationsSet,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetActiveContractionSubtype,Err) ! CHANGED
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the dependent field
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL CMISSFieldTypeInitialise(MaterialField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  !Set Costa constants ! TODO: CORRECT
  DO I=1,7
   CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,I,1.0_CMISSDP,Err)
  END DO
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldVVariableType,CMISSFieldValuesSetType,1,2.0_CMISSDP,Err) ! activate at time 2. TODO: inhomogeneous


  ! create independent field
  CALL CMISSEquationsSetIndependentCreateStart(EquationsSet,IndependentFieldUserNumber,IndependentField,Err)
  CALL CMISSEquationsSetIndependentCreateFinish(EquationsSet,Err)

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
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,-8.0_CMISSDP,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSet,BoundaryConditions,Err)


  DO z_idx=1,NumberGlobalZElements+1
    DO y_idx=1,NumberGlobalYElements+1
      !x=0 nodes
      NodeNumber1=1+(y_idx-1)*(NumberGlobalXElements+1)+(z_idx-1)*(NumberGlobalXElements+1)*(NumberGlobalYElements+1)
      CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber1,1,Node1Domain,Err)
      IF(Node1Domain==ComputationalNodeNumber) THEN
       DO D=1,1 ! just fix left face x position
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,NodeNumber1,D, &
          & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
       ENDDO
      ENDIF
    ENDDO !y_idx
  ENDDO !z_idx


  CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSet,Err)

  !Define the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemElasticityClass,CMISSProblemFiniteElasticityType, &
    & CMISSProblemQuasistaticFiniteElasticitySubtype,Err) ! CHANGED TO CMISSProblemQuasistaticFiniteElasticitySubtype
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
   CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)

   CALL CMISSControlLoopTimesSet(ControlLoop, START_TIME - DT, END_TIME, DT, Err) ! set begin/end timings  . START AT -DT TO SOLVE FOR 0 AS WELL
   CALL CMISSControlLoopTimeOutputSet(ControlLoop,1,Err)    !Set the output timing
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(Solver,CMISSSolverNewtonJacobianFDCalculated,Err)
  CALL CMISSSolverNewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearDirectSolveType,Err)
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)   
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Solve problem
  CALL CMISSProblemSolve(Problem,Err)

  !Output solution  
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"ActiveContraction","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"ActiveContraction","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM ActiveContractionExample

