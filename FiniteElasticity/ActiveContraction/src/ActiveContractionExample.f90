!> \file
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

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField


  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=1.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldGPUserNumber=5 ! temp/test
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumber=6

  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2


  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  REAL(CMISSDP), PARAMETER :: START_TIME = 0.0, END_TIME = 1000.0, DT = 1  ! ms

  LOGICAL, PARAMETER :: TEST_GAUSS_POINT_FIELD = .FALSE.

  !Program types

  !Program variables

  INTEGER(CMISSIntg), PARAMETER, DIMENSION(1:27) :: ROTATE_ELEM = (/ 1, 2, 3,10,11,12,19,20,21, &
                                                                  &  4, 5, 6,13,14,15,22,23,24, &  ! swap xi2 and xi3 directions for fiber angle
                                                                  &  7, 8, 9,16,17,18,25,26,27 /)  ! swap xi2 and xi3 directions for fiber angle

  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: D, E, N, I

  REAL(CMISSDP) :: TMP
  REAL(CMISSDP), DIMENSION(1:7) :: COSTA_PARAMS =  (/ 0.2, 30.0, 12.0, 14.0, 14.0, 10.0, 18.0 /) ! a bff bfs bfn bss bsn bnn

  INTEGER(CMISSIntg), dimension(:,:), allocatable :: Elements
  REAL(CMISSDP)     , dimension(:,:), allocatable :: Nodes
  REAL(CMISSDP)     , dimension(:,:), allocatable :: DirichletConditions
  REAL(CMISSDP)     , dimension(:,:), allocatable :: Fibers
  REAL(CMISSDP)     , dimension(:,:), allocatable :: ActivationTimes

  !CMISS variables
  TYPE(CMISSBasisType) :: QuadraticBasis, LinearBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,FibreField,MaterialField,DependentField, GPfield, IndependentField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSMeshElementsType) :: QuadraticElements,LinearElements
  TYPE(CMISSNodesType) :: CMNodes
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


!  open(unit = 2, file = "./data/hollowcylq-221.in")
!  open(unit = 3, file = "./data/hollowcylq-221.in.gpactiv")
  open(unit = 2, file = "./data/lvq-842.in")
  open(unit = 3, file = "./data/lvq-842.in.gpactiv")
  call read_mesh(2, Elements, Nodes, DirichletConditions, Fibers)
  call read_activation_times(3, ActivationTimes)
  close(2)
  close(3)

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

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange
  CALL CMISSBasisTypeInitialise(LinearBasis,Err)
  CALL CMISSBasisCreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasisCreateFinish(LinearBasis,Err)

  CALL CMISSBasisTypeInitialise(QuadraticBasis,Err)
  CALL CMISSBasisCreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasisInterpolationXiSet(QuadraticBasis,(/CMISSBasisQuadraticLagrangeInterpolation, &
    & CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
  CALL CMISSBasisCreateFinish(QuadraticBasis,Err)

  !Create a mesh with two components, Quadratic for geometry and fibers and linear lagrange
  !for hydrostatic pressure and material properties
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,3,Mesh,Err) ! dim = 3
  CALL CMISSMeshNumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  CALL CMISSMeshNumberOfElementsSet(Mesh,size(Elements,2),Err) ! num elts
  
  !define nodes for the mesh
  CALL CMISSNodesTypeInitialise(CMNodes,Err)
  CALL CMISSNodesCreateStart(Region,size(Nodes,2),CMNodes,Err) ! num nodes
  CALL CMISSNodesCreateFinish(CMNodes,Err)
  !Quadratic component : from file
  CALL CMISSMeshElementsTypeInitialise(QuadraticElements,Err)
  CALL CMISSMeshElementsCreateStart(Mesh,QuadraticMeshComponentNumber,QuadraticBasis,QuadraticElements,Err)
  DO E=1,size(Elements,2)
    CALL CMISSMeshElementsNodesSet(QuadraticElements,E, Elements(ROTATE_ELEM,E),Err)
  ENDDO
  CALL CMISSMeshElementsCreateFinish(QuadraticElements,Err)
  !linear Lagrange component: numbers do not need to be continuous from 1? -> use quadratic nodeno on corners
  CALL CMISSMeshElementsTypeInitialise(LinearElements,Err)
  CALL CMISSMeshElementsCreateStart(Mesh,LinearMeshComponentNumber,LinearBasis,LinearElements,Err)
  DO E=1,size(Elements,2)
    CALL CMISSMeshElementsNodesSet(LinearElements,E, Elements( (/1,3,7,9,19,21,25,27/),E),Err)
  ENDDO
  CALL CMISSMeshElementsCreateFinish(LinearElements,Err)

  !finish mesh creation
  CALL CMISSMeshCreateFinish(Mesh,Err)



  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry)
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)  
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,1,Err) ! 1 var
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,3,Err)   ! 3 components of geom field
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Set node positions
  DO N=1,size(Nodes,2)
  DO D=1,3
    CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,N,D,Nodes(D,N),Err)
  ENDDO
  ENDDO


  !Create a fibre field and attach it to the geometric field  
  CALL CMISSFieldTypeInitialise(FibreField,Err)
  CALL CMISSFieldCreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSFieldTypeSet(FibreField,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreField,CMISSFieldUVariableType,3,Err)   ! 1 var, 3 components -> angles!
  DO D=1,3
    CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,D,QuadraticMeshComponentNumber,Err) ! quadratic interp
    CALL CMISSFieldComponentInterpolationSet(FibreField,CMISSFieldUVariableType,D,CMISSFieldNodeBasedInterpolation,Err) ! node based
  ENDDO
  CALL CMISSFieldCreateFinish(FibreField,Err)

  !Set fiber directions
  DO N=1,size(Nodes,2)
  DO D=1,3
    Fibers(D,N) = 0 ! input file gives unit vectors and opencmiss expects angles. TODO: fix.
    CALL CMISSFieldParameterSetUpdateNode(FibreField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,N,D,Fibers(D,N),Err)
  ENDDO
  ENDDO



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
  DO E=1,size(Elements,2)
  DO I=1,8
    CALL CMISSFieldParameterSetGetGaussPoint(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,E,I,1,TMP,Err)
    CALL CMISSFieldParameterSetUpdateGaussPoint(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,E,I,1,TMP+D,Err)
    CALL CMISSFieldParameterSetGetGaussPoint(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,E,I,2,TMP,Err)
    CALL CMISSFieldParameterSetUpdateGaussPoint(GPfield,CMISSFieldUVariableType,CMISSFieldValuesSetType,E,I,2,TMP+D,Err)
    D=D+1
  ENDDO
  ENDDO

  D=0;
  DO E=1,size(Elements,2)
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
    CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
CALL CMISSEquationsSetCreateStart(EquationSetUserNumber,Region,FibreField,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetActiveContractionSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
   ! CHANGED
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the dependent field
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL CMISSFieldTypeInitialise(MaterialField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  !Set Costa material parameters
  DO I=1,7
   CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,I,COSTA_PARAMS(I),Err)
  END DO

!  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldVVariableType,CMISSFieldValuesSetType,1,2.0_CMISSDP,Err) ! activate at time 2. TODO: inhomogeneous
  ! inhomogeneous activation times from file
  DO E=1,size(Elements,2)
  DO I=1,27
    CALL CMISSFieldParameterSetUpdateGaussPoint(MaterialField,CMISSFieldVVariableType,CMISSFieldValuesSetType, &
    & E, I, 1, ActivationTimes(E,ROTATE_ELEM(I)), Err) ! rotating an element with 27 nodes can be done in the same way as the Gauss points
  ENDDO
  ENDDO

  ! create independent field
  CALL CMISSFieldTypeInitialise(IndependentField,Err)
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
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,-0.0_CMISSDP,Err) ! -8?

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
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(Solver,CMISSSolverNewtonJacobianAnalyticCalculated,Err) ! faster than CMISSSolverNewtonJacobianFDCalculated ?
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

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)


  DO I=1,size(DirichletConditions,2)
    N = INT(DirichletConditions(1,I))
    D = INT(DirichletConditions(2,I))
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,N,D,&
         & CMISSBoundaryConditionFixed, Nodes(D,N) + DirichletConditions(3,I),Err)  ! current + offset
  ENDDO


  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

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


contains
  subroutine read_mesh(fp, elements, node_coords, fixed_nodes, fibers)
    INTEGER(CMISSIntg), intent(in)    :: fp  !< file 'pointer'
    INTEGER(CMISSIntg), dimension(:,:), allocatable, intent(inout) :: elements !< element topology
    REAL(CMISSDP)     , dimension(:,:), allocatable, intent(inout) :: node_coords  !< initial positions etc
    REAL(CMISSDP)     , dimension(:,:), allocatable, intent(inout) :: fixed_nodes  !< dirichlet boundary conditions
    REAL(CMISSDP)     , dimension(:,:), allocatable, intent(inout) :: fibers       !< unit vectors for fiber dir

    INTEGER(CMISSIntg) :: number_of_elts, number_of_dims, number_of_fixednodes, maxnodenr, i,j, eltno, nodeno

    character(len=256) :: desc_str
 
    maxnodenr = 0

    read (fp,*) desc_str
    write (*,*) "'number_of_elements':", desc_str

    read (fp,*) number_of_elts
    write (*,*) "# of elements:", number_of_elts

    allocate(elements(1:27,1:number_of_elts))

    read (fp,*) desc_str
    write (*,*) "'topology':", desc_str
 
    do i = 1, number_of_elts
      read (fp,*)  eltno, elements(:,i)
      write (*,*) "nodes of element ", eltno, " are : ", elements(:,i)
      do j=1,27
        maxnodenr = max(maxnodenr,elements(j,i))
      end do
    end do
   

    read (fp,*) desc_str
    read (fp,*) number_of_dims
    write (*,*) "max node nr:", maxnodenr, " dimensions: ", number_of_dims
    write (*,*) "'initial_positions':", desc_str
 
    allocate(node_coords(1:number_of_dims,1:maxnodenr))
    do i = 1, maxnodenr
      read (fp,*) nodeno, node_coords(:,i)  ! expect input in order
      write (*,*) "initial value of node ", nodeno, "=", i ," is ", node_coords(:,i) 
    end do

    read (fp,*) desc_str
    write (*,*) "'fixed_nodes':", desc_str
    read (fp,*) number_of_fixednodes
    write (*,*) "# of dirichlet bc:", number_of_fixednodes
    allocate(fixed_nodes(3,number_of_fixednodes))
    do i = 1, number_of_fixednodes
      read (fp,*) fixed_nodes(:,i)
      write (*,*) "dirichlet bc: ", fixed_nodes(:,i) 
    end do
  
    read (fp,*) desc_str ! traction
    read (fp,*) desc_str ! 0 0
    read (fp,*) desc_str
    write (*,*) "'fibers':", desc_str

    allocate(fibers(3,1:maxnodenr))
    do i = 1, maxnodenr
      read (fp,*) nodeno, fibers(:,i)  ! expect input in order
      write (*,*) "fiber dir at ", nodeno, "=", i ," is ", fibers(:,i) 
    end do


  end subroutine read_mesh    


  subroutine read_activation_times(fp, activtime)
    INTEGER(CMISSIntg), intent(in)    :: fp  !< file 'pointer'
    REAL(CMISSDP), dimension(:,:), allocatable, intent(inout) :: activtime !< elements x 27 array of activation times
    INTEGER(CMISSIntg)   :: num_elt, num_gp, I
    read (fp,*) num_elt, num_gp
    write(*,*) 'reading activation times on ',num_elt,' elements / ',num_gp,' Gauss points'
    allocate( activtime(1:num_elt,1:num_gp))
    do I=1,num_elt
      read (fp,*) activtime(I,:)
      write(*,*) 'element ',I,' : ',activtime(I,:)
    enddo
  end subroutine read_activation_times

END PROGRAM ActiveContractionExample

