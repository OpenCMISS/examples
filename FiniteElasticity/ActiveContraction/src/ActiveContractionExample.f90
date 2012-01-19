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

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  !Set all diganostic levels on for testing
  !CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_RESIDUAL_EVALUATE"/),Err)

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
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange
  CALL CMISSBasis_Initialise(LinearBasis,Err)
  CALL CMISSBasis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasis_CreateFinish(LinearBasis,Err)

  CALL CMISSBasis_Initialise(QuadraticBasis,Err)
  CALL CMISSBasis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasis_InterpolationXiSet(QuadraticBasis,(/CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & (/CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME/),Err)
  CALL CMISSBasis_CreateFinish(QuadraticBasis,Err)

  !Create a mesh with two components, Quadratic for geometry and fibers and linear lagrange
  !for hydrostatic pressure and material properties
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,3,Mesh,Err) ! dim = 3
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  CALL CMISSMesh_NumberOfElementsSet(Mesh,size(Elements,2),Err) ! num elts
  
  !define nodes for the mesh
  CALL CMISSNodes_Initialise(CMNodes,Err)
  CALL CMISSNodes_CreateStart(Region,size(Nodes,2),CMNodes,Err) ! num nodes
  CALL CMISSNodes_CreateFinish(CMNodes,Err)
  !Quadratic component : from file
  CALL CMISSMeshElements_Initialise(QuadraticElements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,QuadraticMeshComponentNumber,QuadraticBasis,QuadraticElements,Err)
  DO E=1,size(Elements,2)
    CALL CMISSMeshElements_NodesSet(QuadraticElements,E, Elements(ROTATE_ELEM,E),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(QuadraticElements,Err)
  !linear Lagrange component: numbers do not need to be continuous from 1? -> use quadratic nodeno on corners
  CALL CMISSMeshElements_Initialise(LinearElements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,LinearMeshComponentNumber,LinearBasis,LinearElements,Err)
  DO E=1,size(Elements,2)
    CALL CMISSMeshElements_NodesSet(LinearElements,E, Elements( (/1,3,7,9,19,21,25,27/),E),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(LinearElements,Err)

  !finish mesh creation
  CALL CMISSMesh_CreateFinish(Mesh,Err)



  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry)
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)  
  CALL CMISSField_NumberOfVariablesSet(GeometricField,1,Err) ! 1 var
  CALL CMISSField_NumberOfComponentsSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,Err)   ! 3 components of geom field
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Set node positions
  DO N=1,size(Nodes,2)
  DO D=1,3
    CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,N,D,Nodes(D, &
      & N),Err)
  ENDDO
  ENDDO


  !Create a fibre field and attach it to the geometric field  
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_NumberOfComponentsSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,3,Err)   ! 1 var, 3 components -> angles!
  DO D=1,3
    CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,D,QuadraticMeshComponentNumber,Err) ! quadratic interp
    CALL CMISSField_ComponentInterpolationSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,D,CMISS_FIELD_NODE_BASED_INTERPOLATION,Err) ! node based
  ENDDO
  CALL CMISSField_CreateFinish(FibreField,Err)

  !Set fiber directions
  DO N=1,size(Nodes,2)
  DO D=1,3
    Fibers(D,N) = 0 ! input file gives unit vectors and opencmiss expects angles. TODO: fix.
    CALL CMISSField_ParameterSetUpdateNode(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,N,D,Fibers(D,N), &
      & Err)
  ENDDO
  ENDDO



 ! create the gauss point based field, for testing
  IF(TEST_GAUSS_POINT_FIELD) THEN
  WRITE(*,*) '---------<TESTING GAUSS POINT FIELD>---------'
  CALL CMISSField_Initialise(GPfield,Err)
  CALL CMISSField_CreateStart(FieldGPUserNumber,Region,GPfield,Err)
  CALL CMISSField_TypeSet(GPfield,CMISS_FIELD_GENERAL_TYPE,Err) ! ?
  CALL CMISSField_MeshDecompositionSet(GPfield,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(GPfield,GeometricField,Err)

  CALL CMISSField_NumberOfComponentsSet(GPfield,CMISS_FIELD_U_VARIABLE_TYPE,2,Err)

  ! FOR COMPONENTS
  CALL CMISSField_ComponentInterpolationSet(GPfield,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! GP based
  CALL CMISSField_ComponentInterpolationSet(GPfield,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! GP based

  CALL CMISSField_CreateFinish(GPfield,Err)

  CALL CMISSField_ComponentValuesInitialise(GPfield,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,3.14_CMISSDP,Err) ! init!
  CALL CMISSField_ComponentValuesInitialise(GPfield,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,2.17_CMISSDP,Err) ! init!
  CALL CMISSField_ComponentValuesInitialise(GPfield,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,4.14_CMISSDP,Err) ! set to const

  ! test gauss point field
  D=0;
  DO E=1,size(Elements,2)
  DO I=1,8
    CALL CMISSField_ParameterSetGetGaussPoint(GPfield,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,E,I,1,TMP,Err)
    CALL CMISSField_ParameterSetUpdateGaussPoint(GPfield,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,E,I,1,TMP+D,Err)
    CALL CMISSField_ParameterSetGetGaussPoint(GPfield,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,E,I,2,TMP,Err)
    CALL CMISSField_ParameterSetUpdateGaussPoint(GPfield,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,E,I,2,TMP+D,Err)
    D=D+1
  ENDDO
  ENDDO

  D=0;
  DO E=1,size(Elements,2)
  DO I=1,8
    CALL CMISSField_ParameterSetGetGaussPoint(GPfield,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,E,I,1,TMP,Err)
    WRITE(*,*) 'COMPONENT 1 ELEMENT ',E,', GP ', I, ' = ',TMP, ' Should = ', 3.14 + 1 + D
    CALL CMISSField_ParameterSetGetGaussPoint(GPfield,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,E,I,2,TMP,Err)
    WRITE(*,*) 'COMPONENT 2 ELEMENT ',E,', GP ', I, ' = ',TMP, ' Should = ', 2.17 + D
    D=D+1;
  ENDDO
  ENDDO
  WRITE(*,*) '---------</TESTING GAUSS POINT FIELD>---------'
  ENDIF  ! TEST GP FIELD

  !Create the equations_set
    CALL CMISSField_Initialise(EquationsSetField,Err)
CALL CMISSEquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
   ! CHANGED
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL CMISSField_Initialise(MaterialField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set Costa material parameters
  DO I=1,7
   CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,I, &
     & COSTA_PARAMS(I),Err)
  END DO

!  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2.0_CMISSDP,Err) ! activate at time 2. TODO: inhomogeneous
  ! inhomogeneous activation times from file
  DO E=1,size(Elements,2)
  DO I=1,27
    CALL CMISSField_ParameterSetUpdateGaussPoint(MaterialField,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & E, I, 1, ActivationTimes(E,ROTATE_ELEM(I)), Err) ! rotating an element with 27 nodes can be done in the same way as the Gauss points
  ENDDO
  ENDDO

  ! create independent field
  CALL CMISSField_Initialise(IndependentField,Err)
  CALL CMISSEquationsSet_IndependentCreateStart(EquationsSet,IndependentFieldUserNumber,IndependentField,Err)
  CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSet,Err)


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
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,-0.0_CMISSDP, &
    & Err) ! -8?

  !Define the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_ELASTICITY_CLASS,CMISS_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMISS_PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE,Err) ! CHANGED TO CMISS_PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
   CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)

   CALL CMISSControlLoop_TimesSet(ControlLoop, START_TIME - DT, END_TIME, DT, Err) ! set begin/end timings  . START AT -DT TO SOLVE FOR 0 AS WELL
   CALL CMISSControlLoop_TimeOutputSet(ControlLoop,1,Err)    !Set the output timing
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err) ! faster than CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED ?
  CALL CMISSSolver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)   
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)


  DO I=1,size(DirichletConditions,2)
    N = INT(DirichletConditions(1,I))
    D = INT(DirichletConditions(2,I))
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,N,D,&
         & CMISS_BOUNDARY_CONDITION_FIXED, Nodes(D,N) + DirichletConditions(3,I),Err)  ! current + offset
  ENDDO


  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Output solution  
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"ActiveContraction","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"ActiveContraction","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

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

