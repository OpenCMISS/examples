!> \file
!> \author Sander Land
!> \brief This is an example program to solve active contraction based finite elasticity equation using OpenCMISS calls.
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

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP

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
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=7

  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2


  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  REAL(CMISSRP), PARAMETER :: START_TIME = 0.0, END_TIME = 1000.0, DT = 1  ! ms

  LOGICAL, PARAMETER :: TEST_GAUSS_POINT_FIELD = .FALSE.

  !Program types

  !Program variables

  INTEGER(CMISSIntg), PARAMETER, DIMENSION(1:27) :: ROTATE_ELEM = [ 1, 2, 3,10,11,12,19,20,21, &
                                                                  &  4, 5, 6,13,14,15,22,23,24, &  ! swap xi2 and xi3 directions for fiber angle
                                                                  &  7, 8, 9,16,17,18,25,26,27 ]  ! swap xi2 and xi3 directions for fiber angle

  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: D, E, N, I

  REAL(CMISSRP) :: TMP
  REAL(CMISSRP), DIMENSION(1:7) :: COSTA_PARAMS =  [ 0.2, 30.0, 12.0, 14.0, 14.0, 10.0, 18.0 ] ! a bff bfs bfn bss bsn bnn

  INTEGER(CMISSIntg), dimension(:,:), allocatable :: Elements
  REAL(CMISSRP)     , dimension(:,:), allocatable :: Nodes
  REAL(CMISSRP)     , dimension(:,:), allocatable :: DirichletConditions
  REAL(CMISSRP)     , dimension(:,:), allocatable :: Fibers
  REAL(CMISSRP)     , dimension(:,:), allocatable :: ActivationTimes

  !CMISS variables
  TYPE(cmfe_BasisType) :: QuadraticBasis, LinearBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,FibreField,MaterialField,DependentField, GPfield, IndependentField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_MeshElementsType) :: QuadraticElements,LinearElements
  TYPE(cmfe_NodesType) :: CMNodes
  TYPE(cmfe_ControlLoopType) :: ControlLoop

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: Err

  !Intialise cmiss
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  !Set all diganostic levels on for testing
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_RESIDUAL_EVALUATE"],Err)

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)


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
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange
  CALL cmfe_Basis_Initialise(LinearBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)

  CALL cmfe_Basis_Initialise(QuadraticBasis,Err)
  CALL cmfe_Basis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL cmfe_Basis_InterpolationXiSet(QuadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_CreateFinish(QuadraticBasis,Err)

  !Create a mesh with two components, Quadratic for geometry and fibers and linear lagrange
  !for hydrostatic pressure and material properties
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,3,Mesh,Err) ! dim = 3
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,size(Elements,2),Err) ! num elts
  
  !define nodes for the mesh
  CALL cmfe_Nodes_Initialise(CMNodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,size(Nodes,2),CMNodes,Err) ! num nodes
  CALL cmfe_Nodes_CreateFinish(CMNodes,Err)
  !Quadratic component : from file
  CALL cmfe_MeshElements_Initialise(QuadraticElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,QuadraticMeshComponentNumber,QuadraticBasis,QuadraticElements,Err)
  DO E=1,size(Elements,2)
    CALL cmfe_MeshElements_NodesSet(QuadraticElements,E, Elements(ROTATE_ELEM,E),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(QuadraticElements,Err)
  !linear Lagrange component: numbers do not need to be continuous from 1? -> use quadratic nodeno on corners
  CALL cmfe_MeshElements_Initialise(LinearElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,LinearMeshComponentNumber,LinearBasis,LinearElements,Err)
  DO E=1,size(Elements,2)
    CALL cmfe_MeshElements_NodesSet(LinearElements,E, Elements( [1,3,7,9,19,21,25,27],E),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(LinearElements,Err)

  !finish mesh creation
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)



  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)  
  CALL cmfe_Field_NumberOfVariablesSet(GeometricField,1,Err) ! 1 var
  CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,Err)   ! 3 components of geom field
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Set node positions
  DO N=1,size(Nodes,2)
  DO D=1,3
    CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,N,D,Nodes(D, &
      & N),Err)
  ENDDO
  ENDDO


  !Create a fibre field and attach it to the geometric field  
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,Err)   ! 1 var, 3 components -> angles!
  DO D=1,3
    CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,D,QuadraticMeshComponentNumber,Err) ! quadratic interp
    CALL cmfe_Field_ComponentInterpolationSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,D,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err) ! node based
  ENDDO
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !Set fiber directions
  DO N=1,size(Nodes,2)
  DO D=1,3
    Fibers(D,N) = 0 ! input file gives unit vectors and opencmiss expects angles. TODO: fix.
    CALL cmfe_Field_ParameterSetUpdateNode(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,N,D,Fibers(D,N), &
      & Err)
  ENDDO
  ENDDO



 ! create the gauss point based field, for testing
  IF(TEST_GAUSS_POINT_FIELD) THEN
  WRITE(*,*) '---------<TESTING GAUSS POINT FIELD>---------'
  CALL cmfe_Field_Initialise(GPfield,Err)
  CALL cmfe_Field_CreateStart(FieldGPUserNumber,Region,GPfield,Err)
  CALL cmfe_Field_TypeSet(GPfield,CMFE_FIELD_GENERAL_TYPE,Err) ! ?
  CALL cmfe_Field_MeshDecompositionSet(GPfield,Decomposition,Err)        
  CALL cmfe_Field_GeometricFieldSet(GPfield,GeometricField,Err)

  CALL cmfe_Field_NumberOfComponentsSet(GPfield,CMFE_FIELD_U_VARIABLE_TYPE,2,Err)

  ! FOR COMPONENTS
  CALL cmfe_Field_ComponentInterpolationSet(GPfield,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! GP based
  CALL cmfe_Field_ComponentInterpolationSet(GPfield,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! GP based

  CALL cmfe_Field_CreateFinish(GPfield,Err)

  CALL cmfe_Field_ComponentValuesInitialise(GPfield,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,3.14_CMISSRP,Err) ! init!
  CALL cmfe_Field_ComponentValuesInitialise(GPfield,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,2.17_CMISSRP,Err) ! init!
  CALL cmfe_Field_ComponentValuesInitialise(GPfield,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,4.14_CMISSRP,Err) ! set to const

  ! test gauss point field
  D=0;
  DO E=1,size(Elements,2)
  DO I=1,8
    CALL cmfe_Field_ParameterSetGetGaussPoint(GPfield,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,I,E,1,TMP,Err)
    CALL cmfe_Field_ParameterSetUpdateGaussPoint(GPfield,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,I,E,1,TMP+D,Err)
    CALL cmfe_Field_ParameterSetGetGaussPoint(GPfield,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,I,E,2,TMP,Err)
    CALL cmfe_Field_ParameterSetUpdateGaussPoint(GPfield,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,I,E,2,TMP+D,Err)
    D=D+1
  ENDDO
  ENDDO

  D=0;
  DO E=1,size(Elements,2)
  DO I=1,8
    CALL cmfe_Field_ParameterSetGetGaussPoint(GPfield,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,I,E,1,TMP,Err)
    WRITE(*,*) 'COMPONENT 1 ELEMENT ',E,', GP ', I, ' = ',TMP, ' Should = ', 3.14 + 1 + D
    CALL cmfe_Field_ParameterSetGetGaussPoint(GPfield,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,I,E,2,TMP,Err)
    WRITE(*,*) 'COMPONENT 2 ELEMENT ',E,', GP ', I, ' = ',TMP, ' Should = ', 2.17 + D
    D=D+1;
  ENDDO
  ENDDO
  WRITE(*,*) '---------</TESTING GAUSS POINT FIELD>---------'
  ENDIF  ! TEST GP FIELD

  !Create the equations_set
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  ! CHANGED
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set Costa material parameters
  DO I=1,7
   CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,I, &
     & COSTA_PARAMS(I),Err)
  END DO

!  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2.0_CMISSRP,Err) ! activate at time 2. TODO: inhomogeneous
  ! inhomogeneous activation times from file
  DO E=1,size(Elements,2)
  DO I=1,27
    CALL cmfe_Field_ParameterSetUpdateGaussPoint(MaterialField,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & I, E, 1, ActivationTimes(E,ROTATE_ELEM(I)), Err) ! rotating an element with 27 nodes can be done in the same way as the Gauss points
  ENDDO
  ENDDO

  ! create independent field
  CALL cmfe_Field_Initialise(IndependentField,Err)
  CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSet,IndependentFieldUserNumber,IndependentField,Err)
  CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSet,Err)


  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)   

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,-0.0_CMISSRP, &
    & Err) ! -8?

  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE],Problem,Err)
   ! CHANGED TO CMFE_PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
   CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)

   CALL cmfe_ControlLoop_TimesSet(ControlLoop, START_TIME - DT, END_TIME, DT, Err) ! set begin/end timings  . START AT -DT TO SOLVE FOR 0 AS WELL
   CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,1,Err)    !Set the output timing
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(Solver,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err) ! faster than CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED ?
  CALL cmfe_Solver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)   
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)


  DO I=1,size(DirichletConditions,2)
    N = INT(DirichletConditions(1,I))
    D = INT(DirichletConditions(2,I))
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,N,D,&
         & CMFE_BOUNDARY_CONDITION_FIXED, Nodes(D,N) + DirichletConditions(3,I),Err)  ! current + offset
  ENDDO


  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !Output solution  
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"ActiveContraction","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"ActiveContraction","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP


contains
  subroutine read_mesh(fp, elements, node_coords, fixed_nodes, fibers)
    INTEGER(CMISSIntg), intent(in)    :: fp  !< file 'pointer'
    INTEGER(CMISSIntg), dimension(:,:), allocatable, intent(inout) :: elements !< element topology
    REAL(CMISSRP)     , dimension(:,:), allocatable, intent(inout) :: node_coords  !< initial positions etc
    REAL(CMISSRP)     , dimension(:,:), allocatable, intent(inout) :: fixed_nodes  !< dirichlet boundary conditions
    REAL(CMISSRP)     , dimension(:,:), allocatable, intent(inout) :: fibers       !< unit vectors for fiber dir

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
    REAL(CMISSRP), dimension(:,:), allocatable, intent(inout) :: activtime !< elements x 27 array of activation times
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

