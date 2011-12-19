!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve an Analytic Laplace equation using OpenCMISS calls.
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

!> \example ClassicalField/Laplace/ANALYTIC_LINEAR_ELASTICITY/src/ANALYTIC_LINEAR_ELASTICITYExample.f90
!! Example illustrating the use of OpenCMISS to solve the Laplace problem and check with its Analytic Solution.
!! 
!! \htmlinclude ClassicalField/Laplace/ANALYTIC_LINEAR_ELASTICITY/history.html
!< 

!> Main program
PROGRAM ANALYTIC_LINEAR_ELASTICITYEXAMPLE

  USE MPI
  USE OPENCMISS
  USE TEST_FRAMEWORK_ROUTINES

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField


  !Test program parameters

  REAL(CMISSDP), PARAMETER :: ORIGIN(3)=(/0.0_CMISSDP,-1.0_CMISSDP,-1.0_CMISSDP/)
  REAL(CMISSDP), PARAMETER :: LENGTH=5.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: HEIGHT=2.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: NumberOfDomains=1

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber = 1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariables=2

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariables=1

  INTEGER(CMISSIntg), PARAMETER :: FieldAnalyticUserNumber=4

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  REAL(CMISSDP), PARAMETER ::   ZERO = 0.0_CMISSDP

  !Program types

  TYPE(CMISSRegionType) :: WorldRegion
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem

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

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_FINITE_ELEMENT_CALCULATE"/),Err)

  CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT(2,2,2,"TriLinearLagrange")

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

CONTAINS


  !
  !================================================================================================================================
  !  
    !>Check if the convergence of linear langrange interpolation is expected.
  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT(NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements,OutputFilename)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NumberGlobalXElements !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NumberGlobalYElements !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NumberGlobalZElements !<increment interval number of elements per axis
    CHARACTER(LEN=*), INTENT(IN) :: OutputFilename !<The Error condition string
    !Local Variables
    TYPE(CMISSFieldType) :: DependentField

    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC(NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements, &
      & CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION,DependentField)

    CALL CMISSAnalyticAnalysisOutput(DependentField,OutputFilename,Err)
    
    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber, &
      & GeneratedMeshUserNumber,ProblemUserNumber)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT
  
  !
  !================================================================================================================================
  !   
    
  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC(NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements, &
    & InterpolationSpecifications,DependentField)
    !Argument variables 
    INTEGER(CMISSIntg), INTENT(IN) :: NumberGlobalXElements !<number of elements on x axis
    INTEGER(CMISSIntg), INTENT(IN) :: NumberGlobalYElements !<number of elements on y axis
    INTEGER(CMISSIntg), INTENT(IN) :: NumberGlobalZElements !<number of elements on z axis
    INTEGER(CMISSIntg), INTENT(IN) :: InterpolationSpecifications !<the interpolation specifications
    TYPE(CMISSFieldType) :: DependentField

    !Program variables
    REAL(CMISSDP) :: MeshDimensions(3),MaterialParameters(6)
    INTEGER(CMISSIntg) :: AnalyticFunction,Interpolation(3),NumberOfGaussPoints(3),EquationSetSubtype
    INTEGER(CMISSIntg) :: FieldGeometryNumberOfComponents,FieldDependentNumberOfComponents,NumberOfElements(3)
    INTEGER(CMISSIntg) :: MPI_IERROR
    INTEGER(CMISSIntg) :: EquationsSetIndex,FieldComponentIndex,FieldMaterialNumberOfComponents,NumberOfXi
    INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber

    !CMISS variables

    TYPE(CMISSBasisType) :: Basis
    TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
    TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
    TYPE(CMISSDecompositionType) :: Decomposition
    TYPE(CMISSEquationsType) :: Equations
    TYPE(CMISSEquationsSetType) :: EquationsSet
    TYPE(CMISSFieldType) :: AnalyticField,GeometricField,MaterialField
    TYPE(CMISSMeshType) :: Mesh
    TYPE(CMISSProblemType) :: Problem
    TYPE(CMISSRegionType) :: Region
    TYPE(CMISSSolverType) :: Solver
    TYPE(CMISSSolverEquationsType) :: SolverEquations
    TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions

!    IF((NumberGlobalYElements == 0) .AND. (NumberGlobalZElements == 0)) THEN
!      NumberOfXi = 1
!      EquationSetSubtype = CMISS_EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE
!      AnalyticFunction=CMISS_EQUATIONS_SET_LINEAR_ELASTICITY_ONE_DIM_1
!      !Prescribe material properties Area,E1
!      FieldMaterialNumberOfComponents = 2 !Young's Modulus & Poisson's Ratio
!      MaterialParameters = (/WIDTH*HEIGHT,10.0E3_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
!    ELSEIF (NumberGlobalZElements == 0) THEN
!      NumberOfXi = 2
!      EquationSetSubtype = CMISS_EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE
!      AnalyticFunction=CMISS_EQUATIONS_SET_LINEAR_ELASTICITY_TWO_DIM_1
!      !Prescribe material properties h,E1,v12
!      FieldMaterialNumberOfComponents = 3 !Young's Modulus & Poisson's Ratio
!      MaterialParameters = (/HEIGHT,10.0E3_CMISSDP,0.3_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
!    ELSE
      NumberOfXi = 3
      EquationSetSubtype = CMISS_EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE
      AnalyticFunction=CMISS_EQUATIONS_SET_LINEAR_ELASTICITY_THREE_DIM_2
      !Prescribe material properties E1,E2,E3 & v13,v23,v12
      FieldMaterialNumberOfComponents = 6 !Young's Modulus & Poisson's Ratio
      MaterialParameters = (/10.0E3_CMISSDP,10.0E3_CMISSDP,10.0E3_CMISSDP,0.45_CMISSDP,0.45_CMISSDP,0.45_CMISSDP/)
!    ENDIF
    Interpolation = (/InterpolationSpecifications,InterpolationSpecifications,InterpolationSpecifications/)
    NumberOfElements = (/NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements/)
    MeshDimensions = (/LENGTH,WIDTH,HEIGHT/)
    NumberOfGaussPoints = (/4,4,4/)
    FieldGeometryNumberOfComponents=NumberOfXi
    FieldDependentNumberOfComponents=NumberOfXi

    !Get the number of computational nodes and this computational node number
    CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
    CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

    !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
    CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

    !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
    CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
    CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
    CALL CMISSCoordinateSystem_TypeSet(CoordinateSystem,CMISS_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
    CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,NumberOfXi,Err)
    CALL CMISSCoordinateSystem_OriginSet(CoordinateSystem,ORIGIN,Err)
    CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

    !Create a region and assign the CS to the region
    CALL CMISSRegion_Initialise(Region,Err)
    CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
    CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
    CALL CMISSRegion_CreateFinish(Region,Err)

    CALL CMISSBasis_Initialise(Basis,Err)
    CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
    CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CALL CMISSBasis_NumberOfXiSet(Basis,NumberOfXi,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis,Interpolation(1:NumberOfXi),Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,NumberOfGaussPoints(1:NumberOfXi),Err)
    CALL CMISSBasis_CreateFinish(Basis,Err)

    !Start the creation of a generated Mesh in the Region
    CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
    CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
    CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,1,Err)
    CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)

    !Define the Mesh on the Region
    CALL CMISSGeneratedMesh_OriginSet(GeneratedMesh,ORIGIN(1:NumberOfXi),Err)
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,MeshDimensions(1:NumberOfXi),Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,NumberOfElements(1:NumberOfXi),Err)
    CALL CMISSMesh_Initialise(Mesh,Err)
    CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,1,Mesh,Err)

    !Create a decomposition
    CALL CMISSDecomposition_Initialise(Decomposition,Err)
    CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
    CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
    CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

    !Create a field to put the geometry (defualt is geometry)
    CALL CMISSField_Initialise(GeometricField,Err)
    CALL CMISSField_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
    CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
    CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)  
    CALL CMISSField_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
    CALL CMISSField_NumberOfComponentsSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
    DO FieldComponentIndex=1,FieldGeometryNumberOfComponents
      CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,FieldComponentIndex,1,Err)
    ENDDO !FieldComponentIndex
    CALL CMISSField_CreateFinish(GeometricField,Err)

    !Update the geometric field parameters
    CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

    !Create a dependent field with two variables and three components
    CALL CMISSField_Initialise(DependentField,Err)
    CALL CMISSField_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
    CALL CMISSField_TypeSet(DependentField,CMISS_FIELD_GENERAL_TYPE,Err)  
    CALL CMISSField_MeshDecompositionSet(DependentField,Decomposition,Err)
    CALL CMISSField_GeometricFieldSet(DependentField,GeometricField,Err) 
    CALL CMISSField_DependentTypeSet(DependentField,CMISS_FIELD_DEPENDENT_TYPE,Err) 
    CALL CMISSField_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
    CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
    CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
    DO FieldComponentIndex=1,FieldDependentNumberOfComponents
      CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,FieldComponentIndex,1,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,FieldComponentIndex,1,Err)
    ENDDO !FieldComponentIndex
    CALL CMISSField_CreateFinish(DependentField,Err)

    !Create a material field and attach it to the geometric field
    CALL CMISSField_Initialise(MaterialField,Err)
    CALL CMISSField_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
    CALL CMISSField_TypeSet(MaterialField,CMISS_FIELD_MATERIAL_TYPE,Err)
    CALL CMISSField_MeshDecompositionSet(MaterialField,Decomposition,Err)
    CALL CMISSField_GeometricFieldSet(MaterialField,GeometricField,Err)
    CALL CMISSField_NumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
    CALL CMISSField_NumberOfComponentsSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)
    DO FieldComponentIndex=1,FieldMaterialNumberOfComponents
      CALL CMISSField_ComponentMeshComponentSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,FieldComponentIndex,1,Err)
    ENDDO !FieldComponentIndex
    CALL CMISSField_CreateFinish(MaterialField,Err)

    !Set isotropic elasticity material parameters - Young's Modulus & Poisson's Ratio
    DO FieldComponentIndex=1,FieldMaterialNumberOfComponents
      CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
        & FieldComponentIndex, &
        & MaterialParameters(FieldComponentIndex),Err)
    ENDDO !FieldComponentIndex

    !Create a Elasticity Class, Linear Elasticity type, no subtype, EquationsSet
    CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
      CALL CMISSField_Initialise(EquationsSetField,Err)
CALL CMISSEquationsSet_CreateStart(EquationSetUserNumber,Region,GeometricField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_LINEAR_ELASTICITY_TYPE,EquationSetSubtype,EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
    
    CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

    CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
    CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

    CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
    CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

    !Create the Equations set analtyic field variables
    CALL CMISSField_Initialise(AnalyticField,Err)
    CALL CMISSEquationsSet_AnalyticCreateStart(EquationsSet,AnalyticFunction,FieldAnalyticUserNumber,AnalyticField,Err)
    CALL CMISSEquationsSet_AnalyticCreateFinish(EquationsSet,Err)

    !Create the equations set equations
    CALL CMISSEquations_Initialise(Equations,Err)
    CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
    CALL CMISSEquations_SparsityTypeSet(EQUATIONS,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
                                                !CMISS_EQUATIONS_SPARSE_MATRICES=1 !<Use sparse matrices for the equations.
                                                !CMISS_EQUATIONS_FULL_MATRICES=2 !<Use fully populated matrices for the equations. 
    CALL CMISSEquations_OutputTypeSet(EQUATIONS,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
                                              !CMISS_EQUATIONS_NO_OUTPUT !<No output from the equations.
                                              !CMISS_EQUATIONS_TIMING_OUTPUT !<Timing information output.
                                              !CMISS_EQUATIONS_MATRIX_OUTPUT !<All below and equation matrices output.
                                              !CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT !<All below and Element matrices output.
    CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)
    
    !Define the problem
    CALL CMISSProblem_Initialise(Problem,Err)
    CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
    CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_ELASTICITY_CLASS,CMISS_PROBLEM_LINEAR_ELASTICITY_TYPE, &
      & CMISS_PROBLEM_NO_SUBTYPE,Err)
    CALL CMISSProblem_CreateFinish(Problem,Err)

    !Create the problem control loop
    CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
    CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

    !Start the creation of the Problem Solvers
    !Create the problem Solvers
    CALL CMISSSolver_Initialise(Solver,Err)
    CALL CMISSProblem_SolversCreateStart(Problem,Err)
    CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
    CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_MATRIX_OUTPUT,Err)
                                        !CMISS_SOLVER_NO_OUTPUT !<No output from the Solver routines. \see OPENCMISS_SolverOutputTypes,OPENCMISS
                                        !CMISS_SOLVER_PROGRESS_OUTPUT !<Progress output from Solver routines.
                                        !CMISS_SOLVER_TIMING_OUTPUT !<Timing output from the Solver routines plus below.
                                        !CMISS_SOLVER_SOLVER_OUTPUT !<Solver specific output from the Solver routines plus below.
                                        !CMISS_SOLVER_MATRIX_OUTPUT !<Solver matrices output from the Solver routines plus below.
    CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_PETSC_LIBRARY,Err)
                                          !CMISS_SOLVER_CMISS_LIBRARY     !<CMISS (internal) Solver library.
                                          !CMISS_SOLVER_PETSC_LIBRARY     !<PETSc Solver library.
                                          !CMISS_SOLVER_MUMPS_LIBRARY     !<MUMPS Solver library.
                                          !CMISS_SOLVER_SUPERLU_LIBRARY   !<SuperLU Solver library.
                                          !CMISS_SOLVER_SPOOLES_LIBRARY !<SPOOLES Solver library.
                                          !CMISS_SOLVER_UMFPACK_LIBRARY   !<UMFPACK Solver library.
                                          !CMISS_SOLVER_LUSOL_LIBRARY     !<LUSOL Solver library.
                                          !CMISS_SOLVER_ESSL_LIBRARY      !<ESSL Solver library.
                                          !CMISS_SOLVER_LAPACK_LIBRARY    !<LAPACK Solver library.
                                          !CMISS_SOLVER_TAO_LIBRARY       !<TAO Solver library.
                                          !CMISS_SOLVER_HYPRE_LIBRARY     !<Hypre Solver library.
    CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
                                        !CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE    !<Direct linear Solver type.
                                        !CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE !<Iterative linear Solver type.
    CALL CMISSProblem_SolversCreateFinish(Problem,Err)

    !Create the problem Solver equations
    CALL CMISSSolver_Initialise(Solver,Err)
    CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
    CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)   
    CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
    CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
    CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
                                                            !CMISS_SOLVER_SPARSE_MATRICES !<Use sparse Solver matrices.
                                                            !CMISS_SOLVER_FULL_MATRICES !<Use fully populated Solver matrices.
    CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
    CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

    !Prescribe boundary conditions
    CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
    CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
    CALL CMISSSolverEquations_BoundaryConditionsAnalytic(SolverEquations,Err)
    CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

    !=SOLVE Problem==================================================================================================================
    !Solve the Problem
    CALL CMISSProblem_Solve(Problem,Err)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC

  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber, &
    & GeneratedMeshUserNumber,ProblemUserNumber)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: CoordinateSystemUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: RegionUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: BasisUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: GeneratedMeshUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: ProblemUserNumber

    CALL CMISSProblem_Destroy(ProblemUserNumber,Err)
    CALL CMISSGeneratedMesh_Destroy(RegionUserNumber,GeneratedMeshUserNumber,Err)
    CALL CMISSBasis_Destroy(BasisUserNumber,Err)
    CALL CMISSRegion_Destroy(RegionUserNumber,Err)
    CALL CMISSCoordinateSystem_Destroy(CoordinateSystemUserNumber,Err)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN

END PROGRAM ANALYTIC_LINEAR_ELASTICITYEXAMPLE 

