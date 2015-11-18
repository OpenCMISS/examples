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
#ifndef NOMPIMOD
  USE MPI
#endif
  USE OpenCMISS
  USE OpenCMISS_Iron
  USE TEST_FRAMEWORK_ROUTINES

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  REAL(CMISSRP), PARAMETER :: ORIGIN(3)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP]
  REAL(CMISSRP), PARAMETER :: LENGTH=20.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=20.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: HEIGHT=5.0_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: NumberOfDomains=1

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber = 1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariables=2

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariables=1

  INTEGER(CMISSIntg), PARAMETER :: FieldAnalyticUserNumber=4

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program types

  TYPE(cmfe_RegionType) :: WorldRegion
  TYPE(cmfe_CoordinateSystemType) :: WorldCoordinateSystem

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
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_FINITE_ELEMENT_CALCULATE"],Err)

  CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT(1,0,0,"LinearLagrange")
  CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT(1,1,0,"BiLinearLagrange")
  CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT(1,1,1,"TriLinearLagrange")
  !CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_QUADRATIC_LAGRANGE_EXPORT(1,0,0,"QuadraticLagrange")
  CALL cmfe_Finalise(Err)

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
    TYPE(cmfe_FieldType) :: DependentField

    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC(NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements, &
      & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,DependentField)

    CALL cmfe_AnalyticAnalysis_Output(DependentField,OutputFilename,Err)
    
    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber, &
      & GeneratedMeshUserNumber,ProblemUserNumber)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT

  !
  !================================================================================================================================
  !  
    !>Check if the convergence of quadratic langrange interpolation is expected.
  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_QUADRATIC_LAGRANGE_EXPORT(NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements,OutputFilename)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NumberGlobalXElements !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NumberGlobalYElements !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NumberGlobalZElements !<increment interval number of elements per axis
    CHARACTER(LEN=*), INTENT(IN) :: OutputFilename !<The Error condition string
    !Local Variables
    TYPE(cmfe_FieldType) :: DependentField

    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC(NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements, &
      & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,DependentField)

    CALL cmfe_AnalyticAnalysis_Output(DependentField,OutputFilename,Err)
    
    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber, &
      & GeneratedMeshUserNumber,ProblemUserNumber)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_QUADRATIC_LAGRANGE_EXPORT

  !
  !================================================================================================================================
  !  
    !>Check if the convergence of cubic langrange interpolation is expected.
  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_CUBIC_LAGRANGE_EXPORT(NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements,OutputFilename)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NumberGlobalXElements !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NumberGlobalYElements !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NumberGlobalZElements !<increment interval number of elements per axis
    CHARACTER(LEN=*), INTENT(IN) :: OutputFilename !<The Error condition string
    !Local Variables
    TYPE(cmfe_FieldType) :: DependentField

    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC(NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements, &
      & CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION,DependentField)

    CALL cmfe_AnalyticAnalysis_Output(DependentField,OutputFilename,Err)
    
    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber, &
      & GeneratedMeshUserNumber,ProblemUserNumber)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_CUBIC_LAGRANGE_EXPORT

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
    TYPE(cmfe_FieldType) :: DependentField

    !Program variables
    REAL(CMISSRP) :: MeshDimensions(3),MaterialParameters(6)
    INTEGER(CMISSIntg) :: AnalyticFunction,Interpolation(3),NumberOfGaussPoints(3),EquationSetSubtype
    INTEGER(CMISSIntg) :: FieldGeometryNumberOfComponents,FieldDependentNumberOfComponents,NumberOfElements(3)
    INTEGER(CMISSIntg) :: MPI_IERROR
    INTEGER(CMISSIntg) :: EquationsSetIndex,FieldComponentIndex,FieldMaterialNumberOfComponents,NumberOfXi
    INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber

    !CMISS variables

    TYPE(cmfe_BasisType) :: Basis
    TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
    TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
    TYPE(cmfe_DecompositionType) :: Decomposition
    TYPE(cmfe_EquationsType) :: Equations
    TYPE(cmfe_EquationsSetType) :: EquationsSet
    TYPE(cmfe_FieldType) :: AnalyticField,EquationsSetField,GeometricField,MaterialField
    TYPE(cmfe_MeshType) :: Mesh
    TYPE(cmfe_ProblemType) :: Problem
    TYPE(cmfe_RegionType) :: Region
    TYPE(cmfe_SolverType) :: Solver
    TYPE(cmfe_SolverEquationsType) :: SolverEquations
    TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions

    IF((NumberGlobalYElements == 0) .AND. (NumberGlobalZElements == 0)) THEN
      NumberOfXi = 1
      EquationSetSubtype = CMFE_EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE
      AnalyticFunction=CMFE_EQUATIONS_SET_LINEAR_ELASTICITY_ONE_DIM_1
      !Prescribe material properties Area,E1
      FieldMaterialNumberOfComponents = 2 !Young's Modulus & Poisson's Ratio
      MaterialParameters = [WIDTH*HEIGHT,10.0E3_CMISSRP,0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP]
    ELSEIF (NumberGlobalZElements == 0) THEN
      NumberOfXi = 2
      EquationSetSubtype = CMFE_EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE
      AnalyticFunction=CMFE_EQUATIONS_SET_LINEAR_ELASTICITY_TWO_DIM_1
      !Prescribe material properties h,E1,v12
      FieldMaterialNumberOfComponents = 3 !Young's Modulus & Poisson's Ratio
      MaterialParameters = [HEIGHT,10.0E3_CMISSRP,0.3_CMISSRP,0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP]
    ELSE
      NumberOfXi = 3
      EquationSetSubtype = CMFE_EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE
      AnalyticFunction=CMFE_EQUATIONS_SET_LINEAR_ELASTICITY_THREE_DIM_1
      !Prescribe material properties E1,E2,E3 & v13,v23,v12
      FieldMaterialNumberOfComponents = 6 !Young's Modulus & Poisson's Ratio
      MaterialParameters = [10.0E3_CMISSRP,10.0E3_CMISSRP,10.0E3_CMISSRP,0.3_CMISSRP,0.3_CMISSRP,0.3_CMISSRP]
    ENDIF
    Interpolation = [InterpolationSpecifications,InterpolationSpecifications,InterpolationSpecifications]
    NumberOfElements = [NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements]
    MeshDimensions = [LENGTH,WIDTH,HEIGHT]
    NumberOfGaussPoints = [4,4,4]
    FieldGeometryNumberOfComponents=NumberOfXi
    FieldDependentNumberOfComponents=NumberOfXi

    !Get the number of computational nodes and this computational node number
    CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
    CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

    !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
    CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

    !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
    CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
    CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
    CALL cmfe_CoordinateSystem_TypeSet(CoordinateSystem,CMFE_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,NumberOfXi,Err)
    CALL cmfe_CoordinateSystem_OriginSet(CoordinateSystem,ORIGIN,Err)
    CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

    !Create a region and assign the CS to the region
    CALL cmfe_Region_Initialise(Region,Err)
    CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
    CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
    CALL cmfe_Region_CreateFinish(Region,Err)

    CALL cmfe_Basis_Initialise(Basis,Err)
    CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
    CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CALL cmfe_Basis_NumberOfXiSet(Basis,NumberOfXi,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis,Interpolation(1:NumberOfXi),Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,NumberOfGaussPoints(1:NumberOfXi),Err)
    CALL cmfe_Basis_CreateFinish(Basis,Err)

    !Start the creation of a generated Mesh in the Region
    CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
    CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
    CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,1,Err)
    CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)

    !Define the Mesh on the Region
    CALL cmfe_GeneratedMesh_OriginSet(GeneratedMesh,ORIGIN(1:NumberOfXi),Err)
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,MeshDimensions(1:NumberOfXi),Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,NumberOfElements(1:NumberOfXi),Err)
    CALL cmfe_Mesh_Initialise(Mesh,Err)
    CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,1,Mesh,Err)

    !Create a decomposition
    CALL cmfe_Decomposition_Initialise(Decomposition,Err)
    CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
    CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
    CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

    !Create a field to put the geometry (defualt is geometry)
    CALL cmfe_Field_Initialise(GeometricField,Err)
    CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
    CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
    CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)  
    CALL cmfe_Field_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
    CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
    DO FieldComponentIndex=1,FieldGeometryNumberOfComponents
      CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,FieldComponentIndex,1,Err)
    ENDDO !FieldComponentIndex
    CALL cmfe_Field_CreateFinish(GeometricField,Err)

    !Update the geometric field parameters
    CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

    !Create a dependent field with two variables and three components
    CALL cmfe_Field_Initialise(DependentField,Err)
    CALL cmfe_Field_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
    CALL cmfe_Field_TypeSet(DependentField,CMFE_FIELD_GENERAL_TYPE,Err)  
    CALL cmfe_Field_MeshDecompositionSet(DependentField,Decomposition,Err)
    CALL cmfe_Field_GeometricFieldSet(DependentField,GeometricField,Err) 
    CALL cmfe_Field_DependentTypeSet(DependentField,CMFE_FIELD_DEPENDENT_TYPE,Err) 
    CALL cmfe_Field_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
    CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
    CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
    DO FieldComponentIndex=1,FieldDependentNumberOfComponents
      CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,FieldComponentIndex,1,Err)
      CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,FieldComponentIndex,1,Err)
    ENDDO !FieldComponentIndex
    CALL cmfe_Field_CreateFinish(DependentField,Err)

    !Create a material field and attach it to the geometric field
    CALL cmfe_Field_Initialise(MaterialField,Err)
    CALL cmfe_Field_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
    CALL cmfe_Field_TypeSet(MaterialField,CMFE_FIELD_MATERIAL_TYPE,Err)
    CALL cmfe_Field_MeshDecompositionSet(MaterialField,Decomposition,Err)
    CALL cmfe_Field_GeometricFieldSet(MaterialField,GeometricField,Err)
    CALL cmfe_Field_NumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
    CALL cmfe_Field_NumberOfComponentsSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)
    DO FieldComponentIndex=1,FieldMaterialNumberOfComponents
      CALL cmfe_Field_ComponentMeshComponentSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,FieldComponentIndex,1,Err)
    ENDDO !FieldComponentIndex
    CALL cmfe_Field_CreateFinish(MaterialField,Err)

    !Set isotropic elasticity material parameters - Young's Modulus & Poisson's Ratio
    DO FieldComponentIndex=1,FieldMaterialNumberOfComponents
      CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & FieldComponentIndex, &
        & MaterialParameters(FieldComponentIndex),Err)
    ENDDO !FieldComponentIndex

    !Create a Elasticity Class, Linear Elasticity type, no subtype, EquationsSet
    CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
    CALL cmfe_Field_Initialise(EquationsSetField,Err)
    CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
      & CMFE_EQUATIONS_SET_LINEAR_ELASTICITY_TYPE,EquationSetSubtype],EquationsSetFieldUserNumber,EquationsSetField, &
      & EquationsSet,Err)
    
    CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

    CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
    CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

    CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
    CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

    !Create the Equations set analtyic field variables
    CALL cmfe_Field_Initialise(AnalyticField,Err)
    CALL cmfe_EquationsSet_AnalyticCreateStart(EquationsSet,AnalyticFunction,FieldAnalyticUserNumber,AnalyticField,Err)
    CALL cmfe_EquationsSet_AnalyticCreateFinish(EquationsSet,Err)

    !Create the equations set equations
    CALL cmfe_Equations_Initialise(Equations,Err)
    CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
    CALL cmfe_Equations_SparsityTypeSet(EQUATIONS,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
                                                !CMFE_EQUATIONS_SPARSE_MATRICES=1 !<Use sparse matrices for the equations.
                                                !CMFE_EQUATIONS_FULL_MATRICES=2 !<Use fully populated matrices for the equations. 
    CALL cmfe_Equations_OutputTypeSet(EQUATIONS,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
                                              !CMFE_EQUATIONS_NO_OUTPUT !<No output from the equations.
                                              !CMFE_EQUATIONS_TIMING_OUTPUT !<Timing information output.
                                              !CMFE_EQUATIONS_MATRIX_OUTPUT !<All below and equation matrices output.
                                              !CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT !<All below and Element matrices output.
    CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)
    
    !Define the problem
    CALL cmfe_Problem_Initialise(Problem,Err)
    CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_LINEAR_ELASTICITY_TYPE, &
      & CMFE_PROBLEM_NO_SUBTYPE],Problem,Err)
    CALL cmfe_Problem_CreateFinish(Problem,Err)

    !Create the problem control loop
    CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
    CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

    !Start the creation of the Problem Solvers
    !Create the problem Solvers
    CALL cmfe_Solver_Initialise(Solver,Err)
    CALL cmfe_Problem_SolversCreateStart(Problem,Err)
    CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
    CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
                                        !CMFE_SOLVER_NO_OUTPUT !<No output from the Solver routines. \see OPENCMISS_SolverOutputTypes,OPENCMISS
                                        !CMFE_SOLVER_PROGRESS_OUTPUT !<Progress output from Solver routines.
                                        !CMFE_SOLVER_TIMING_OUTPUT !<Timing output from the Solver routines plus below.
                                        !CMFE_SOLVER_SOLVER_OUTPUT !<Solver specific output from the Solver routines plus below.
                                        !CMFE_SOLVER_MATRIX_OUTPUT !<Solver matrices output from the Solver routines plus below.
    CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_PETSC_LIBRARY,Err)
                                          !CMFE_SOLVER_CMISS_LIBRARY     !<CMISS (internal) Solver library.
                                          !CMFE_SOLVER_PETSC_LIBRARY     !<PETSc Solver library.
                                          !CMFE_SOLVER_MUMPS_LIBRARY     !<MUMPS Solver library.
                                          !CMFE_SOLVER_SUPERLU_LIBRARY   !<SuperLU Solver library.
                                          !CMFE_SOLVER_SPOOLES_LIBRARY !<SPOOLES Solver library.
                                          !CMFE_SOLVER_UMFPACK_LIBRARY   !<UMFPACK Solver library.
                                          !CMFE_SOLVER_LUSOL_LIBRARY     !<LUSOL Solver library.
                                          !CMFE_SOLVER_ESSL_LIBRARY      !<ESSL Solver library.
                                          !CMFE_SOLVER_LAPACK_LIBRARY    !<LAPACK Solver library.
                                          !CMFE_SOLVER_TAO_LIBRARY       !<TAO Solver library.
                                          !CMFE_SOLVER_HYPRE_LIBRARY     !<Hypre Solver library.
    CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
                                        !CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE    !<Direct linear Solver type.
                                        !CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE !<Iterative linear Solver type.
    CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

    !Create the problem Solver equations
    CALL cmfe_Solver_Initialise(Solver,Err)
    CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
    CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)   
    CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
    CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
    CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
                                                            !CMFE_SOLVER_SPARSE_MATRICES !<Use sparse Solver matrices.
                                                            !CMFE_SOLVER_FULL_MATRICES !<Use fully populated Solver matrices.
    CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
    CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

    !Prescribe boundary conditions
    CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
    CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
    CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(SolverEquations,Err)
    CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

    !=SOLVE Problem==================================================================================================================
    !Solve the Problem
    CALL cmfe_Problem_Solve(Problem,Err)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC

  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber, &
    & GeneratedMeshUserNumber,ProblemUserNumber)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: CoordinateSystemUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: RegionUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: BasisUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: GeneratedMeshUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: ProblemUserNumber

    CALL cmfe_Problem_Destroy(ProblemUserNumber,Err)
    CALL cmfe_GeneratedMesh_Destroy(RegionUserNumber,GeneratedMeshUserNumber,Err)
    CALL cmfe_Basis_Destroy(BasisUserNumber,Err)
    CALL cmfe_Region_Destroy(RegionUserNumber,Err)
    CALL cmfe_CoordinateSystem_Destroy(CoordinateSystemUserNumber,Err)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN

END PROGRAM ANALYTIC_LINEAR_ELASTICITYEXAMPLE 

