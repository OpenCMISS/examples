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

  REAL(CMISSDP), PARAMETER :: ORIGIN(3)=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
  REAL(CMISSDP), PARAMETER :: LENGTH=20.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=20.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: HEIGHT=5.0_CMISSDP

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

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  !CALL CMISSDiagnosticsSetOn(CMISSFromDiagType,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_FINITE_ELEMENT_CALCULATE"/),Err)

  CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT(1,0,0,"LinearLagrange")
  CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT(1,1,0,"BiLinearLagrange")
  CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT(1,1,1,"TriLinearLagrange")
  !CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_QUADRATIC_LAGRANGE_EXPORT(1,0,0,"QuadraticLagrange")
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
      & CMISSBasisLinearLagrangeInterpolation,DependentField)

    CALL CMISSAnalyticAnalysisOutput(DependentField,OutputFilename,Err)
    
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
    TYPE(CMISSFieldType) :: DependentField

    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC(NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements, &
      & CMISSBasisQuadraticLagrangeInterpolation,DependentField)

    CALL CMISSAnalyticAnalysisOutput(DependentField,OutputFilename,Err)
    
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
    TYPE(CMISSFieldType) :: DependentField

    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC(NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements, &
      & CMISSBasisCubicLagrangeInterpolation,DependentField)

    CALL CMISSAnalyticAnalysisOutput(DependentField,OutputFilename,Err)
    
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

    IF((NumberGlobalYElements == 0) .AND. (NumberGlobalZElements == 0)) THEN
      NumberOfXi = 1
      EquationSetSubtype = CMISSEquationsSetOneDimensionalSubtype
      AnalyticFunction=CMISSEquationsSetLinearElasticityOneDim1
      !Prescribe material properties Area,E1
      FieldMaterialNumberOfComponents = 2 !Young's Modulus & Poisson's Ratio
      MaterialParameters = (/WIDTH*HEIGHT,10.0E3_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
    ELSEIF (NumberGlobalZElements == 0) THEN
      NumberOfXi = 2
      EquationSetSubtype = CMISSEquationsSetPlaneStressSubtype
      AnalyticFunction=CMISSEquationsSetLinearElasticityTwoDim1
      !Prescribe material properties h,E1,v12
      FieldMaterialNumberOfComponents = 3 !Young's Modulus & Poisson's Ratio
      MaterialParameters = (/HEIGHT,10.0E3_CMISSDP,0.3_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
    ELSE
      NumberOfXi = 3
      EquationSetSubtype = CMISSEquationsSetThreeDimensionalSubtype
      AnalyticFunction=CMISSEquationsSetLinearElasticityThreeDim1
      !Prescribe material properties E1,E2,E3 & v13,v23,v12
      FieldMaterialNumberOfComponents = 6 !Young's Modulus & Poisson's Ratio
      MaterialParameters = (/10.0E3_CMISSDP,10.0E3_CMISSDP,10.0E3_CMISSDP,0.3_CMISSDP,0.3_CMISSDP,0.3_CMISSDP/)
    ENDIF
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
    CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
    CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
    CALL CMISSCoordinateSystemTypeSet(CoordinateSystem,CMISSCoordinateRectangularCartesianType,Err)
    CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NumberOfXi,Err)
    CALL CMISSCoordinateSystemOriginSet(CoordinateSystem,ORIGIN,Err)
    CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

    !Create a region and assign the CS to the region
    CALL CMISSRegionTypeInitialise(Region,Err)
    CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
    CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
    CALL CMISSRegionCreateFinish(Region,Err)

    CALL CMISSBasisTypeInitialise(Basis,Err)
    CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
    CALL CMISSBasisTypeSet(Basis,CMISSBasisLagrangeHermiteTPType,Err)
    CALL CMISSBasisNumberOfXiSet(Basis,NumberOfXi,Err)
    CALL CMISSBasisInterpolationXiSet(Basis,Interpolation(1:NumberOfXi),Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis,NumberOfGaussPoints(1:NumberOfXi),Err)
    CALL CMISSBasisCreateFinish(Basis,Err)

    !Start the creation of a generated Mesh in the Region
    CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
    CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
    CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,1,Err)
    CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)

    !Define the Mesh on the Region
    CALL CMISSGeneratedMeshOriginSet(GeneratedMesh,ORIGIN(1:NumberOfXi),Err)
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,MeshDimensions(1:NumberOfXi),Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,NumberOfElements(1:NumberOfXi),Err)
    CALL CMISSMeshTypeInitialise(Mesh,Err)
    CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,1,Mesh,Err)

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
    CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)  
    CALL CMISSFieldNumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
    CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)  
    DO FieldComponentIndex=1,FieldGeometryNumberOfComponents
      CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,FieldComponentIndex,1,Err)
    ENDDO !FieldComponentIndex
    CALL CMISSFieldCreateFinish(GeometricField,Err)

    !Update the geometric field parameters
    CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

    !Create a dependent field with two variables and three components
    CALL CMISSFieldTypeInitialise(DependentField,Err)
    CALL CMISSFieldCreateStart(FieldDependentUserNumber,Region,DependentField,Err)
    CALL CMISSFieldTypeSet(DependentField,CMISSFieldGeneralType,Err)  
    CALL CMISSFieldMeshDecompositionSet(DependentField,Decomposition,Err)
    CALL CMISSFieldGeometricFieldSet(DependentField,GeometricField,Err) 
    CALL CMISSFieldDependentTypeSet(DependentField,CMISSFieldDependentType,Err) 
    CALL CMISSFieldNumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
    CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldUVariableType,FieldDependentNumberOfComponents,Err)
    CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldDelUDelNVariableType,FieldDependentNumberOfComponents,Err)
    DO FieldComponentIndex=1,FieldDependentNumberOfComponents
      CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,FieldComponentIndex,1,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,FieldComponentIndex,1,Err)
    ENDDO !FieldComponentIndex
    CALL CMISSFieldCreateFinish(DependentField,Err)

    !Create a material field and attach it to the geometric field
    CALL CMISSFieldTypeInitialise(MaterialField,Err)
    CALL CMISSFieldCreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
    CALL CMISSFieldTypeSet(MaterialField,CMISSFieldMaterialType,Err)
    CALL CMISSFieldMeshDecompositionSet(MaterialField,Decomposition,Err)
    CALL CMISSFieldGeometricFieldSet(MaterialField,GeometricField,Err)
    CALL CMISSFieldNumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
    CALL CMISSFieldNumberOfComponentsSet(MaterialField,CMISSFieldUVariableType,FieldMaterialNumberOfComponents,Err)
    DO FieldComponentIndex=1,FieldMaterialNumberOfComponents
      CALL CMISSFieldComponentMeshComponentSet(MaterialField,CMISSFieldUVariableType,FieldComponentIndex,1,Err)
    ENDDO !FieldComponentIndex
    CALL CMISSFieldCreateFinish(MaterialField,Err)

    !Set isotropic elasticity material parameters - Young's Modulus & Poisson's Ratio
    DO FieldComponentIndex=1,FieldMaterialNumberOfComponents
      CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,FieldComponentIndex, &
        & MaterialParameters(FieldComponentIndex),Err)
    ENDDO !FieldComponentIndex

    !Create a Elasticity Class, Linear Elasticity type, no subtype, EquationsSet
    CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
      CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
CALL CMISSEquationsSetCreateStart(EquationSetUserNumber,Region,GeometricField,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetLinearElasticityType,EquationSetSubtype,EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
    
    CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

    CALL CMISSEquationsSetDependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
    CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

    CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
    CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

    !Create the Equations set analtyic field variables
    CALL CMISSFieldTypeInitialise(AnalyticField,Err)
    CALL CMISSEquationsSetAnalyticCreateStart(EquationsSet,AnalyticFunction,FieldAnalyticUserNumber,AnalyticField,Err)
    CALL CMISSEquationsSetAnalyticCreateFinish(EquationsSet,Err)

    !Create the equations set equations
    CALL CMISSEquationsTypeInitialise(Equations,Err)
    CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
    CALL CMISSEquationsSparsityTypeSet(EQUATIONS,CMISSEquationsSparseMatrices,Err)
                                                !CMISSEquationsSparseMatrices=1 !<Use sparse matrices for the equations.
                                                !CMISSEquationsFullMatrices=2 !<Use fully populated matrices for the equations. 
    CALL CMISSEquationsOutputTypeSet(EQUATIONS,CMISSEquationsElementMatrixOutput,Err)
                                              !CMISSEquationsNoOutput !<No output from the equations.
                                              !CMISSEquationsTimingOutput !<Timing information output.
                                              !CMISSEquationsMatrixOutput !<All below and equation matrices output.
                                              !CMISSEquationsElementMatrixOutput !<All below and Element matrices output.
    CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)
    
    !Define the problem
    CALL CMISSProblemTypeInitialise(Problem,Err)
    CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
    CALL CMISSProblemSpecificationSet(Problem,CMISSProblemElasticityClass,CMISSProblemLinearElasticityType, &
      & CMISSProblemNoSubtype,Err)
    CALL CMISSProblemCreateFinish(Problem,Err)

    !Create the problem control loop
    CALL CMISSProblemControlLoopCreateStart(Problem,Err)
    CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

    !Start the creation of the Problem Solvers
    !Create the problem Solvers
    CALL CMISSSolverTypeInitialise(Solver,Err)
    CALL CMISSProblemSolversCreateStart(Problem,Err)
    CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
    CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
                                        !CMISSSolverNoOutput !<No output from the Solver routines. \see OPENCMISS_SolverOutputTypes,OPENCMISS
                                        !CMISSSolverProgressOutput !<Progress output from Solver routines.
                                        !CMISSSolverTimingOutput !<Timing output from the Solver routines plus below.
                                        !CMISSSolverSolverOutput !<Solver specific output from the Solver routines plus below.
                                        !CMISSSolverSolverMatrixOutput !<Solver matrices output from the Solver routines plus below.
    CALL CMISSSolverLibraryTypeSet(Solver,CMISSSolverPETScLibrary,Err)
                                          !CMISSSolverCMISSLibrary     !<CMISS (internal) Solver library.
                                          !CMISSSolverPETScLibrary     !<PETSc Solver library.
                                          !CMISSSolverMUMPSLibrary     !<MUMPS Solver library.
                                          !CMISSSolverSuperLULibrary   !<SuperLU Solver library.
                                          !CMISSSolverSpoolesLULibrary !<SPOOLES Solver library.
                                          !CMISSSolverUMFPACKLibrary   !<UMFPACK Solver library.
                                          !CMISSSolverLUSOLLibrary     !<LUSOL Solver library.
                                          !CMISSSolverESSLLibrary      !<ESSL Solver library.
                                          !CMISSSolverLAPACKLibrary    !<LAPACK Solver library.
                                          !CMISSSolverTAOLibrary       !<TAO Solver library.
                                          !CMISSSolverHypreLibrary     !<Hypre Solver library.
    CALL CMISSSolverLinearTypeSet(Solver,CMISSSolverLinearDirectSolveType,Err)
                                        !CMISSSolverLinearDirectSolveType    !<Direct linear Solver type.
                                        !CMISSSolverLinearIterativeSolveType !<Iterative linear Solver type.
    CALL CMISSProblemSolversCreateFinish(Problem,Err)

    !Create the problem Solver equations
    CALL CMISSSolverTypeInitialise(Solver,Err)
    CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
    CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)   
    CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
    CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
    CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
                                                            !CMISSSolverEquationsSparseMatrices !<Use sparse Solver matrices.
                                                            !CMISSSolverEquationsFullMatrices !<Use fully populated Solver matrices.
    CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
    CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

    !Prescribe boundary conditions
    CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
    CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
    CALL CMISSProblemSolverEquationsBoundaryConditionsAnalytic(SolverEquations,Err)
    CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

    !=SOLVE Problem==================================================================================================================
    !Solve the Problem
    CALL CMISSProblemSolve(Problem,Err)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC

  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber, &
    & GeneratedMeshUserNumber,ProblemUserNumber)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: CoordinateSystemUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: RegionUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: BasisUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: GeneratedMeshUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: ProblemUserNumber

    CALL CMISSProblemDestroy(ProblemUserNumber,Err)
    CALL CMISSGeneratedMeshDestroy(RegionUserNumber,GeneratedMeshUserNumber,Err)
    CALL CMISSBasisDestroy(BasisUserNumber,Err)
    CALL CMISSRegionDestroy(RegionUserNumber,Err)
    CALL CMISSCoordinateSystemDestroy(CoordinateSystemUserNumber,Err)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN

END PROGRAM ANALYTIC_LINEAR_ELASTICITYEXAMPLE 

