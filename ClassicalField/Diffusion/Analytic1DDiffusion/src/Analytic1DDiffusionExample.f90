!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a 1D diffusion equation and compare it to the analytic solution using OpenCMISS calls.
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

!> \example ClassicalField/Diffusion/Analytic1DDiffusion/src/Analytic1DDiffusionExample.f90
!! Example program to solve a 1D diffusion equation using OpenCMISS calls.
!!
!! \htmlinclude ClassicalField/Diffusion/Analytic1DDiffusion/history.html
!<

!> Main program
PROGRAM ANALYTIC1DDIFFUSIONEXAMPLE

  USE OPENCMISS

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER ::    PI=3.141592653589793238462643383279502884197_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: NUMBER_GLOBAL_X_ELEMENTS=6
  REAL(CMISSDP), PARAMETER :: LENGTH=3.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: END_TIME=0.1_CMISSDP
  REAL(CMISSDP), PARAMETER :: TIME_STEP=0.01_CMISSDP
  REAL(CMISSDP), PARAMETER :: A=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: B=PI/2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: C=0.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: K=1.0_CMISSDP
  
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=13
  !Program types
  
  !Program variables
  
  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,EquationsSetField,MaterialsField,AnalyticField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver, LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: EquationsSetIndex
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

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  CALL CMISSOutputSetOn("Diffusion1DAnalytic",Err)
  
  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 1D
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,1,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)
  
  !Start the creation of the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Label the Region
  CALL CMISSRegionLabelSet(Region,"Region",Err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)
  
  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasisTypeInitialise(Basis,Err)
  CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
  !Set the basis to be a linear Lagrange basis
  CALL CMISSBasisNumberOfXiSet(Basis,1,Err)
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(BASIS,Err)
  
  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[LENGTH],Err)
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS],Err)
  !Finish the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
  
  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,1,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,1,Err)
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)
       
  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)
  
  !Create the equations_set
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetClassicalFieldClass, &
    & CMISSEquationsSetDiffusionEquationType,CMISSEquationsSetNoSourceDiffusionSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field variables
  CALL CMISSFieldTypeInitialise(MaterialsField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)
  !Set the conductivity
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,K,Err)
 
  !Create the equations set analytic field variables
  CALL CMISSFieldTypeInitialise(AnalyticField,Err)
  CALL CMISSEquationsSetAnalyticCreateStart(EquationsSet,CMISSEquationsSetDiffusionOneDim1, & 
    & AnalyticFieldUserNumber,AnalyticField,Err)
  !Finish the equations set analytic field variables
  CALL CMISSEquationsSetAnalyticCreateFinish(EquationsSet,Err)
  !Set the analytic field parameters
  !Set the multiplicative constant. 
  CALL CMISSFieldComponentValuesInitialise(AnalyticField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,A,Err)
  !Set the phase. 
  CALL CMISSFieldComponentValuesInitialise(AnalyticField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,B,Err)
  !Set the offset. 
  CALL CMISSFieldComponentValuesInitialise(AnalyticField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,C,Err)
  !Set the length to be the length of the mesh
  CALL CMISSFieldComponentValuesInitialise(AnalyticField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,LENGTH,Err)
  
  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)

  !Create the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a No Source Diffusion problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemClassicalFieldClass,CMISSProblemDiffusionEquationType, &
    & CMISSProblemlinearSourceDiffusionSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Create the problem control
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  !Get the control loop
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,0.0_CMISSDP,END_TIME,TIME_STEP,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
  CALL CMISSSolverDynamicLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolverOutputTypeSet(LinearSolver,CMISSSolverProgressOutput,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)  
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Create the equations set boundary conditions
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  CALL CMISSProblemSolverEquationsBoundaryConditionsAnalytic(SolverEquations,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)

  !Output Analytic analysis
  !CALL CMISSEquationsSetAnalyticTimeSet(EquationsSet,END_TIME,Err)
  !CALL CMISSEquationsSetAnalyticEvaluate(EquationsSet,Err)
  CALL CMISSAnalyticAnalysisOutput(DependentField,"Diffusion1DAnalytic",Err)

  !Output fields
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"Diffusion1DAnalytic","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"Diffusion1DAnalytic","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  !Finalise and quit
  CALL CMISSFinalise(Err)
  WRITE(*,'(A)') "Program successfully completed."

  STOP
  
END PROGRAM ANALYTIC1DDIFFUSIONEXAMPLE
