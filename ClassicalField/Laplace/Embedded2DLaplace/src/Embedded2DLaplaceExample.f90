!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a Laplace equation on a 2D
!> bicubic Hermite surface embedded in 3D space using OpenCMISS calls.
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

!> Main program
PROGRAM Embedded2DLaplaceExample 

  USE OPENCMISS
  USE MPI


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=3.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
 
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: numberOfArguments,argumentLength,status
  INTEGER(CMISSIntg) :: numberGlobalXElements,numberGlobalYElements
  CHARACTER(LEN=255) :: commandArgument,filename

  !CMISS variables

  TYPE(CMISSBasisType) :: basis
  TYPE(CMISSBoundaryConditionsType) :: boundaryConditions
  TYPE(CMISSCoordinateSystemType) :: coordinateSystem,worldCoordinateSystem
  TYPE(CMISSDecompositionType) :: decomposition
  TYPE(CMISSEquationsType) :: equations
  TYPE(CMISSEquationsSetType) :: equationsSet
  TYPE(CMISSFieldType) :: geometricField,equationsSetField,dependentField
  TYPE(CMISSFieldsType) :: fields
  TYPE(CMISSGeneratedMeshType) :: generatedMesh  
  TYPE(CMISSMeshType) :: mesh
  TYPE(CMISSNodesType) :: nodes
  TYPE(CMISSProblemType) :: problem
  TYPE(CMISSRegionType) :: region,worldRegion
  TYPE(CMISSSolverType) :: solver
  TYPE(CMISSSolverEquationsType) :: solverEquations
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(CMISSIntg) :: equationsSetIndex
  INTEGER(CMISSIntg) :: firstNodeNumber,lastNodeNumber
  INTEGER(CMISSIntg) :: firstNodeDomain,lastNodeDomain
  INTEGER(CMISSIntg) :: err
  
  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(numberOfArguments >= 2) THEN
    !If we have enough arguments then use the first one for setting up the problem. The subsequent arguments may be used to
    !pass flags to, say, PETSc.
    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 1.")
    READ(commandArgument(1:argumentLength),*) numberGlobalXElements
    CALL GET_COMMAND_ARGUMENT(2,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 2.")
    READ(commandArgument(1:argumentLength),*) numberGlobalYElements
  ELSE
    !If there are not enough arguments default the problem specification 
    numberGlobalXElements=4
    numberGlobalYElements=4
  ENDIF
  
  !Intialise OpenCMISS
  CALL CMISSInitialise(worldCoordinateSystem,worldRegion,err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,err)

  CALL CMISSRandomSeedsSet(9999,err)
  
  CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["COORDINATE_METRICS_CALCULATE"],err)

  WRITE(filename,'(A,"_",I0,"_",I0)') "Embedded2DLaplace",numberGlobalXElements,numberGlobalYElements
  
  CALL CMISSOutputSetOn(filename,err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(numberOfComputationalNodes,err)
  CALL CMISSComputationalNodeNumberGet(computationalNodeNumber,err)
    
  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(coordinateSystem,err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,coordinateSystem,err)
  !Set the coordinate system to be 3D
  CALL CMISSCoordinateSystem_DimensionSet(coordinateSystem,3,err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(coordinateSystem,err)

  !Start the creation of the region
  CALL CMISSRegion_Initialise(region,err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,worldRegion,region,err)
  CALL CMISSRegion_LabelSet(region,"LaplaceRegion",err)
  !Set the regions coordinate system to the 3D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(region,err)

  !Start the creation of a cubic Hermite basis 
  CALL CMISSBasis_Initialise(basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,basis,err)
  CALL CMISSBasis_TypeSet(basis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CALL CMISSBasis_NumberOfXiSet(basis,2,Err)
  CALL CMISSBasis_InterpolationXiSet(basis,[CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION,CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION],err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(basis,[4,4],err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(basis,err)
   
  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_Initialise(generatedMesh,err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(generatedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(generatedMesh,basis,err)   
  !Define the mesh on the region
  CALL CMISSGeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
  CALL CMISSGeneratedMesh_NumberOfElementsSet(generatedMesh,[numberGlobalXElements,numberGlobalYElements],err)
  CALL CMISSGeneratedMesh_BaseVectorsSet(generatedMesh,RESHAPE([WIDTH,HEIGHT,LENGTH,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP],[3,2]),err)
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(mesh,err)
  CALL CMISSGeneratedMesh_CreateFinish(generatedMesh,MeshUserNumber,mesh,err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(decomposition,err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,mesh,decomposition,err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,err)
  CALL CMISSDecomposition_NumberOfDomainsSet(decomposition,numberOfComputationalNodes,err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(decomposition,err)
 
  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(geometricField,err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,region,geometricField,err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(geometricField,decomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(geometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,err)
  CALL CMISSField_ComponentMeshComponentSet(geometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,err)
  CALL CMISSField_ComponentMeshComponentSet(geometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,1,err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)
  
  !Get the node domains
  FirstNodeNumber=1
  CALL CMISSNodes_Initialise(nodes,err)
  CALL CMISSRegion_NodesGet(region,nodes,err)
  CALL CMISSNodes_NumberOfNodesGet(nodes,lastNodeNumber,err)
  CALL CMISSDecomposition_NodeDomainGet(decomposition,firstNodeNumber,1,firstNodeDomain,err)
  CALL CMISSDecomposition_NodeDomainGet(decomposition,lastNodeNumber,1,lastNodeDomain,err)
  
  !Create the Standard Laplace Equations set
  CALL CMISSEquationsSet_Initialise(equationsSet,err)
  CALL CMISSField_Initialise(equationsSetField,err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,region,geometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMISS_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE,EquationsSetFieldUserNumber, &
    & equationsSetField,equationsSet,err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(equationsSet,err)

  !Create the equations set dependent field variables
  CALL CMISSField_Initialise(dependentField,err)
  CALL CMISSEquationsSet_DependentCreateStart(equationsSet,DependentFieldUserNumber,dependentField,err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(equationsSet,err)

  !Initialise the field with an initial guess
  CALL CMISSField_ComponentValuesInitialise(dependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.5_CMISSDP, &
    & err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(equations,err)
  CALL CMISSEquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type
  !CALL CMISSEquations_SparsityTypeSet(equations,CMISS_EQUATIONS_SPARSE_MATRICES,err)
  CALL CMISSEquations_SparsityTypeSet(equations,CMISS_EQUATIONS_FULL_MATRICES,err)
  !Set the equations set output
  !CALL CMISSEquations_OutputTypeSet(equations,CMISS_EQUATIONS_NO_OUTPUT,err)
  !CALL CMISSEquations_OutputTypeSet(equations,CMISS_EQUATIONS_TIMING_OUTPUT,err)
  !CALL CMISSEquations_OutputTypeSet(equations,CMISS_EQUATIONS_MATRIX_OUTPUT,err)
  CALL CMISSEquations_OutputTypeSet(equations,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(equationsSet,err)
  
  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(problem,err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,problem,err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblem_SpecificationSet(problem,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS,CMISS_PROBLEM_LAPLACE_EQUATION_TYPE, &
    & CMISS_PROBLEM_STANDARD_LAPLACE_SUBTYPE,err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(problem,err)

  !Start the creation of the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(problem,err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(problem,err)
 
  !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(solver,Err)
  CALL CMISSProblem_SolversCreateStart(problem,err)
  CALL CMISSProblem_SolverGet(problem,CMISS_CONTROL_LOOP_NODE,1,solver,err)
  !CALL CMISSSolver_OutputTypeSet(solver,CMISS_SOLVER_NO_OUTPUT,err)
  !CALL CMISSSolver_OutputTypeSet(solver,CMISS_SOLVER_PROGRESS_OUTPUT,err)
  !CALL CMISSSolver_OutputTypeSet(solver,CMISS_SOLVER_TIMING_OUTPUT,err)
  CALL CMISSSolver_OutputTypeSet(solver,CMISS_SOLVER_SOLVER_OUTPUT,err)
  !CALL CMISSSolver_OutputTypeSet(solver,CMISS_SOLVER_MATRIX_OUTPUT,err)
  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(problem,err)

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(solver,Err)
  CALL CMISSSolverEquations_Initialise(solverEquations,err)
  CALL CMISSProblem_SolverEquationsCreateStart(problem,err)
  !Get the solve equations
  CALL CMISSProblem_SolverGet(problem,CMISS_CONTROL_LOOP_NODE,1,solver,err)
  CALL CMISSSolver_SolverEquationsGet(solver,solverEquations,err)
  !Set the solver equations sparsity
  !CALL CMISSSolverEquations_SparsityTypeSet(solverEquations,CMISS_SOLVER_SPARSE_MATRICES,err)
  CALL CMISSSolverEquations_SparsityTypeSet(solverEquations,CMISS_SOLVER_FULL_MATRICES,err)  
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(problem,err)

  !Start the creation of the equations set boundary conditions
  CALL CMISSBoundaryConditions_Initialise(boundaryConditions,err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  !Set the first node to 0.0 and the last node to 1.0
  IF(firstNodeDomain==computationalNodeNumber) THEN
    CALL CMISSBoundaryConditions_SetNode(boundaryConditions,dependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,firstNodeNumber,1, &
      & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,err)
  ENDIF
  IF(lastNodeDomain==computationalNodeNumber) THEN
    CALL CMISSBoundaryConditions_SetNode(boundaryConditions,dependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,lastNodeNumber,1, &
      & CMISS_BOUNDARY_CONDITION_FIXED,1.0_CMISSDP,err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Solve the problem
  CALL CMISSProblem_Solve(problem,err)

  !Export results
  CALL CMISSFields_Initialise(fields,err)
  CALL CMISSFields_Create(region,Fields,err)
  CALL CMISSFields_NodesExport(fields,"Embedded2DLaplace","FORTRAN",err)
  CALL CMISSFields_ElementsExport(fields,"Embedded2DLaplace","FORTRAN",err)
  CALL CMISSFields_Finalise(fields,err)
  
  !Finialise CMISS
  CALL CMISSFinalise(err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
CONTAINS

  SUBROUTINE HandleError(errorString)

    CHARACTER(LEN=*), INTENT(IN) :: errorString

    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP

  END SUBROUTINE HandleError
    
END PROGRAM Embedded2DLaplaceExample
