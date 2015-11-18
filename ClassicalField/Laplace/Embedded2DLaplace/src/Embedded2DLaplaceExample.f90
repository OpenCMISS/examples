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

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=3.0_CMISSRP

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

  TYPE(cmfe_BasisType) :: basis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem,worldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) :: geometricField,equationsSetField,dependentField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh  
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_NodesType) :: nodes
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: solver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  
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
  CALL cmfe_Initialise(worldCoordinateSystem,worldRegion,err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)

  CALL cmfe_RandomSeedsSet(9999,err)
  
  CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["COORDINATE_METRICS_CALCULATE"],err)

  WRITE(filename,'(A,"_",I0,"_",I0)') "Embedded2DLaplace",numberGlobalXElements,numberGlobalYElements
  
  CALL cmfe_OutputSetOn(filename,err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(numberOfComputationalNodes,err)
  CALL cmfe_ComputationalNodeNumberGet(computationalNodeNumber,err)
    
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,coordinateSystem,err)
  !Set the coordinate system to be 3D
  CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,3,err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !Start the creation of the region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,worldRegion,region,err)
  CALL cmfe_Region_LabelSet(region,"LaplaceRegion",err)
  !Set the regions coordinate system to the 3D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(region,err)

  !Start the creation of a cubic Hermite basis 
  CALL cmfe_Basis_Initialise(basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,basis,err)
  CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CALL cmfe_Basis_NumberOfXiSet(basis,2,Err)
  CALL cmfe_Basis_InterpolationXiSet(basis,[CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION],err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[4,4],err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(basis,err)
   
  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,basis,err)   
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberGlobalXElements,numberGlobalYElements],err)
  CALL cmfe_GeneratedMesh_BaseVectorsSet(generatedMesh,RESHAPE([WIDTH,HEIGHT,LENGTH,1.0_CMISSRP,1.0_CMISSRP,1.0_CMISSRP],[3,2]),err)
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MeshUserNumber,mesh,err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,mesh,decomposition,err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(decomposition,numberOfComputationalNodes,err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)
 
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(geometricField,err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,region,geometricField,err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(geometricField,decomposition,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,err)
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)
  
  !Get the node domains
  FirstNodeNumber=1
  CALL cmfe_Nodes_Initialise(nodes,err)
  CALL cmfe_Region_NodesGet(region,nodes,err)
  CALL cmfe_Nodes_NumberOfNodesGet(nodes,lastNodeNumber,err)
  CALL cmfe_Decomposition_NodeDomainGet(decomposition,firstNodeNumber,1,firstNodeDomain,err)
  CALL cmfe_Decomposition_NodeDomainGet(decomposition,lastNodeNumber,1,lastNodeDomain,err)
  
  !Create the Standard Laplace Equations set
  CALL cmfe_EquationsSet_Initialise(equationsSet,err)
  CALL cmfe_Field_Initialise(equationsSetField,err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,region,geometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMFE_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE],EquationsSetFieldUserNumber, &
    & equationsSetField,equationsSet,err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(dependentField,err)
  CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,DependentFieldUserNumber,dependentField,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)

  !Initialise the field with an initial guess
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.5_CMISSRP, &
    & err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type
  !CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_FULL_MATRICES,err)
  !Set the equations set output
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_TIMING_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_MATRIX_OUTPUT,err)
  CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)
  
  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(problem,err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS,CMFE_PROBLEM_LAPLACE_EQUATION_TYPE, &
    & CMFE_PROBLEM_STANDARD_LAPLACE_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(problem,err)

  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)
 
  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(solver,Err)
  CALL cmfe_Problem_SolversCreateStart(problem,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_NO_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_TIMING_OUTPUT,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_SOLVER_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_MATRIX_OUTPUT,err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(problem,err)

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(solver,Err)
  CALL cmfe_SolverEquations_Initialise(solverEquations,err)
  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
  !Get the solve equations
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
  !Set the solver equations sparsity
  !CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)
  CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_FULL_MATRICES,err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

  !Start the creation of the equations set boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  !Set the first node to 0.0 and the last node to 1.0
  IF(firstNodeDomain==computationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,firstNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
  ENDIF
  IF(lastNodeDomain==computationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,lastNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Solve the problem
  CALL cmfe_Problem_Solve(problem,err)

  !Export results
  CALL cmfe_Fields_Initialise(fields,err)
  CALL cmfe_Fields_Create(region,Fields,err)
  CALL cmfe_Fields_NodesExport(fields,"Embedded2DLaplace","FORTRAN",err)
  CALL cmfe_Fields_ElementsExport(fields,"Embedded2DLaplace","FORTRAN",err)
  CALL cmfe_Fields_Finalise(fields,err)
  
  !Finialise CMISS
  CALL cmfe_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
CONTAINS

  SUBROUTINE HandleError(errorString)

    CHARACTER(LEN=*), INTENT(IN) :: errorString

    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP

  END SUBROUTINE HandleError
    
END PROGRAM Embedded2DLaplaceExample
