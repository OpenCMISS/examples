/*
 * \file
 * \author Chris Bradley
 * \brief This is an example program to solve Laplace's equation using OpenCMISS calls from C.
 *
 * \section LICENSE
 *
 * Version: MPL 1.1/GPL 2.0/LGPL 2.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 * License for the specific language governing rights and limitations
 * under the License.
 *
 * The Original Code is OpenCMISS
 *
 * The Initial Developer of the Original Code is University of Auckland,
 * Auckland, New Zealand and University of Oxford, Oxford, United
 * Kingdom. Portions created by the University of Auckland and University
 * of Oxford are Copyright (C) 2007 by the University of Auckland and
 * the University of Oxford. All Rights Reserved.
 *
 * Contributor(s):
 *
 * Alternatively, the contents of this file may be used under the terms of
 * either the GNU General Public License Version 2 or later (the "GPL"), or
 * the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
 * in which case the provisions of the GPL or the LGPL are applicable instead
 * of those above. If you wish to allow use of your version of this file only
 * under the terms of either the GPL or the LGPL, and not to allow others to
 * use your version of this file under the terms of the MPL, indicate your
 * decision by deleting the provisions above and replace them with the notice
 * and other provisions required by the GPL or the LGPL. If you do not delete
 * the provisions above, a recipient may use your version of this file under
 * the terms of any one of the MPL, the GPL or the LGPL.
 *
 */
#include <stdlib.h>
#include <stdio.h>

#include "opencmiss.h"

#define STRING_SIZE 20

#define HEIGHT 1.0
#define WIDTH 2.0
#define LENGTH 3.0

#define NUMBER_GLOBAL_X_ELEMENTS 5
#define NUMBER_GLOBAL_Y_ELEMENTS 5
#define NUMBER_GLOBAL_Z_ELEMENTS 5

#define COORDINATE_SYSTEM_USER_NUMBER 1
#define REGION_USER_NUMBER 2
#define BASIS_USER_NUMBER 3
#define GENERATED_MESH_USER_NUMBER 4
#define MESH_USER_NUMBER 5
#define DECOMPOSITION_USER_NUMBER 6
#define GEOMETRIC_FIELD_USER_NUMBER 7
#define DEPENDENT_FIELD_USER_NUMBER 8
#define EQUATIONS_SET_USER_NUMBER 9
#define PROBLEM_USER_NUMBER 10
#define EQUATIONS_SET_FIELD_USER_NUMBER 11

#define MAX_COORDINATES 3

#define CHECK_ERROR(S) \
  if(Err != CMISS_NO_ERROR) { \
    if(Err == CMISS_ERROR_CONVERTING_POINTER) { \
      fprintf(stderr,"Error: %s: Error converting pointer.\n",(S)); \
    } \
    else if(Err == CMISS_POINTER_IS_NULL) { \
      fprintf(stderr,"Error: %s: Pointer is null.\n",(S)); \
    } \
    else if(Err == CMISS_POINTER_NOT_NULL) { \
      fprintf(stderr,"Error: %s: Pointer is not null.\n",(S)); \
    } \
    else if(Err == CMISS_COULD_NOT_ALLOCATE_POINTER) { \
      fprintf(stderr,"Error: %s: Could not allocate pointer.\n",(S)); \
    } \
    exit(Err); \
  }

int main()
{
  CMISSBasisType Basis = (CMISSBasisType)NULL;
  CMISSBoundaryConditionsType BoundaryConditions=(CMISSBoundaryConditionsType)NULL;
  CMISSCoordinateSystemType CoordinateSystem=(CMISSCoordinateSystemType)NULL,WorldCoordinateSystem=(CMISSCoordinateSystemType)NULL;
  CMISSDecompositionType Decomposition=(CMISSDecompositionType)NULL;
  CMISSEquationsType Equations=(CMISSEquationsType)NULL;
  CMISSEquationsSetType EquationsSet=(CMISSEquationsSetType)NULL;
  CMISSFieldType GeometricField=(CMISSFieldType)NULL,DependentField=(CMISSFieldType)NULL,EquationsSetField=(CMISSFieldType)NULL;
  CMISSFieldsType Fields=(CMISSFieldsType)NULL;
  CMISSGeneratedMeshType GeneratedMesh=(CMISSGeneratedMeshType)NULL;
  CMISSMeshType Mesh=(CMISSMeshType)NULL;
  CMISSProblemType Problem=(CMISSProblemType)NULL;
  CMISSRegionType Region=(CMISSRegionType)NULL,WorldRegion=(CMISSRegionType)NULL;
  CMISSSolverType Solver=(CMISSSolverType)NULL;
  CMISSSolverEquationsType SolverEquations=(CMISSSolverEquationsType)NULL;

  int NumberOfComputationalNodes,ComputationalNodeNumber;
  int EquationsSetIndex;
  int FirstNodeNumber,LastNodeNumber;
  int FirstNodeDomain,LastNodeDomain;

  int NumberXiElements[MAX_COORDINATES];
  int ControlLoopIdentifier[1];
  double MeshExtent[MAX_COORDINATES];

  int Err;

  ControlLoopIdentifier[0]=CMISS_CONTROL_LOOP_NODE;

  Err = CMISSCoordinateSystem_Initialise(&WorldCoordinateSystem);
  CHECK_ERROR("Initialising world coordinate system");
  Err = CMISSRegion_Initialise(&WorldRegion);
  CHECK_ERROR("Initialising world region");
  Err = CMISSInitialise(WorldCoordinateSystem,WorldRegion);
  CHECK_ERROR("Initialising CMISS");
  Err = CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR);

  Err = CMISSComputationalNumberOfNodesGet(&NumberOfComputationalNodes);
  Err = CMISSComputationalNodeNumberGet(&ComputationalNodeNumber);

  /* Start the creation of a new RC coordinate system */
  Err = CMISSCoordinateSystem_Initialise(&CoordinateSystem);
  Err = CMISSCoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,CoordinateSystem);
  if(NUMBER_GLOBAL_Z_ELEMENTS==0)
    {
      /* Set the coordinate system to be 2D */
      Err = CMISSCoordinateSystem_DimensionSet(CoordinateSystem,2);
    }
  else
    {
      /* Set the coordinate system to be 3D */
      Err = CMISSCoordinateSystem_DimensionSet(CoordinateSystem,3);
    }
  /* Finish the creation of the coordinate system */
  Err = CMISSCoordinateSystem_CreateFinish(CoordinateSystem);

  /* Start the creation of the region */
  Err = CMISSRegion_Initialise(&Region);
  Err = CMISSRegion_CreateStart(REGION_USER_NUMBER,WorldRegion,Region);
  /* Set the regions coordinate system to the 2D RC coordinate system that we have created */
  Err = CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem);
  /* Finish the creation of the region */
  Err = CMISSRegion_CreateFinish(Region);

  /* Start the creation of a basis (default is trilinear lagrange) */
  Err = CMISSBasis_Initialise(&Basis);
  Err = CMISSBasis_CreateStart(BASIS_USER_NUMBER,Basis);
  if(NUMBER_GLOBAL_Z_ELEMENTS==0)
    {
      /* Set the basis to be a bilinear Lagrange basis */
      Err = CMISSBasis_NumberOfXiSet(Basis,2);
    }
  else
    {
      /* Set the basis to be a trilinear Lagrange basis */
      Err = CMISSBasis_NumberOfXiSet(Basis,3);
    }
  /* Finish the creation of the basis */
  Err = CMISSBasis_CreateFinish(Basis);

  /* Start the creation of a generated mesh in the region */
  Err = CMISSGeneratedMesh_Initialise(&GeneratedMesh);
  Err = CMISSGeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,Region,GeneratedMesh);
  /* Set up a regular x*y*z mesh */
  Err = CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE);
  /* Set the default basis */
  Err = CMISSGeneratedMesh_BasisSet(GeneratedMesh,1,&Basis);
  CHECK_ERROR("Setting mesh basis");
  /* Define the mesh on the region */
  MeshExtent[0]=WIDTH;
  MeshExtent[1]=HEIGHT;
  NumberXiElements[0]=NUMBER_GLOBAL_X_ELEMENTS;
  NumberXiElements[1]=NUMBER_GLOBAL_Y_ELEMENTS;
  if(NUMBER_GLOBAL_Z_ELEMENTS!=0)
    {
      MeshExtent[2]=LENGTH;
      NumberXiElements[2]=NUMBER_GLOBAL_Z_ELEMENTS;
    }
  Err = CMISSGeneratedMesh_ExtentSet(GeneratedMesh,MAX_COORDINATES,MeshExtent);
  Err = CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,MAX_COORDINATES,NumberXiElements);
  /* Finish the creation of a generated mesh in the region */
  Err = CMISSMesh_Initialise(&Mesh);
  /* Finish the creation of a generated mesh in the region */
  Err = CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MESH_USER_NUMBER,Mesh);

  /* Create a decomposition */
  Err = CMISSDecomposition_Initialise(&Decomposition);
  Err = CMISSDecomposition_CreateStart(DECOMPOSITION_USER_NUMBER,Mesh,Decomposition);
  /* Set the decomposition to be a general decomposition with the specified number of domains */
  Err = CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE);
  Err = CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes);
  /* Finish the decomposition */
  Err = CMISSDecomposition_CreateFinish(Decomposition);

  /* Start to create a default (geometric) field on the region */
  Err = CMISSField_Initialise(&GeometricField);
  Err = CMISSField_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,Region,GeometricField);
  /* Set the decomposition to use */
  Err = CMISSField_MeshDecompositionSet(GeometricField,Decomposition);
  /* Set the domain to be used by the field components. */
  Err = CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,1);
  Err = CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1);
  if(NUMBER_GLOBAL_Z_ELEMENTS!=0)
    {
      Err = CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,1);
    }
  /* Finish creating the field */
  Err = CMISSField_CreateFinish(GeometricField);

  /* Update the geometric field parameters */
  Err = CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField);

  /* Create the equations_set */
  Err = CMISSEquationsSet_Initialise(&EquationsSet);
  Err = CMISSField_Initialise(&EquationsSetField);
  Err = CMISSEquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,Region,GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, \
    CMISS_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMISS_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE,EQUATIONS_SET_FIELD_USER_NUMBER, \
    EquationsSetField,EquationsSet);
  CHECK_ERROR("Creating equations set");
  /* Set the equations set to be a standard Laplace problem */
  Err = CMISSEquationsSet_SpecificationSet(EquationsSet,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, \
    CMISS_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMISS_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE);
  /* Finish creating the equations set */
  Err = CMISSEquationsSet_CreateFinish(EquationsSet);

  /* Create the equations set dependent field variables */
  Err = CMISSField_Initialise(&DependentField);
  Err = CMISSEquationsSet_DependentCreateStart(EquationsSet,DEPENDENT_FIELD_USER_NUMBER,DependentField);
  /* Finish the equations set dependent field variables */
  Err = CMISSEquationsSet_DependentCreateFinish(EquationsSet);

  /* Create the equations set equations */
  Err = CMISSEquations_Initialise(&Equations);
  Err = CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations);
  /* Set the equations matrices sparsity type */
  Err = CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES);
  /* Set the equations set output */
  /* Err = CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT); */
  Err = CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_TIMING_OUTPUT);
  /* Err = CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_MATRIX_OUTPUT); */
  /* Err = CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT); */
  /* Finish the equations set equations */
  Err = CMISSEquationsSet_EquationsCreateFinish(EquationsSet);

  /* Start the creation of a problem. */
  Err = CMISSProblem_Initialise(&Problem);
  Err = CMISSProblem_CreateStart(PROBLEM_USER_NUMBER,Problem);
  /* Set the problem to be a standard Laplace problem */
  Err = CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS,CMISS_PROBLEM_LAPLACE_EQUATION_TYPE, \
    CMISS_PROBLEM_STANDARD_LAPLACE_SUBTYPE);
  /* Finish the creation of a problem. */
  Err = CMISSProblem_CreateFinish(Problem);

  /* Start the creation of the problem control loop */
  Err = CMISSProblem_ControlLoopCreateStart(Problem);
  /* Finish creating the problem control loop */
  Err = CMISSProblem_ControlLoopCreateFinish(Problem);

  /* Start the creation of the problem solvers */
  Err = CMISSSolver_Initialise(&Solver);
  Err = CMISSProblem_SolversCreateStart(Problem);
  Err = CMISSProblem_SolverGet(Problem,1,ControlLoopIdentifier,1,Solver);
  /* Err = CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT); */
  /* Err = CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT); */
  /* Err = CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_TIMING_OUTPUT); */
  /* Err = CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_SOLVER_OUTPUT); */
  Err = CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_MATRIX_OUTPUT);
  Err = CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE);
  Err = CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_MUMPS_LIBRARY);
  /* Finish the creation of the problem solver */
  Err = CMISSProblem_SolversCreateFinish(Problem);

  /* Start the creation of the problem solver equations */
  Solver=(CMISSSolverType)NULL;
  Err = CMISSSolver_Initialise(&Solver);
  Err = CMISSSolverEquations_Initialise(&SolverEquations);
  Err = CMISSProblem_SolverEquationsCreateStart(Problem);
  /* Get the solve equations */
  Err = CMISSProblem_SolverGet(Problem,1,ControlLoopIdentifier,1,Solver);
  Err = CMISSSolver_SolverEquationsGet(Solver,SolverEquations);
  /* Set the solver equations sparsity */
  Err = CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES);
  /* Err = CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_FULL_MATRICES);  */
  /* Add in the equations set */
  Err = CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,&EquationsSetIndex);
  /* Finish the creation of the problem solver equations */
  Err = CMISSProblem_SolverEquationsCreateFinish(Problem);

  /* Start the creation of the equations set boundary conditions */
  Err = CMISSBoundaryConditions_Initialise(&BoundaryConditions);
  Err = CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions);
  /* Set the first node to 0.0 and the last node to 1.0 */
  FirstNodeNumber=1;
  if(NUMBER_GLOBAL_Z_ELEMENTS==0)
    {
      LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1);
    }
  else
    {
      LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1);
    }
  Err = CMISSDecomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,&FirstNodeDomain);
  Err = CMISSDecomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,&LastNodeDomain);
  if(FirstNodeDomain==ComputationalNodeNumber)
    {
      Err = CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, \
        CMISS_BOUNDARY_CONDITION_FIXED,0.0);
    }
  if(LastNodeDomain==ComputationalNodeNumber)
    {
      Err = CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, \
        CMISS_BOUNDARY_CONDITION_FIXED,1.0);
    }
  /* Finish the creation of the equations set boundary conditions */
  Err = CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations);

  /* Solve the problem */
  Err = CMISSProblem_Solve(Problem);

  Err = CMISSFinalise();

  return Err;
}
