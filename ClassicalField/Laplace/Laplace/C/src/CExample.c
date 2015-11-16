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

#include "opencmiss/iron.h"

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
  if(Err != CMFE_NO_ERROR) { \
    if(Err == CMFE_ERROR_CONVERTING_POINTER) { \
      fprintf(stderr,"Error: %s: Error converting pointer.\n",(S)); \
    } \
    else if(Err == CMFE_POINTER_IS_NULL) { \
      fprintf(stderr,"Error: %s: Pointer is null.\n",(S)); \
    } \
    else if(Err == CMFE_POINTER_NOT_NULL) { \
      fprintf(stderr,"Error: %s: Pointer is not null.\n",(S)); \
    } \
    else if(Err == CMFE_COULD_NOT_ALLOCATE_POINTER) { \
      fprintf(stderr,"Error: %s: Could not allocate pointer.\n",(S)); \
    } \
    exit(Err); \
  }

int main()
{
  cmfe_BasisType Basis = (cmfe_BasisType)NULL;
  cmfe_BoundaryConditionsType BoundaryConditions=(cmfe_BoundaryConditionsType)NULL;
  cmfe_CoordinateSystemType CoordinateSystem=(cmfe_CoordinateSystemType)NULL,WorldCoordinateSystem=(cmfe_CoordinateSystemType)NULL;
  cmfe_DecompositionType Decomposition=(cmfe_DecompositionType)NULL;
  cmfe_EquationsType Equations=(cmfe_EquationsType)NULL;
  cmfe_EquationsSetType EquationsSet=(cmfe_EquationsSetType)NULL;
  cmfe_FieldType GeometricField=(cmfe_FieldType)NULL,DependentField=(cmfe_FieldType)NULL,EquationsSetField=(cmfe_FieldType)NULL;
  cmfe_FieldsType Fields=(cmfe_FieldsType)NULL;
  cmfe_GeneratedMeshType GeneratedMesh=(cmfe_GeneratedMeshType)NULL;
  cmfe_MeshType Mesh=(cmfe_MeshType)NULL;
  cmfe_ProblemType Problem=(cmfe_ProblemType)NULL;
  cmfe_RegionType Region=(cmfe_RegionType)NULL,WorldRegion=(cmfe_RegionType)NULL;
  cmfe_SolverType Solver=(cmfe_SolverType)NULL;
  cmfe_SolverEquationsType SolverEquations=(cmfe_SolverEquationsType)NULL;

  int NumberOfComputationalNodes,ComputationalNodeNumber;
  int EquationsSetIndex;
  int FirstNodeNumber,LastNodeNumber;
  int FirstNodeDomain,LastNodeDomain;

  int NumberXiElements[MAX_COORDINATES];
  int ControlLoopIdentifier[1];
  double MeshExtent[MAX_COORDINATES];

  int EquationsSetSpecification[3];
  int ProblemSpecification[3];

  int Err;

  ControlLoopIdentifier[0]=CMFE_CONTROL_LOOP_NODE;

  Err = cmfe_CoordinateSystem_Initialise(&WorldCoordinateSystem);
  CHECK_ERROR("Initialising world coordinate system");
  Err = cmfe_Region_Initialise(&WorldRegion);
  CHECK_ERROR("Initialising world region");
  Err = cmfe_Initialise(WorldCoordinateSystem,WorldRegion);
  CHECK_ERROR("Initialising OpenCMISS-Iron");
  Err = cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR);

  Err = cmfe_ComputationalNumberOfNodesGet(&NumberOfComputationalNodes);
  Err = cmfe_ComputationalNodeNumberGet(&ComputationalNodeNumber);

  /* Start the creation of a new RC coordinate system */
  Err = cmfe_CoordinateSystem_Initialise(&CoordinateSystem);
  Err = cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,CoordinateSystem);
  if(NUMBER_GLOBAL_Z_ELEMENTS==0)
    {
      /* Set the coordinate system to be 2D */
      Err = cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2);
    }
  else
    {
      /* Set the coordinate system to be 3D */
      Err = cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,3);
    }
  /* Finish the creation of the coordinate system */
  Err = cmfe_CoordinateSystem_CreateFinish(CoordinateSystem);

  /* Start the creation of the region */
  Err = cmfe_Region_Initialise(&Region);
  Err = cmfe_Region_CreateStart(REGION_USER_NUMBER,WorldRegion,Region);
  /* Set the regions coordinate system to the 2D RC coordinate system that we have created */
  Err = cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem);
  /* Finish the creation of the region */
  Err = cmfe_Region_CreateFinish(Region);

  /* Start the creation of a basis (default is trilinear lagrange) */
  Err = cmfe_Basis_Initialise(&Basis);
  Err = cmfe_Basis_CreateStart(BASIS_USER_NUMBER,Basis);
  if(NUMBER_GLOBAL_Z_ELEMENTS==0)
    {
      /* Set the basis to be a bilinear Lagrange basis */
      Err = cmfe_Basis_NumberOfXiSet(Basis,2);
    }
  else
    {
      /* Set the basis to be a trilinear Lagrange basis */
      Err = cmfe_Basis_NumberOfXiSet(Basis,3);
    }
  /* Finish the creation of the basis */
  Err = cmfe_Basis_CreateFinish(Basis);

  /* Start the creation of a generated mesh in the region */
  Err = cmfe_GeneratedMesh_Initialise(&GeneratedMesh);
  Err = cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,Region,GeneratedMesh);
  /* Set up a regular x*y*z mesh */
  Err = cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE);
  /* Set the default basis */
  Err = cmfe_GeneratedMesh_BasisSet(GeneratedMesh,1,&Basis);
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
  Err = cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,MAX_COORDINATES,MeshExtent);
  Err = cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,MAX_COORDINATES,NumberXiElements);
  /* Finish the creation of a generated mesh in the region */
  Err = cmfe_Mesh_Initialise(&Mesh);
  /* Finish the creation of a generated mesh in the region */
  Err = cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MESH_USER_NUMBER,Mesh);

  /* Create a decomposition */
  Err = cmfe_Decomposition_Initialise(&Decomposition);
  Err = cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,Mesh,Decomposition);
  /* Set the decomposition to be a general decomposition with the specified number of domains */
  Err = cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE);
  Err = cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes);
  /* Finish the decomposition */
  Err = cmfe_Decomposition_CreateFinish(Decomposition);

  /* Start to create a default (geometric) field on the region */
  Err = cmfe_Field_Initialise(&GeometricField);
  Err = cmfe_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,Region,GeometricField);
  /* Set the decomposition to use */
  Err = cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition);
  /* Set the domain to be used by the field components. */
  Err = cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1);
  Err = cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1);
  if(NUMBER_GLOBAL_Z_ELEMENTS!=0)
    {
      Err = cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1);
    }
  /* Finish creating the field */
  Err = cmfe_Field_CreateFinish(GeometricField);

  /* Update the geometric field parameters */
  Err = cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField);

  /* Create the equations_set */
  Err = cmfe_EquationsSet_Initialise(&EquationsSet);
  Err = cmfe_Field_Initialise(&EquationsSetField);
  EquationsSetSpecification[0] = CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS;
  EquationsSetSpecification[1] = CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE;
  EquationsSetSpecification[2] = CMFE_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE;
  Err = cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,Region,GeometricField, \
    3,EquationsSetSpecification,EQUATIONS_SET_FIELD_USER_NUMBER, \
    EquationsSetField,EquationsSet);
  CHECK_ERROR("Creating equations set");
  /* Finish creating the equations set */
  Err = cmfe_EquationsSet_CreateFinish(EquationsSet);

  /* Create the equations set dependent field variables */
  Err = cmfe_Field_Initialise(&DependentField);
  Err = cmfe_EquationsSet_DependentCreateStart(EquationsSet,DEPENDENT_FIELD_USER_NUMBER,DependentField);
  /* Finish the equations set dependent field variables */
  Err = cmfe_EquationsSet_DependentCreateFinish(EquationsSet);

  /* Create the equations set equations */
  Err = cmfe_Equations_Initialise(&Equations);
  Err = cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations);
  /* Set the equations matrices sparsity type */
  Err = cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES);
  /* Set the equations set output */
  /* Err = cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT); */
  Err = cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT);
  /* Err = cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT); */
  /* Err = cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT); */
  /* Finish the equations set equations */
  Err = cmfe_EquationsSet_EquationsCreateFinish(EquationsSet);

  /* Start the creation of a problem, setting the problem to be a standard Laplace problem. */
  Err = cmfe_Problem_Initialise(&Problem);
  ProblemSpecification[0] = CMFE_PROBLEM_CLASSICAL_FIELD_CLASS;
  ProblemSpecification[1] = CMFE_PROBLEM_LAPLACE_EQUATION_TYPE;
  ProblemSpecification[2] = CMFE_PROBLEM_STANDARD_LAPLACE_SUBTYPE;
  Err = cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,3,ProblemSpecification,Problem);
  /* Finish the creation of a problem. */
  Err = cmfe_Problem_CreateFinish(Problem);

  /* Start the creation of the problem control loop */
  Err = cmfe_Problem_ControlLoopCreateStart(Problem);
  /* Finish creating the problem control loop */
  Err = cmfe_Problem_ControlLoopCreateFinish(Problem);

  /* Start the creation of the problem solvers */
  Err = cmfe_Solver_Initialise(&Solver);
  Err = cmfe_Problem_SolversCreateStart(Problem);
  Err = cmfe_Problem_SolverGet(Problem,1,ControlLoopIdentifier,1,Solver);
  /* Err = cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT); */
  /* Err = cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT); */
  /* Err = cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT); */
  /* Err = cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT); */
  Err = cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT);
  Err = cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE);
  Err = cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_MUMPS_LIBRARY);
  /* Finish the creation of the problem solver */
  Err = cmfe_Problem_SolversCreateFinish(Problem);

  /* Start the creation of the problem solver equations */
  Solver=(cmfe_SolverType)NULL;
  Err = cmfe_Solver_Initialise(&Solver);
  Err = cmfe_SolverEquations_Initialise(&SolverEquations);
  Err = cmfe_Problem_SolverEquationsCreateStart(Problem);
  /* Get the solve equations */
  Err = cmfe_Problem_SolverGet(Problem,1,ControlLoopIdentifier,1,Solver);
  Err = cmfe_Solver_SolverEquationsGet(Solver,SolverEquations);
  /* Set the solver equations sparsity */
  Err = cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES);
  /* Err = cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES);  */
  /* Add in the equations set */
  Err = cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,&EquationsSetIndex);
  /* Finish the creation of the problem solver equations */
  Err = cmfe_Problem_SolverEquationsCreateFinish(Problem);

  /* Start the creation of the equations set boundary conditions */
  Err = cmfe_BoundaryConditions_Initialise(&BoundaryConditions);
  Err = cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions);
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
  Err = cmfe_Decomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,&FirstNodeDomain);
  Err = cmfe_Decomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,&LastNodeDomain);
  if(FirstNodeDomain==ComputationalNodeNumber)
    {
      Err = cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, \
        CMFE_BOUNDARY_CONDITION_FIXED,0.0);
    }
  if(LastNodeDomain==ComputationalNodeNumber)
    {
      Err = cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, \
        CMFE_BOUNDARY_CONDITION_FIXED,1.0);
    }
  /* Finish the creation of the equations set boundary conditions */
  Err = cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations);

  /* Solve the problem */
  Err = cmfe_Problem_Solve(Problem);

  Err = cmfe_Finalise();

  return Err;
}
