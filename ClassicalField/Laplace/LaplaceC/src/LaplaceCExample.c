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

int main() 
{
  CMISSBasisType Basis=(CMISSBasisType)NULL;
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
  
  ControlLoopIdentifier[1]=CMISSControlLoopNode;

  if(CMISSInitialise(&WorldCoordinateSystem,&WorldRegion) == CMISSNoError)
    {
      Err = CMISSErrorHandlingModeSet(CMISSTrapError);

      Err = CMISSComputationalNumberOfNodesGet(&NumberOfComputationalNodes);
      Err = CMISSComputationalNodeNumberGet(&ComputationalNodeNumber);

      /* Start the creation of a new RC coordinate system */
      Err = CMISSCoordinateSystemTypeInitialise(&CoordinateSystem);
      Err = CMISSCoordinateSystemCreateStart(COORDINATE_SYSTEM_USER_NUMBER,&CoordinateSystem);
      if(NUMBER_GLOBAL_Z_ELEMENTS==0) 
	{
	  /* Set the coordinate system to be 2D */
	  Err = CMISSCoordinateSystemDimensionSet(CoordinateSystem,2);
	}
      else
	{
	  /* Set the coordinate system to be 3D */
	  Err = CMISSCoordinateSystemDimensionSet(CoordinateSystem,3);
	}
      /* Finish the creation of the coordinate system */
      Err = CMISSCoordinateSystemCreateFinish(CoordinateSystem);

      /* Start the creation of the region */
      Err = CMISSRegionTypeInitialise(&Region);
      Err = CMISSRegionCreateStart(REGION_USER_NUMBER,WorldRegion,&Region);
      /* Set the regions coordinate system to the 2D RC coordinate system that we have created */
      Err = CMISSRegionCoordinateSystemSet(Region,CoordinateSystem);
      /* Finish the creation of the region */
      Err = CMISSRegionCreateFinish(Region);

      /* Start the creation of a basis (default is trilinear lagrange) */
      Err = CMISSBasisTypeInitialise(&Basis);
      Err = CMISSBasisCreateStart(BASIS_USER_NUMBER,&Basis);
      if(NUMBER_GLOBAL_Z_ELEMENTS==0)
	{
	  /* Set the basis to be a bilinear Lagrange basis */
	  Err = CMISSBasisNumberOfXiSet(Basis,2);
	}
      else
	{
	  /* Set the basis to be a trilinear Lagrange basis */
	  Err = CMISSBasisNumberOfXiSet(Basis,3);
	}
      /* Finish the creation of the basis */
      Err = CMISSBasisCreateFinish(Basis);
   
      /* Start the creation of a generated mesh in the region */
      Err = CMISSGeneratedMeshTypeInitialise(&GeneratedMesh);
      Err = CMISSGeneratedMeshCreateStart(GENERATED_MESH_USER_NUMBER,Region,&GeneratedMesh);
      /* Set up a regular x*y*z mesh */
      Err = CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType);
      /* Set the default basis */
      Err = CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis);   
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
      Err = CMISSGeneratedMeshExtentSet(GeneratedMesh,MAX_COORDINATES,MeshExtent);
      Err = CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,MAX_COORDINATES,NumberXiElements);
      /* Finish the creation of a generated mesh in the region */
      Err = CMISSMeshTypeInitialise(&Mesh);
      Err = CMISSGeneratedMeshCreateFinish(GeneratedMesh,MESH_USER_NUMBER,&Mesh);

      /* Create a decomposition */
      Err = CMISSDecompositionTypeInitialise(&Decomposition);
      Err = CMISSDecompositionCreateStart(DECOMPOSITION_USER_NUMBER,Mesh,&Decomposition);
      /* Set the decomposition to be a general decomposition with the specified number of domains */
      Err = CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType);
      Err = CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfComputationalNodes);
      /* Finish the decomposition */
      Err = CMISSDecompositionCreateFinish(Decomposition);
  
      /* Start to create a default (geometric) field on the region */
      Err = CMISSFieldTypeInitialise(&GeometricField);
      Err = CMISSFieldCreateStart(GEOMETRIC_FIELD_USER_NUMBER,Region,&GeometricField);
      /* Set the decomposition to use */
      Err = CMISSFieldMeshDecompositionSet(GeometricField,Decomposition);
      /* Set the domain to be used by the field components. */
      Err = CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,1);
      Err = CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,1);
      if(NUMBER_GLOBAL_Z_ELEMENTS!=0)
	{
	  Err = CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,1);
	}
      /* Finish creating the field */
      Err = CMISSFieldCreateFinish(GeometricField);

      /* Update the geometric field parameters */
      Err = CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh);
  
      /* Create the equations_set */
      Err = CMISSEquationsSetTypeInitialise(&EquationsSet);
      Err = CMISSEquationsSetCreateStart(EQUATIONS_SET_USER_NUMBER,Region,GeometricField,CMISSEquationsSetClassicalFieldClass,CMISSEquationsSetLaplaceEquationType,CMISSEquationsSetStandardLaplaceSubtype,EQUATIONS_SET_FIELD_USER_NUMBER,EquationsSetField,&EquationsSet);
      /* Set the equations set to be a standard Laplace problem */
      //Err = CMISSEquationsSetSpecificationSet(EquationsSet,CMISSEquationsSetClassicalFieldClass,CMISSEquationsSetLaplaceEquationType,CMISSEquationsSetStandardLaplaceSubtype);
      /* Finish creating the equations set */
      Err = CMISSEquationsSetCreateFinish(EquationsSet);
      
      /* Create the equations set dependent field variables */
      Err = CMISSFieldTypeInitialise(&DependentField);
      Err = CMISSEquationsSetDependentCreateStart(EquationsSet,DEPENDENT_FIELD_USER_NUMBER,&DependentField);
      /* Finish the equations set dependent field variables */
      Err = CMISSEquationsSetDependentCreateFinish(EquationsSet);

      /* Create the equations set equations */
      Err = CMISSEquationsTypeInitialise(&Equations);
      Err = CMISSEquationsSetEquationsCreateStart(EquationsSet,&Equations);
      /* Set the equations matrices sparsity type */
      Err = CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices);
      /* Set the equations set output */
      /* Err = CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput); */
      Err = CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput);
      /* Err = CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput); */
      /* Err = CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput); */
      /* Finish the equations set equations */
      Err = CMISSEquationsSetEquationsCreateFinish(EquationsSet);

      /* Start the creation of the equations set boundary conditions */
      Err = CMISSBoundaryConditionsTypeInitialise(&BoundaryConditions);
      Err = CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSet,&BoundaryConditions);
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
      Err = CMISSDecompositionNodeDomainGet(Decomposition,FirstNodeNumber,1,&FirstNodeDomain);
      Err = CMISSDecompositionNodeDomainGet(Decomposition,LastNodeNumber,1,&LastNodeDomain);
      if(FirstNodeDomain==ComputationalNodeNumber) 
	{
	  Err = CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,FirstNodeNumber,1,CMISSBoundaryConditionFixed,0.0);
	}
      if(LastNodeDomain==ComputationalNodeNumber)
	{
	  Err = CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,LastNodeNumber,1,CMISSBoundaryConditionFixed,1.0);
	}
      /* Finish the creation of the equations set boundary conditions */
      Err = CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSet);
  
      /* Start the creation of a problem. */
      Err = CMISSProblemTypeInitialise(&Problem);
      Err = CMISSProblemCreateStart(PROBLEM_USER_NUMBER,&Problem);
      /* Set the problem to be a standard Laplace problem */
      Err = CMISSProblemSpecificationSet(Problem,CMISSProblemClassicalFieldClass,CMISSProblemLaplaceEquationType,CMISSProblemStandardLaplaceSubtype);
      /* Finish the creation of a problem. */
      Err = CMISSProblemCreateFinish(Problem);

      /* Start the creation of the problem control loop */
      Err = CMISSProblemControlLoopCreateStart(Problem);
      /* Finish creating the problem control loop */
      Err = CMISSProblemControlLoopCreateFinish(Problem);
 
      /* Start the creation of the problem solvers */
      Err = CMISSSolverTypeInitialise(&Solver);
      Err = CMISSProblemSolversCreateStart(Problem);
      Err = CMISSProblemSolverGet(Problem,1,ControlLoopIdentifier,1,&Solver);
      /* Err = CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput); */
      /* Err = CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput); */
      /* Err = CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput); */
      /* Err = CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput); */
      Err = CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput);
      Err = CMISSSolverLinearTypeSet(Solver,CMISSSolverLinearDirectSolveType);
      Err = CMISSSolverLibraryTypeSet(Solver,CMISSSolverMUMPSLibrary);
      /* Finish the creation of the problem solver */
      Err = CMISSProblemSolversCreateFinish(Problem);

      /* Start the creation of the problem solver equations */
      Err = CMISSSolverTypeInitialise(&Solver);
      Err = CMISSSolverEquationsTypeInitialise(&SolverEquations);
      Err = CMISSProblemSolverEquationsCreateStart(Problem);
      /* Get the solve equations */
      Err = CMISSProblemSolverGet(Problem,1,ControlLoopIdentifier,1,&Solver);
      Err = CMISSSolverSolverEquationsGet(Solver,&SolverEquations);
      /* Set the solver equations sparsity */
      Err = CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices);
      /* Err = CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices);  */
      /* Add in the equations set */
      Err = CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,&EquationsSetIndex);
      /* Finish the creation of the problem solver equations */
      Err = CMISSProblemSolverEquationsCreateFinish(Problem);

      /* Solve the problem */
      Err = CMISSProblemSolve(Problem);

      Err = CMISSFinalise();
    }

  return Err;
}
