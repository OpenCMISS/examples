!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a finite elasticity equation using OpenCMISS calls.
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
!> Contributor(s): Jack Lee
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

!> \example FiniteElasticity/UniAxialExtension/src/UniAxialExtensionExample.f90
!! Example program to solve a finite elasticity equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/UniAxialExtension/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/UniAxialExtension/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM CYLINDERINFLATIONEXAMPLE

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

  !\todo: don't hard code, read in + default
  REAL(CMISSRP), PARAMETER :: INNER_PRESSURE=1.0_CMISSRP !Positive is compressive
  REAL(CMISSRP), PARAMETER :: OUTER_PRESSURE=0.0_CMISSRP !Positive is compressive
  REAL(CMISSRP), PARAMETER :: LAMBDA=1.1_CMISSRP
  REAL(CMISSRP), PARAMETER :: TSI=0.0_CMISSRP    !Not yet working. Leave at 0
  REAL(CMISSRP), PARAMETER :: INNER_RAD=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: OUTER_RAD=1.2_CMISSRP
  REAL(CMISSRP), PARAMETER :: HEIGHT=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: C1=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: C2=6.0_CMISSRP
  INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalXElements=2
  INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalYElements=8
  INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalZElements=2

!   !Standard test parameters (don't remove or change)
!   REAL(CMISSRP), PARAMETER :: INNER_PRESSURE=1.0_CMISSRP !Positive is compressive
!   REAL(CMISSRP), PARAMETER :: OUTER_PRESSURE=0.0_CMISSRP !Positive is compressive
!   REAL(CMISSRP), PARAMETER :: LAMBDA=1.1_CMISSRP
!   REAL(CMISSRP), PARAMETER :: TSI=0.0_CMISSRP    !Not yet used
!   REAL(CMISSRP), PARAMETER :: INNER_RAD=1.0_CMISSRP
!   REAL(CMISSRP), PARAMETER :: OUTER_RAD=1.2_CMISSRP
!   REAL(CMISSRP), PARAMETER :: HEIGHT=2.0_CMISSRP
!   REAL(CMISSRP), PARAMETER :: C1=2.0_CMISSRP
!   REAL(CMISSRP), PARAMETER :: C2=6.0_CMISSRP
!   INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalXElements=2 !\todo: don't hardcode?
!   INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalYElements=8
!   INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalZElements=2
!  Increment loop of 2

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponents=2

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariables=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponents=4

  INTEGER(CMISSIntg), PARAMETER :: FieldAnalyticUserNumber=6

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program types


  !Program variables
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  !CMISS variables

  TYPE(cmfe_BasisType) :: QuadraticBasis, LinearBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,FibreField,MaterialField
  TYPE(cmfe_FieldType) :: DependentField,EquationsSetField,AnalyticField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  !TYPE(cmfe_NodesType) :: Nodes
  !TYPE(cmfe_MeshElementsType) :: QuadraticElements,LinearElements
  TYPE(cmfe_ControlLoopType) :: ControlLoop

  !Other variables
  INTEGER(CMISSIntg) :: NE,E  

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
  CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_FINITE_ELEMENT_CALCULATE"],Err)

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  write(*,*) "NumberOfDomains=",NumberOfComputationalNodes
  NumberOfDomains=NumberOfComputationalNodes !1

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
!   CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR) ! NOW A PARAMETER
!   CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR) ! NOW A PARAMETER
!   CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR) ! NOW A PARAMETER 
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  
  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_TypeSet(CoordinateSystem,CMFE_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL cmfe_CoordinateSystem_OriginSet(CoordinateSystem,[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP],Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the CS to the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange
  CALL cmfe_Basis_Initialise(LinearBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis, &
    & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)

  CALL cmfe_Basis_Initialise(QuadraticBasis,Err)
  CALL cmfe_Basis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL cmfe_Basis_InterpolationXiSet(QuadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err) !Enable 3D interpolation on faces
  CALL cmfe_Basis_CreateFinish(QuadraticBasis,Err)

  !Start the creation of a generated cylinder mesh
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up an cylinder mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_CYLINDER_MESH_TYPE,Err)
  !Set the quadratic basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,[QuadraticBasis,LinearBasis],Err)
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[INNER_RAD, OUTER_RAD, HEIGHT],Err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements],Err)
  
  !Finish the creation of generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL cmfe_RandomSeedsSet(0_CMISSIntg,Err) !To keep the automatic decomposition same each time
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Automatic decomposition
!   CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
!   CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  !Manual decomposition
  IF(NumberOfDomains>1) THEN
    CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_USER_DEFINED_TYPE,Err)
    !Set all elements but last one to first domain
    CALL cmfe_Mesh_NumberOfElementsGet(Mesh,NE,Err)
    do E=1,NE/2
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,E,0,Err)
    enddo
    do E=NE/2+1,NE
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,E,1,Err)
    enddo
    CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  ENDIF
  CALL cmfe_Decomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)  
  CALL cmfe_Field_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_ScalingTypeSet(GeometricField,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Create a fibre field and attach it to the geometric field  
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_ScalingTypeSet(FibreField,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !Create a material field and attach it to the geometric field  
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL cmfe_Field_TypeSet(MaterialField,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialField,Decomposition,Err)        
  CALL cmfe_Field_GeometricFieldSet(MaterialField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL cmfe_Field_ScalingTypeSet(MaterialField,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_CreateFinish(MaterialField,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,C1,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,C2,Err)

  !Create the equations_set
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field with 2 variables and 4 components (3 displacement, 1 pressure)
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL cmfe_Field_TypeSet(DependentField,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL cmfe_Field_DependentTypeSet(DependentField,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ScalingTypeSet(DependentField,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL cmfe_Field_CreateFinish(DependentField,Err)

  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set up analytic field
  CALL cmfe_Field_Initialise(AnalyticField,Err)
  CALL cmfe_EquationsSet_AnalyticCreateStart(EquationsSet,CMFE_EQUATIONS_SET_FINITE_ELASTICITY_CYLINDER, &
    & FieldAnalyticUserNumber,AnalyticField,Err)
  !Finish the equations set analytic field variables
  CALL cmfe_EquationsSet_AnalyticCreateFinish(EquationsSet,Err)

  !Set the analytic parameters
  CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSet,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_PIN_IDX,INNER_PRESSURE, &
    & Err)
  CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSet,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_POUT_IDX, &
    & OUTER_PRESSURE, &
    & Err)
  CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSet,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_LAMBDA_IDX,LAMBDA,Err)
  CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSet,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_TSI_IDX,TSI,Err)
  CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSet,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_RIN_IDX,INNER_RAD,Err)
  CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSet,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_ROUT_IDX,OUTER_RAD,Err)
  CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSet,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C1_IDX,C1,Err)
  CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSet,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C2_IDX,C2,Err)

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
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
    & -14.0_CMISSRP, &
    & Err)

  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_NO_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoop,1,Err)  ! this one sets the increment loop counter
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)
  
  !Create the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(Solver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)  !Slower
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(Solver,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(Solver,1.0E-7_CMISSRP,Err)
  CALL cmfe_Solver_NewtonSolutionToleranceSet(Solver,1.0E-7_CMISSRP,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)   
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

!! MANUAL BC ASSIGNMENT - REPLACED BY THE ANALYTIC BC ROUTINE BELOW
!   !Prescribe boundary conditions (absolute nodal parameters)
!   CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
!   CALL cmfe_EquationsSetBoundaryConditionsCreateStart(EquationsSet,BoundaryConditions,Err)
! 
!   !Get surfaces - will fix two nodes on bottom face, pressure conditions inside
!   CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_CYLINDER_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi,Err)
!   CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_CYLINDER_INNER_SURFACE,InnerSurfaceNodes,InnerNormalXi,Err)
!   CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_CYLINDER_OUTER_SURFACE,OuterSurfaceNodes,OuterNormalXi,Err)
! 
!   !Set all inner surface nodes to inner pressure
!   DO NN=1,SIZE(InnerSurfaceNodes,1)
!       CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,InnerSurfaceNodes(NN), &
!       & abs(InnerNormalXi),CMFE_BOUNDARY_CONDITION_PRESSURE_INCREMENTED,INNER_PRESSURE,Err)   ! INNER_PRESSURE
!     IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER PRESSURE TO NODE", InnerSurfaceNodes(NN)
!   ENDDO
! 
!   !Set all outer surface nodes to outer pressure
!   DO NN=1,SIZE(OuterSurfaceNodes,1)
!     CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,OuterSurfaceNodes(NN), &
!       & abs(OuterNormalXi),CMFE_BOUNDARY_CONDITION_PRESSURE_INCREMENTED,OUTER_PRESSURE,Err)
!     IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER PRESSURE TO NODE", OuterSurfaceNodes(NN)
!   ENDDO
! 
!   !Set all top nodes fixed in z plane at the set height
!   DO NN=1,SIZE(TopSurfaceNodes,1)
!     CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CMFE_FIELD_U_VARIABLE_TYPE,1,TopSurfaceNodes(NN), &
!       & 3,CMFE_BOUNDARY_CONDITION_FIXED,deformedHeight,Err)
!     IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING FIXED CONDITION TO NODE", TopSurfaceNodes(NN)
!   ENDDO
! 
!   !Set all bottom nodes fixed in z plane
!   DO NN=1,SIZE(BottomSurfaceNodes,1)
!     CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CMFE_FIELD_U_VARIABLE_TYPE,1,BottomSurfaceNodes(NN), &
!       & 3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
!     IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING FIXED CONDITION TO NODE", BottomSurfaceNodes(NN)
!   ENDDO
! 
!   !Set two nodes on the bottom surface to axial displacement only
!   X_FIXED=.FALSE.
!   Y_FIXED=.FALSE.
!   DO NN=1,SIZE(BottomSurfaceNodes,1)
!     IF (.NOT.X_FIXED) THEN
!       CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
!         & 1,BottomSurfaceNodes(NN),1,xValue,Err)
!       IF(abs(xValue)<1e-5_CMISSRP) THEN
!         !Constrain it in x direction
!         CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CMFE_FIELD_U_VARIABLE_TYPE,1,BottomSurfaceNodes(NN),1, &
!           & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
!         X_FIXED=.TRUE.
!         WRITE(*,*) "CyliderInflationExample: SUCCESSFULLY CONSTRAINED IN X DIRECTION NODE",BottomSurfaceNodes(NN)
!       ENDIF
!     ENDIF
!     IF(.NOT.Y_FIXED) THEN
!       CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
!         & 1,BottomSurfaceNodes(NN),2,yValue,Err)
!       IF(abs(yValue)<1e-5_CMISSRP) THEN
!         !Constrain it in y direction
!         CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CMFE_FIELD_U_VARIABLE_TYPE,1,BottomSurfaceNodes(NN),2, &
!           & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
!         Y_FIXED=.TRUE.
!         WRITE(*,*) "CyliderInflationExample: SUCCESSFULLY CONSTRAINED IN Y DIRECTION NODE",BottomSurfaceNodes(NN)
!       ENDIF
!     ENDIF
!     IF (X_FIXED.AND.Y_FIXED) EXIT
!   ENDDO
!   !Check
!   IF(.NOT.X_FIXED .OR. .NOT.Y_FIXED) THEN
!     Write(*,*) "Couldn't fix bottom surface. No node lies on x or y axis, try changing number of elements"// &
!       & " in theta coordinate"
!     STOP
!   ENDIF
! 
!   CALL cmfe_EquationsSetBoundaryConditionsCreateFinish(EquationsSet,Err)

!! END OF MANUAL BC ASSIGNMENT

  !Set the bc using the analytic solution routine
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(SolverEquations,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !Output Analytic analysis
  Call cmfe_AnalyticAnalysis_Output(DependentField,"outputs/CylinderInflation",Err)

  !Output solution  
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"outputs/CylinderInflation","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"outputs/CylinderInflation","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM CYLINDERINFLATIONEXAMPLE

