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
PROGRAM CUBOIDGENERICEXAMPLE

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
  CHARACTER(LEN=100),ALLOCATABLE :: ALL_ARGS(:)
  INTEGER(CMISSIntg) :: Basis1,Basis2
  REAL(CMISSRP) :: X_LENG,Y_LENG,Z_LENG   ! dimension of cuboid
  INTEGER(CMISSIntg)  ::   NumberGlobalXElements
  INTEGER(CMISSIntg)  ::   NumberGlobalYElements
  INTEGER(CMISSIntg)  ::   NumberGlobalZElements
  INTEGER(CMISSIntg)  ::   incrementSteps
  REAL(CMISSRP),ALLOCATABLE :: BC(:,:)
  REAL(CMISSRP), PARAMETER :: C1=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: C2=6.0_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: CubicBasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: DisplacementMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureMeshComponentNumber=2

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

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program types


  !Program variables
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  !CMISS variables

  TYPE(cmfe_BasisType) :: CubicBasis, QuadraticBasis, LinearBasis, Bases(2)
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,FibreField,MaterialField
  TYPE(cmfe_FieldType) :: DependentField,EquationsSetField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoop

  !Other variables
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face1Nodes(:),Face2Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face3Nodes(:),Face4Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face5Nodes(:),Face6Nodes(:)
  INTEGER(CMISSIntg),POINTER :: Nodes(:)
  INTEGER(CMISSIntg) :: FaceXi(6)
  INTEGER(CMISSIntg) :: I,VariableType

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

  !Parse command line options: basis_1, basis_2, num_els_x, .. z, x_leng,.. z_leng, BC
  CALL GET_ALL_ARGS(ALL_ARGS)
  CALL GET_ARGUMENT(ALL_ARGS,"basis_displacement",ARG_INT=Basis1)
  CALL GET_ARGUMENT(ALL_ARGS,"basis_pressure",ARG_INT=Basis2)
  CALL GET_ARGUMENT(ALL_ARGS,"NumElsX",ARG_INT=NumberGlobalXElements)
  CALL GET_ARGUMENT(ALL_ARGS,"NumElsY",ARG_INT=NumberGlobalYElements)
  CALL GET_ARGUMENT(ALL_ARGS,"NumElsZ",ARG_INT=NumberGlobalZElements)
  CALL GET_ARGUMENT(ALL_ARGS,"XLeng",ARG_DP=X_LENG)
  CALL GET_ARGUMENT(ALL_ARGS,"YLeng",ARG_DP=Y_LENG)
  CALL GET_ARGUMENT(ALL_ARGS,"ZLeng",ARG_DP=Z_LENG)
  CALL GET_ARGUMENT(ALL_ARGS,"IncrementSteps",ARG_INT=incrementSteps)
  CALL GET_BC(ALL_ARGS,BC)

  !Set all diganostic levels on for testing
  CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_FINITE_ELEMENT_CALCULATE"],Err)
  !CALL cmfe_DiagnosticsSetOff(Err)

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  write(*,*) "NumberOfDomains=",NumberOfComputationalNodes
  NumberOfDomains=NumberOfComputationalNodes

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR) ! NOW A PARAMETER
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR) ! NOW A PARAMETER
  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR) ! NOW A PARAMETER 
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

  !Define basis functions
  CALL cmfe_Basis_Initialise(LinearBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis, &
    & [CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err)
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)

  CALL cmfe_Basis_Initialise(QuadraticBasis,Err)
  CALL cmfe_Basis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL cmfe_Basis_InterpolationXiSet(QuadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & [CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err)
  CALL cmfe_Basis_CreateFinish(QuadraticBasis,Err)

  CALL cmfe_Basis_Initialise(CubicBasis,Err)
  CALL cmfe_Basis_CreateStart(CubicBasisUserNumber,CubicBasis,Err)
  CALL cmfe_Basis_InterpolationXiSet(CubicBasis,[CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(CubicBasis, &
    & [CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(CubicBasis,.true.,Err) !Enable 3D interpolation on faces
  CALL cmfe_Basis_CreateFinish(CubicBasis,Err)

  SELECT CASE(Basis1)
  CASE(2)
    Bases(1)=QuadraticBasis
  CASE(3)
    Bases(1)=CubicBasis
  CASE DEFAULT
    WRITE(*,*) "Invalid basis type for displacement"
    STOP
  END SELECT

  SELECT CASE(Basis2)
  CASE(1)
    Bases(2)=LinearBasis
  CASE(2)
    Bases(2)=QuadraticBasis
  CASE DEFAULT
    WRITE(*,*) "Invalid basis type for pressure"
    STOP
  END SELECT

  !Start the creation of a generated cylinder mesh
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up an regular (cuboid) mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the quadratic basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Bases,Err)
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[X_LENG,Y_LENG,Z_LENG],Err)
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
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  !Manual decomposition
!   IF(NumberOfDomains>1) THEN
!     CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_USER_DEFINED_TYPE,Err)
!     !Set all elements but last one to first domain
!     CALL cmfe_Mesh_NumberOfElementsGet(Mesh,NE,Err)
!     do E=1,NE/2
!       CALL cmfe_Decomposition_ElementDomainSet(Decomposition,E,0,Err)
!     enddo
!     do E=NE/2+1,NE
!       CALL cmfe_Decomposition_ElementDomainSet(Decomposition,E,1,Err)
!     enddo
!     CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
!   ENDIF
  CALL cmfe_Decomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)  
  CALL cmfe_Field_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,DisplacementMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,DisplacementMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,DisplacementMeshComponentNumber,Err)
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
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,PressureMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,PressureMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,PressureMeshComponentNumber,Err)
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
  CALL cmfe_Field_TypeSet(DependentField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL cmfe_Field_DependentTypeSet(DependentField,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,DisplacementMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,2,DisplacementMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,3,DisplacementMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,PressureMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,DisplacementMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,DisplacementMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,DisplacementMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,PressureMeshComponentNumber,Err)
  CALL cmfe_Field_ScalingTypeSet(DependentField,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_CreateFinish(DependentField,Err)

  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

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
  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoop,incrementSteps,Err)  ! this one sets the increment loop counter
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)
  
  !Create the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(Solver,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
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

  !BC Assignment
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Get surfaces - will fix two nodes on bottom face, pressure conditions inside
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,Face1Nodes,FaceXi(1),Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_BACK_SURFACE,Face2Nodes,FaceXi(2),Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,Face3Nodes,FaceXi(3),Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,Face4Nodes,FaceXi(4),Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_TOP_SURFACE,Face5Nodes,FaceXi(5),Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,Face6Nodes,FaceXi(6),Err)

  DO I=1,SIZE(BC,1)
    IF(INT(BC(I,3))==CMFE_BOUNDARY_CONDITION_PRESSURE.OR.INT(BC(I,3))==CMFE_BOUNDARY_CONDITION_PRESSURE_INCREMENTED) THEN
      VariableType=CMFE_FIELD_DELUDELN_VARIABLE_TYPE
    ELSE
      VariableType=CMFE_FIELD_U_VARIABLE_TYPE
    ENDIF
    SELECT CASE (INT(BC(I,1)))
    CASE (1)
      Nodes=>Face1Nodes
    CASE (2)
      Nodes=>Face2Nodes
    CASE (3)
      Nodes=>Face3Nodes
    CASE (4)
      Nodes=>Face4Nodes
    CASE (5)
      Nodes=>Face5Nodes
    CASE (6)
      Nodes=>Face6Nodes
    END SELECT
    CALL SET_BC(Decomposition,DependentField,GeometricField,BoundaryConditions,VariableType,INT(BC(I,3)),Nodes,INT(BC(I,2)), &
      & BC(I,4),ComputationalNodeNumber)
  ENDDO

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !Output solution  
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"outputs/CuboidGeneric","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"outputs/CuboidGeneric","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  CALL cmfe_Finalise(Err)
  IF(ALLOCATED(BC)) DEALLOCATE(BC)
  IF(ALLOCATED(Face1Nodes)) DEALLOCATE(Face1Nodes)
  IF(ALLOCATED(Face2Nodes)) DEALLOCATE(Face2Nodes)
  IF(ALLOCATED(Face3Nodes)) DEALLOCATE(Face3Nodes)
  IF(ALLOCATED(Face4Nodes)) DEALLOCATE(Face4Nodes)
  IF(ALLOCATED(Face5Nodes)) DEALLOCATE(Face5Nodes)
  IF(ALLOCATED(Face6Nodes)) DEALLOCATE(Face6Nodes)
  IF(ALLOCATED(ALL_ARGS)) DEALLOCATE(ALL_ARGS)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

  CONTAINS

  !-------------------------------------------------------------------
  !> Converts a string to upper case
  SUBROUTINE UPPER_CASE(LWORD,UWORD)
    CHARACTER(LEN=*),INTENT(IN) :: LWORD
    CHARACTER(LEN=LEN(LWORD)),INTENT(OUT) :: UWORD
    INTEGER(CMISSIntg) ::I,IC,NLEN

    NLEN = LEN_TRIM(LWORD)
    UWORD=LWORD
    DO I=1,NLEN
      IC = ICHAR(LWORD(I:I))
      IF (IC >= 97 .and. IC <= 122) UWORD(I:I) = CHAR(IC-32)
    ENDDO
  END SUBROUTINE UPPER_CASE

  !-------------------------------------------------------------------
  !> Just grabs the whole input file
  SUBROUTINE GET_ALL_ARGS(ALL_ARGS)
    CHARACTER(LEN=100),ALLOCATABLE,INTENT(INOUT) :: ALL_ARGS(:)
    !Local variables
    CHARACTER(LEN=100) :: ARGS(100)
    INTEGER(CMISSIntg) :: I,REASON
    CHARACTER(LEN=100) :: WORD

    ARGS=""
    I=0
    DO
      I=I+1
      IF(I>100) THEN
        WRITE(*,*) "Increase ARGS size"
        STOP
      ENDIF
      READ(5,'(A)',IOSTAT=Reason) WORD
      IF(Reason>0)  THEN
        WRITE(*,*) "Error reading options file"
        STOP
      ELSEIF(Reason<0) THEN
        EXIT
      ELSE
        CALL UPPER_CASE(WORD,ARGS(I))
      ENDIF
    ENDDO

    ALLOCATE(ALL_ARGS(I))
    ALL_ARGS=ARGS(1:I)

  END SUBROUTINE GET_ALL_ARGS

  !-------------------------------------------------------------------
  !> Returns the argument of the requested type as a varying string
  SUBROUTINE GET_ARGUMENT(ALL_ARGS,ARG_TYPE,ARG,ARG_INT,ARG_DP)
    CHARACTER(LEN=*), INTENT(IN) :: ALL_ARGS(:)
    CHARACTER(LEN=*), INTENT(IN) :: ARG_TYPE
    CHARACTER(LEN=256), INTENT(OUT), OPTIONAL :: ARG
    INTEGER(CMISSIntg), INTENT(OUT), OPTIONAL :: ARG_INT
    REAL(CMISSRP), INTENT(OUT), OPTIONAL :: ARG_DP
    !Local variables
    CHARACTER(LEN=LEN_TRIM(ARG_TYPE)) :: ARG_TYPE_UPPER
    CHARACTER(LEN=256) :: WORD_UPPER,ARGOUT
    INTEGER(CMISSIntg) :: NARGS,I,LENG,WORD_LENG

    NARGS=SIZE(ALL_ARGS,1)
    LENG=LEN_TRIM(ARG_TYPE) !STRING TO LOOK FOR
    CALL UPPER_CASE(ARG_TYPE,ARG_TYPE_UPPER)

    DO I=1,NARGS
      WORD_UPPER=ALL_ARGS(I)
      WORD_LENG=LEN_TRIM(WORD_UPPER)
      IF(WORD_UPPER(1:1+LENG)=="-"//TRIM(ARG_TYPE_UPPER)) THEN
        IF(WORD_UPPER(2+LENG:2+LENG)=="=") THEN
          ! USING = AS DELIMITER
          ARGOUT=WORD_UPPER(3+LENG:WORD_LENG)
          EXIT
        ELSE
          ! USING A SPACE AS DELIMITER
          ARGOUT=TRIM(WORD_UPPER(3+LENG:WORD_LENG))
          EXIT
        ENDIF
      ENDIF
    ENDDO

    IF(PRESENT(ARG)) CALL UPPER_CASE(ARGOUT,ARG)
    IF(PRESENT(ARG_INT)) READ(ARGOUT,*) ARG_INT
    IF(PRESENT(ARG_DP)) READ(ARGOUT,*) ARG_DP
  END SUBROUTINE GET_ARGUMENT

  !-------------------------------------------------------------------
  !> Sets a given boundary condition to a group of nodes
  SUBROUTINE SET_BC(Decomposition,DependentField,GeometricField,BoundaryConditions,VariableType,BCType,Nodes,Component,Value, &
    & ComputationalNodeNumber)
    TYPE(cmfe_DecompositionType),INTENT(IN) :: Decomposition
    TYPE(cmfe_FieldType),INTENT(IN) :: DependentField
    TYPE(cmfe_FieldType),INTENT(IN) :: GeometricField    
    TYPE(cmfe_BoundaryConditionsType),INTENT(INOUT) :: BoundaryConditions
    INTEGER(CMISSIntg),INTENT(IN) :: VariableType
    INTEGER(CMISSIntg),INTENT(IN) :: BCType
    INTEGER(CMISSIntg),INTENT(IN) :: Nodes(:)
    INTEGER(CMISSIntg),INTENT(IN) :: Component
    REAL(CMISSRP),INTENT(IN) :: Value
    INTEGER(CMISSIntg),INTENT(IN) :: ComputationalNodeNumber
    !Local variables
    INTEGER(CMISSIntg) :: i,node,NodeDomain,BCType2
    REAL(CMISSRP) :: coord

    DO I=1,SIZE(Nodes)
      node=Nodes(I)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(BCType<0) THEN
          !Displacement type condition - get the undeformed geometric value first
          CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node, &
            & Component,coord,Err)
          coord=coord+Value
          IF(BCType==-777) BCType2=CMFE_BOUNDARY_CONDITION_FIXED
          IF(BCType==-888) BCType2=CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,VariableType,1,1,node,Component,BCType2,coord,Err)
        ELSE
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,VariableType,1,1,node,Component,BCType,Value,Err)
        ENDIF
      ENDIF
    ENDDO
    
  END SUBROUTINE SET_BC

  !-------------------------------------------------------------------
  !> Grabs all specified BC from the command line
  SUBROUTINE GET_BC(ALL_ARGS,BC)
    CHARACTER(LEN=*), INTENT(IN) :: ALL_ARGS(:)
    REAL(CMISSRP),ALLOCATABLE,INTENT(OUT) :: BC(:,:)
    !Local variable
    REAL(CMISSRP) :: BCtemp(100,4)
    INTEGER(CMISSIntg) :: NARG,N,I,J,POS1,POS2
    CHARACTER(LEN=100) :: WORD_UPPER
    CHARACTER(LEN=100) :: WORDS(10)
    
    N=0

    !Find out how many BCs there are
    NARG=SIZE(ALL_ARGS,1)
    DO I=1,NARG
      WORD_UPPER=ALL_ARGS(I)
      IF(WORD_UPPER(1:3)=="-BC") THEN
        N=N+1
        !Split the string into words
        POS1=1
        J=0
        DO
          POS2 = INDEX(TRIM(WORD_UPPER(POS1:)), " ")
          IF(POS2 == 0) THEN
            J=J+1
            WORDS(J) = TRIM(WORD_UPPER(POS1:))
            EXIT
          ENDIF
          J=J+1
          WORDS(J) = WORD_UPPER(POS1:POS1+POS2-2)
          POS1=POS2+POS1
        ENDDO
        ! surface type
        SELECT CASE (TRIM(WORDS(2)))
        CASE ("FRONT")
          BCtemp(N,1)=1
        CASE ("BACK")
          BCtemp(N,1)=2
        CASE ("RIGHT")
          BCtemp(N,1)=3
        CASE ("LEFT")
          BCtemp(N,1)=4
        CASE ("TOP")
          BCtemp(N,1)=5
        CASE ("BOTTOM")
          BCtemp(N,1)=6
        CASE DEFAULT
          WRITE(*,*) "Invalid type of surface in bc specification"
          STOP
        END SELECT
        ! component
        READ(WORDS(3),*) BCtemp(N,2)
        ! bc type
        WORD_UPPER=WORDS(4)
        SELECT CASE (TRIM(WORD_UPPER))
        CASE ("DISPLACEMENT")
          BCtemp(N,3)=-777
        CASE ("DISPLACEMENTINCREMENTED")
          BCtemp(N,3)=-888
        CASE ("FIXED")
          BCtemp(N,3)=CMFE_BOUNDARY_CONDITION_FIXED
        CASE ("FIXEDINCREMENTED")
          BCtemp(N,3)=CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED
        CASE ("PRESSURE")
          BCtemp(N,3)=CMFE_BOUNDARY_CONDITION_PRESSURE
        CASE ("PRESSUREINCREMENTED")
          BCtemp(N,3)=CMFE_BOUNDARY_CONDITION_PRESSURE_INCREMENTED
        CASE DEFAULT
          WRITE(*,*) "Invalid type of bc in bc specification"
          STOP
        END SELECT
        ! value
        WORD_UPPER=WORDS(5)
        READ(WORD_UPPER,*) BCtemp(N,4)
      ENDIF
    ENDDO

    ALLOCATE(BC(N,4)) ! face_no,component,bctype,value
    BC=BCtemp(1:N,:)
    
  END SUBROUTINE GET_BC


END PROGRAM CUBOIDGENERICEXAMPLE

