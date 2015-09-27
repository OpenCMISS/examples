!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a generalised Laplace equation using OpenCMISS calls.
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

!> \example ClassicalField/Laplace/GeneralisedLaplace/Fortran/src/FortranExample.f90
!! Example program to solve a Generalised Laplace equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/Laplace/GeneralisedLaplace/history.html
!!
!<

!> Main program
PROGRAM GENERALISEDLAPLACEEXAMPLE

  USE OPENCMISS
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

  REAL(CMISSDP), PARAMETER :: WIDTH=4.0_CMISSDP   ! x-direction
  REAL(CMISSDP), PARAMETER :: HEIGHT=2.0_CMISSDP  ! y-direction
  REAL(CMISSDP), PARAMETER :: LENGTH=0.0_CMISSDP  ! z-direction
  
  REAL(CMISSDP), PARAMETER :: PI=4.0_CMISSDP*DATAN(1.0_CMISSDP)

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=12
  
  INTEGER(CMISSIntg), PARAMETER :: FibreFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: FibreFieldNumberOfVariables=1
  
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
  
  INTEGER(CMISSIntg), PARAMETER :: DerivativeUserNumber=1
 
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
    & INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI,node_idx,component_idx,TotalNumberOfNodes
  INTEGER(CMISSIntg) :: FibreFieldNumberOfComponents  
  REAL(CMISSDP) :: SIGMA11,SIGMA22,SIGMA12,X,Y,VALUE
  REAL(CMISSDP) :: ANALYTICAL_SOL_X_Y,CALCULATED_SOLUTION_X_Y,NODAL_ERROR,TOTAL_ERROR_SQUARED,DX,DY,ERROR_L2
  REAL(CMISSDP) :: FibreFieldAngle(3)
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT,Filename
  
  INTEGER(CMISSIntg),ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: BackSurfaceNodes(:)
  INTEGER(CMISSIntg) :: FrontNormalXi,LeftNormalXi,RightNormalXi,BackNormalXi

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,EquationsSetField,DependentField,MaterialsField,FibreField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber,NodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain,NodeDomain
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

  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS >= 4) THEN
    !If we have enough arguments then use the first four for setting up the problem. The subsequent arguments may be used to
    !pass flags to, say, PETSc.
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_X_ELEMENTS
    IF(NUMBER_GLOBAL_X_ELEMENTS<=0) CALL HANDLE_ERROR("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Y_ELEMENTS
    IF(NUMBER_GLOBAL_Y_ELEMENTS<=0) CALL HANDLE_ERROR("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Z_ELEMENTS
    IF(NUMBER_GLOBAL_Z_ELEMENTS<0) CALL HANDLE_ERROR("Invalid number of Z elements.")
    CALL GET_COMMAND_ARGUMENT(4,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 4.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) INTERPOLATION_TYPE
    IF(INTERPOLATION_TYPE<=0) CALL HANDLE_ERROR("Invalid Interpolation specification.")
  ELSE
    !If there are not enough arguments default the problem specification 
    NUMBER_GLOBAL_X_ELEMENTS=16
    NUMBER_GLOBAL_Y_ELEMENTS=8
    NUMBER_GLOBAL_Z_ELEMENTS=0
!    INTERPOLATION_TYPE=1
    
!    INTERPOLATION_TYPE=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
!    INTERPOLATION_TYPE=CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
    INTERPOLATION_TYPE=CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION    
    
  ENDIF
  
  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  CALL CMISSRandomSeedsSet(9999,Err)
  
  CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"],Err)
  
  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  WRITE(Filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "GeneralisedLaplace",NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS,INTERPOLATION_TYPE
    
  CALL CMISSOutputSetOn(Filename,Err)

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_LabelSet(Region,"GeneralisedLaplaceRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1,2,3,4)
    CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1)
    NUMBER_OF_GAUSS_XI=2
  CASE(2)
    NUMBER_OF_GAUSS_XI=3
  CASE(3,4)
    NUMBER_OF_GAUSS_XI=4
  CASE DEFAULT
    NUMBER_OF_GAUSS_XI=0 !Don't set number of Gauss points for tri/tet
  END SELECT
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bi-interpolation basis
    CALL CMISSBasis_NumberOfXiSet(Basis,2,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ELSE
    !Set the basis to be a tri-interpolation basis
    CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(Basis,Err)
   
  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)
 
  !Destory the mesh now that we have decomposed it
!  CALL CMISSMesh_Destroy(Mesh,Err)
 
  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)
  
  !Start to create a fibre field and attach it to the geometric field
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FibreFieldUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  !Set the decomposition to use  
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreField,FibreFieldNumberOfVariables,Err)
  
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    FibreFieldNumberOfComponents=2
  ELSE
    FibreFieldNumberOfComponents=3
  ENDIF    
  
  CALL CMISSField_NumberOfComponentsSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,FibreFieldNumberOfComponents,Err) 
  CALL CMISSField_VariableLabelSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field  
  CALL CMISSField_CreateFinish(FibreField,Err)

  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSRegion_NodesGet(Region,Nodes,Err)
  CALL CMISSNodes_NumberOfNodesGet(Nodes,TotalNumberOfNodes,Err)

  !Rotation Angles (in radiant!!)
  ! in 2D an entry in Angle(1) means rotated x-axis, 
  !          entry in Angle(2) doesn't make sense, as rotates out of surface ...
  ! in 3D an entry in Angle(1) means rotated around z-axis, entry in Angle(2) means rotated around y-axis
  !          entry in Angle(3) means rotated around x-axis => no change
  ! 45° equivalent to pi/4, 90° equivalent to pi/2
  
!  FibreFieldAngle=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)

  FibreFieldAngle=(/PI/4,0.0_CMISSDP,0.0_CMISSDP/)
!  FibreFieldAngle=(/PI/2,0.0_CMISSDP,0.0_CMISSDP/)

  DO node_idx=1,TotalNumberOfNodes
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      DO component_idx=1,FibreFieldNumberOfComponents
        CALL CMISSField_ParameterSetUpdateNode(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
          & DerivativeUserNumber,node_idx,component_idx,FibreFieldAngle(component_idx),Err)
      ENDDO
    ENDIF
  ENDDO

  !WRITE(*,*) TotalNumberOfNodes
  !WRITE(*,*) FibreFieldAngle(1)
  
  !Start to create material field on the region
  CALL CMISSField_Initialise(MaterialsField,Err)
  CALL CMISSField_CreateStart(MaterialsFieldUserNumber,Region,MaterialsField,Err)
  CALL CMISSField_TypeSet(MaterialsField,CMISS_FIELD_MATERIAL_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(MaterialsField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(MaterialsField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialsField,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    ! symmetric 2x2 tensor => 3 different entries
    ! 1 - 11, 2 - 22, 3 - 12=21
    CALL CMISSField_NumberOfComponentsSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,3,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,3,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  ELSE
    ! symmetric 3x3 tensor => 6 different entries
    ! 1 - 11, 2 - 22, 3 - 33, 4 - 12=21, 5 - 23=32, 6 - 13=31
    CALL CMISSField_NumberOfComponentsSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,6,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,3,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,5,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,6,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  ENDIF
  !default is CMISS_FIELD_NODE_BASED_INTERPOLATION
  CALL CMISSField_VariableLabelSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"Material",Err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(MaterialsField,Err)

  !Set material parameters
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 1,0.5_CMISSDP,Err)  ! 11
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 2,0.1_CMISSDP,Err)  ! 22
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 3,0.0_CMISSDP,Err)  ! 12=21
  ELSE
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 1,0.5_CMISSDP,Err)  ! 11
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 2,0.1_CMISSDP,Err)  ! 22
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 3,0.1_CMISSDP,Err)  ! 33
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 4,0.0_CMISSDP,Err)  ! 12=21
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 5,0.0_CMISSDP,Err)  ! 23=32
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 6,0.0_CMISSDP,Err)  !13=31
  ENDIF   
  
  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
  
  !Create the Generalised Laplace Equations set
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMISS_EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Set the DOFs to be contiguous across components
  CALL CMISSField_DOFOrderTypeSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  CALL CMISSField_DOFOrderTypeSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Initialise the field with an initial guess
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.5_CMISSDP, &
    & Err)

  !Create the equations set material field variables
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)
  
  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a generalised Laplace problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS,CMISS_PROBLEM_LAPLACE_EQUATION_TYPE, &
    & CMISS_PROBLEM_GENERALISED_LAPLACE_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)
 
  !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  
!  CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
!  CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-12_CMISSDP,Err)
!  CALL CMISSSolver_LinearIterativeRelativeToleranceSet(Solver,1.0E-12_CMISSDP,Err)
!  CALL CMISSSolver_LinearIterativeTypeSet(Solver,CMISS_SOLVER_ITERATIVE_BiCGSTAB,Err)
!!  CALL CMISSSolver_LinearIterativeTypeSet(Solver,CMISS_SOLVER_ITERATIVE_CONJUGATE_GRADIENT,Err)
!!  !CALL CMISSSolver_LinearIterativeTypeSet(Solver,CMISS_SOLVER_ITERATIVE_GMRES,Err)
!  CALL CMISSSolver_LinearIterativeMaximumIterationsSet(Solver,1000000000,Err)
  
  CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  
  !CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  !CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_LAPACK_LIBRARY,Err)
  !CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_SUPERLU_LIBRARY,Err)
  !CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_PASTIX_LIBRARY,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Start the creation of the equations set boundary conditions
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  
!  !Set the first node to 0.0 and the last node to 1.0
!  FirstNodeNumber=1
!  CALL CMISSNodes_Initialise(Nodes,Err)
!  CALL CMISSRegion_NodesGet(Region,Nodes,Err)
!  CALL CMISSNodes_NumberOfNodesGet(Nodes,LastNodeNumber,Err)
!  CALL CMISSDecomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
!  CALL CMISSDecomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
!  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
!      & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
!  ENDIF
!  IF(LastNodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
!      & CMISS_BOUNDARY_CONDITION_FIXED,1.0_CMISSDP,Err)
!  ENDIF
  
  !Set the boundaries to analytical value of phi (-> see function below, ref.: phd-thesis Chris Bradley)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_BACK_SURFACE,BackSurfaceNodes,BackNormalXi,Err)
  
!  !only valid as sigma for no rotation, as effective conductivity needed for analytical example
!  CALL CMISSField_ParameterSetGetConstant(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,SIGMA11,Err)
!  CALL CMISSField_ParameterSetGetConstant(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,SIGMA22,Err)
!  CALL CMISSField_ParameterSetGetConstant(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,SIGMA12,Err)
  
  ! front
  DO node_idx=1,SIZE(FrontSurfaceNodes,1)
    NodeNumber=FrontSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      !WRITE(*,*) NodeNumber
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,&
        & X,Err)
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,2,&
        & Y,Err)
!      VALUE=PHI(X,Y,SIGMA11,SIGMA22,SIGMA12)
      VALUE=PHI(X,Y,0.3_CMISSDP,0.3_CMISSDP,0.2_CMISSDP)
!      VALUE=PHI(X,Y,2.0_CMISSDP,1.0_CMISSDP,0.0_CMISSDP)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  ENDDO
  ! left
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      !WRITE(*,*) NodeNumber    
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,&
        & X,Err)
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,2,&
        & Y,Err)
!      VALUE=PHI(X,Y,SIGMA11,SIGMA22,SIGMA12)
      VALUE=PHI(X,Y,0.3_CMISSDP,0.3_CMISSDP,0.2_CMISSDP) 
!      VALUE=PHI(X,Y,2.0_CMISSDP,1.0_CMISSDP,0.0_CMISSDP)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  ENDDO
  ! right
  DO node_idx=1,SIZE(RightSurfaceNodes,1)
    NodeNumber=RightSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      !WRITE(*,*) NodeNumber
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,&
        & X,Err)
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,2,&
        & Y,Err)
!      VALUE=PHI(X,Y,SIGMA11,SIGMA22,SIGMA12)
      VALUE=PHI(X,Y,0.3_CMISSDP,0.3_CMISSDP,0.2_CMISSDP)
!      VALUE=PHI(X,Y,2.0_CMISSDP,1.0_CMISSDP,0.0_CMISSDP)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  ENDDO
  ! back
  DO node_idx=1,SIZE(BackSurfaceNodes,1)
    NodeNumber=BackSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
       !WRITE(*,*) NodeNumber
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,&
        & X,Err)
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,2,&
        & Y,Err)
!      VALUE=PHI(X,Y,SIGMA11,SIGMA22,SIGMA12)
      VALUE=PHI(X,Y,0.3_CMISSDP,0.3_CMISSDP,0.2_CMISSDP)
!      VALUE=PHI(X,Y,2.0_CMISSDP,1.0_CMISSDP,0.0_CMISSDP)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  ENDDO
  
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Export results
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"GeneralisedLaplace","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"GeneralisedLaplace","FORTRAN",Err)
  
!--------------------------------------------------------------------------------------------------------------------------------  
  ! do the calculation of the L2-norm-error
  ! calculate analytical solution and error for each node
  TOTAL_ERROR_SQUARED=0.0_CMISSDP
  WRITE(*,*) "#        node         X-value           Y-value         analytical_solution     calculated_solution     nodal_error"
  
  DO node_idx=1,TotalNumberOfNodes
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,1,&
        & X,Err)
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,2,&
        & Y,Err)
      ANALYTICAL_SOL_X_Y=PHI(X,Y,0.3_CMISSDP,0.3_CMISSDP,0.2_CMISSDP)
      CALL CMISSField_ParameterSetGetNode(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,1,&
        & CALCULATED_SOLUTION_X_Y,Err)
      NODAL_ERROR=ABS(ANALYTICAL_SOL_X_Y-CALCULATED_SOLUTION_X_Y)
      TOTAL_ERROR_SQUARED=TOTAL_ERROR_SQUARED+NODAL_ERROR**2
    ENDIF
    WRITE(*,*) node_idx, X, Y, ANALYTICAL_SOL_X_Y, CALCULATED_SOLUTION_X_Y, NODAL_ERROR
  ENDDO

  DX=WIDTH/NUMBER_GLOBAL_X_ELEMENTS
  DY=HEIGHT/NUMBER_GLOBAL_Y_ELEMENTS
  ERROR_L2=SQRT(DX*DY*TOTAL_ERROR_SQUARED)
  
  WRITE(*,*) "# L2NORM_ERROR: ",ERROR_L2
  WRITE(*,*) DX, DY

!--------------------------------------------------------------------------------------------------------------------------------    
  CALL CMISSFields_Finalise(Fields,Err)
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

!  WRITE(*,*) PHI(1.0_CMISSDP,1.0_CMISSDP,0.3_CMISSDP,0.3_CMISSDP,0.2_CMISSDP)
!  WRITE(*,*) PHI(2.0_CMISSDP,1.0_CMISSDP,0.3_CMISSDP,0.3_CMISSDP,0.2_CMISSDP)
!  WRITE(*,*) PHI(3.0_CMISSDP,1.0_CMISSDP,0.3_CMISSDP,0.3_CMISSDP,0.2_CMISSDP)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
!--------------------------------------------------------------------------------------------------------------------------------  
CONTAINS
  !------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
  !------------------------------------------------------------------------------------------------------------------------------
  ! analytical solution
  REAL(CMISSDP) FUNCTION PHI(X,Y,SIGMA11,SIGMA22,SIGMA12)
    REAL(CMISSDP), INTENT(IN) :: X,Y,SIGMA11,SIGMA22,SIGMA12
    PHI = 2.0_CMISSDP*EXP(X)*EXP(-SIGMA12/SIGMA22*Y)*COS(SQRT(SIGMA11*SIGMA22-SIGMA12*SIGMA12)/SIGMA22*Y)
    RETURN
  END FUNCTION PHI
  !------------------------------------------------------------------------------------------------------------------------------  
END PROGRAM GENERALISEDLAPLACEEXAMPLE
