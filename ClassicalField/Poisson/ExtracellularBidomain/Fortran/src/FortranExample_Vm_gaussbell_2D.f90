!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a special case of a generalised Poisson equation, 
!<  namely the extracellular bidomain equation, using OpenCMISS calls.
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

!> \example ClassicalField/Poisson/ExtracellularBidomain/Fortran/src/FortranExample.f90
!! Example program to solve a special case of a generalised Poisson equation, 
!! namely the extracellular bidomain equation, using OpenCMISS calls.
!! \htmlinclude ClassicalField/Poisson/ExtracellularBidomain/history.html
!!
!<

!> Main program
PROGRAM EXTRACELLULARBIDOMAINEXAMPLE

  USE OPENCMISS
  USE MPI
!  USE CONSTANTS   !for pi

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: WIDTH=8.0_CMISSDP   ! x-direction
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
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=14  
  
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
    & INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI,node_idx,component_idx,TotalNumberOfNodes,idx
  INTEGER(CMISSIntg) :: FibreFieldNumberOfComponents  
!  REAL(CMISSDP) :: SIGMA11,SIGMA22,SIGMA12,X,Y,VALUE
!  REAL(CMISSDP) :: ANALYTICAL_SOL_X_Y,CALCULATED_SOLUTION_X_Y,NODAL_ERROR,TOTAL_ERROR_SQUARED,DX,DY,DZ,ERROR_L2
  REAL(CMISSDP) :: FibreFieldAngle(3),Vm(33)
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT,Filename
  
!  INTEGER(CMISSIntg),ALLOCATABLE :: FrontSurfaceNodes(:)
!  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
!  INTEGER(CMISSIntg),ALLOCATABLE :: RightSurfaceNodes(:)
!  INTEGER(CMISSIntg),ALLOCATABLE :: BackSurfaceNodes(:)
!  INTEGER(CMISSIntg) :: FrontNormalXi,LeftNormalXi,RightNormalXi,BackNormalXi

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,EquationsSetField,DependentField,MaterialsField,FibreField,SourceField
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
    NUMBER_GLOBAL_X_ELEMENTS=32
    NUMBER_GLOBAL_Y_ELEMENTS=4
    NUMBER_GLOBAL_Z_ELEMENTS=0
    
    INTERPOLATION_TYPE=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
!    INTERPOLATION_TYPE=CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
!    INTERPOLATION_TYPE=CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION    
    
  ENDIF
  
  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  CALL CMISSRandomSeedsSet(9999,Err)
  
  CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"],Err)
  
  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  WRITE(Filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "ExtracellularBidomain",NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
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
  CALL CMISSRegion_LabelSet(Region,"ExtracellularBidomainRegion",Err)
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

!--------------------------------------------------------------------------------------------------------------------------------
! FIELDS
!--------------------------------------------------------------------------------------------------------------------------------
!
! GEOMETRIC 
 
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
  
!-------------------------------------------------------------------------
!
! FIBRE  
  
  !Start to create a fibre field and attach it to the geometric field
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FibreFieldUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  !Set the decomposition to use  
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreField,FibreFieldNumberOfVariables,Err)
  CALL CMISSField_VariableTypesSet(FibreField,[CMISS_FIELD_U_VARIABLE_TYPE],Err)  
  
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    FibreFieldNumberOfComponents=2  
  ELSE
    FibreFieldNumberOfComponents=3
  ENDIF    
  
  CALL CMISSField_NumberOfComponentsSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,FibreFieldNumberOfComponents,Err) 
  CALL CMISSField_VariableLabelSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)   ! as 2D = 1, after if
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
  
!  FibreFieldAngle=(/PI/4.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
  FibreFieldAngle=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)

  DO node_idx=1,TotalNumberOfNodes
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      DO component_idx=1,FibreFieldNumberOfComponents
        CALL CMISSField_ParameterSetUpdateNode(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
          & DerivativeUserNumber,node_idx,component_idx,FibreFieldAngle(component_idx),Err)
      ENDDO
    ENDIF
  ENDDO

!-------------------------------------------------------------------------
!
! MATERIAL  
  
  !Start to create a material field on the region
  CALL CMISSField_Initialise(MaterialsField,Err)
  CALL CMISSField_CreateStart(MaterialsFieldUserNumber,Region,MaterialsField,Err)
  CALL CMISSField_TypeSet(MaterialsField,CMISS_FIELD_MATERIAL_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(MaterialsField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(MaterialsField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialsField,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    ! need 2 tensors, sigma_i and sigma_e
    ! symmetric 2x2 tensor => 3 different entries
    ! 1 - 11, 2 - 22, 3 - 12=21
    ! first i then e
    CALL CMISSField_NumberOfComponentsSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,6,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,3,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,5,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,6,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)    
  ELSE
    ! need 2 tensors, sigma_i and sigma_e
    ! symmetric 3x3 tensor => 6 different entries
    ! 1 - 11, 2 - 22, 3 - 33, 4 - 12=21, 5 - 23=32, 6 - 13=31
    ! first i then e 
    CALL CMISSField_NumberOfComponentsSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,12,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,3,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,5,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,6,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,7,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,8,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,9,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,10,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,11,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,12,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  ENDIF
  !default is CMISS_FIELD_NODE_BASED_INTERPOLATION
  CALL CMISSField_VariableLabelSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"Material",Err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(MaterialsField,Err)

  !Set material parameters (sigma_i = 8.93, sigma_e = 6.7)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    ! sigma_i
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 1,8.93E-6_CMISSDP,Err)  ! 11
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
!      & 2,8.93E-6_CMISSDP,Err)  ! 22
      & 2,0.893E-6_CMISSDP,Err)  ! 22  
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 3,0.0_CMISSDP,Err)  ! 12=21
    ! sigma_e
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 4,6.7E-6_CMISSDP,Err)  ! 11
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 5,6.7E-6_CMISSDP,Err)  ! 22      
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 6,0.0_CMISSDP,Err)  ! 12=21      
  ELSE
    ! sigma_i 
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 1,8.93E-6_CMISSDP,Err)  ! 11
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 2,8.93E-6_CMISSDP,Err)  ! 22
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 3,8.93E-6_CMISSDP,Err)  ! 33
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 4,0.0_CMISSDP,Err)  ! 12=21
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 5,0.0_CMISSDP,Err)  ! 23=32
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 6,0.0_CMISSDP,Err)  !13=31
    ! sigma_e 
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 7,6.7E-6_CMISSDP,Err)  ! 11
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 8,6.7E-6_CMISSDP,Err)  ! 22
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 9,6.7E-6_CMISSDP,Err)  ! 33
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 10,0.0_CMISSDP,Err)  ! 12=21
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 11,0.0_CMISSDP,Err)  ! 23=32
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 12,0.0_CMISSDP,Err)  !13=31    
  ENDIF   

!-------------------------------------------------------------------------
  
  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

!-------------------------------------------------------------------------
!
! SOURCE

  !Start to create a source field on the region
  CALL CMISSField_Initialise(SourceField,Err)
  CALL CMISSField_CreateStart(SourceFieldUserNumber,Region,SourceField,Err)
  CALL CMISSField_TypeSet(SourceField,CMISS_FIELD_GENERAL_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(SourceField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(SourceField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(SourceField,1,Err)
  CALL CMISSField_VariableTypesSet(SourceField,[CMISS_FIELD_U_VARIABLE_TYPE],Err)  
  CALL CMISSField_NumberOfComponentsSet(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,1,Err)   
  CALL CMISSField_VariableLabelSet(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,"Vm",Err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(SourceField,Err)

  !initialise source values to -70.0
  CALL CMISSField_ComponentValuesInitialise(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 1,-70.0_CMISSDP,Err)

  !---------------------------------
  !set yource value on specific nodes
  
  !---------------------------------
  !make Vm a Gaussian bell function along x-direction(s)
  !Vm \in [0,100]  
!  Vm=[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0012_CMISSDP,0.0085_CMISSDP,0.0504_CMISSDP,0.2479_CMISSDP, &
!    & 1.0115_CMISSDP,3.4218_CMISSDP,9.5967_CMISSDP,22.3130_CMISSDP,43.0095_CMISSDP,68.7289_CMISSDP,91.0510_CMISSDP,100.0_CMISSDP, &
!    & 91.0510_CMISSDP,68.7289_CMISSDP,43.0095_CMISSDP,22.3130_CMISSDP,9.5967_CMISSDP,3.4218_CMISSDP,1.0115_CMISSDP, &
!    & 0.2479_CMISSDP,0.0504_CMISSDP,0.0085_CMISSDP,0.0012_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP]

  !Vm \in [-70,30]      
  Vm=[-70.0_CMISSDP,-70.0_CMISSDP,-70.0_CMISSDP,-70.0_CMISSDP,-70.0_CMISSDP,-69.9988_CMISSDP,-69.9915_CMISSDP,-69.9496_CMISSDP, &
    & -69.7521_CMISSDP, &
    & -68.9885_CMISSDP,-66.5782_CMISSDP,-60.4033_CMISSDP,-47.687_CMISSDP,-26.9905_CMISSDP,-1.2711_CMISSDP,21.051_CMISSDP, &
    & 30.0_CMISSDP, &
    & 21.051_CMISSDP,-1.2711_CMISSDP,-26.9905_CMISSDP,-47.687_CMISSDP,-60.4033_CMISSDP,-66.5782_CMISSDP,-68.9885_CMISSDP, &
    & -69.7521_CMISSDP,-69.9496_CMISSDP,-69.9915_CMISSDP,-69.9988_CMISSDP,-70.0_CMISSDP,-70.0_CMISSDP,-70.0_CMISSDP,-70.0_CMISSDP, &
    & -70.0_CMISSDP]    
  
  !x=0.0
  idx=0  
  DO node_idx=1,33
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      idx=idx+1    
      CALL CMISSField_ParameterSetUpdateNode(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
        & 1,node_idx,1,Vm(idx),Err)
    ENDIF
  ENDDO

  !x=0.5
  idx=0  
  DO node_idx=34,66
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      idx=idx+1    
      CALL CMISSField_ParameterSetUpdateNode(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
        & 1,node_idx,1,Vm(idx),Err)
    ENDIF
  ENDDO

  !x=1.0
  idx=0
  DO node_idx=67,99
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      idx=idx+1
      CALL CMISSField_ParameterSetUpdateNode(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
        & 1,node_idx,1,Vm(idx),Err)
    ENDIF
  ENDDO
  
  !x=1.5
  idx=0
  DO node_idx=100,132
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      idx=idx+1
      CALL CMISSField_ParameterSetUpdateNode(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
        & 1,node_idx,1,Vm(idx),Err)
    ENDIF
  ENDDO

  !x=2.0
  idx=0
  DO node_idx=133,165
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      idx=idx+1
      CALL CMISSField_ParameterSetUpdateNode(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
        & 1,node_idx,1,Vm(idx),Err)
    ENDIF
  ENDDO

  !---------------------------------
    
!  !here just middle node (or last node)
!!  CALL CMISSDecomposition_NodeDomainGet(Decomposition,365,1,NodeDomain,Err)
!!  IF(NodeDomain==ComputationalNodeNumber) THEN
!!    CALL CMISSField_ParameterSetUpdateNode(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
!!      & 1,365,1,100.0_CMISSDP,Err)
!!  ENDIF
!  
!  CALL CMISSDecomposition_NodeDomainGet(Decomposition,TotalNumberOfNodes,1,NodeDomain,Err)
!  IF(NodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSField_ParameterSetUpdateNode(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
!      & 1,TotalNumberOfNodes,1,5.0_CMISSDP,Err)
!  ENDIF

  !---------------------------------

  ! loop for any nodes
!  nodes = [25,125];
!  DO node_idx=1,2
!    node = nodes(node_idx)
!    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node,1,NodeDomain,Err)
!    IF(NodeDomain==ComputationalNodeNumber) THEN
!      CALL CMISSField_ParameterSetUpdateNode(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
!        & 1,node,1,100.0_CMISSDP,Err)
!    ENDIF
!  ENDDO

!-------------------------------------------------------------------------
!
! DEPENDENT  

  !Create the equations set dependent field variables
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSField_CreateStart(DependentFieldUserNumber,Region,DependentField,Err) 
  CALL CMISSField_TypeSet(DependentField,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_DependentTypeSet(DependentField,CMISS_FIELD_DEPENDENT_TYPE,Err)  
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentField,2,Err)  
  CALL CMISSField_VariableTypesSet(DependentField,(/CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE/),Err)  
    
  CALL CMISSField_VariableLabelSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,"Phi",Err)        
  CALL CMISSField_VariableLabelSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,"del Phi/del n",Err)        
    
  CALL CMISSField_DimensionSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMISSField_DimensionSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
                              
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,Err)
    
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,Err)   
    
!  !Set the DOFs to be contiguous across components
!  CALL CMISSField_DOFOrderTypeSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
!  CALL CMISSField_DOFOrderTypeSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  
  CALL CMISSField_CreateFinish(DependentField,Err)
  
  !Initialise the field with an initial guess
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.5_CMISSDP, &
    & Err)  

!-------------------------------------------------------------------------
!
! EQUATIONSSET  

  !Create the extracellular bidomain Poisson Equations set
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  
  CALL CMISSField_Initialise(EquationsSetField,Err)
  
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_POISSON_EQUATION_TYPE,CMISS_EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE, &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

!-------------------------------------------------------------------------

  !Create the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

!-------------------------------------------------------------------------

  !Create the equations set material field variables
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set material field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

!-------------------------------------------------------------------------

  !Create the equations set source field variables
  CALL CMISSEquationsSet_SourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  !Finish the equations set material field variables
  CALL CMISSEquationsSet_SourceCreateFinish(EquationsSet,Err)

!--------------------------------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)
  
  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a extracellular bidomain Poisson problem --> "subproblem" of linear source subtype
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS,CMISS_PROBLEM_POISSON_EQUATION_TYPE, &
    & CMISS_PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE,Err)
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
  
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSRegion_NodesGet(Region,Nodes,Err)
  CALL CMISSNodes_NumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL CMISSDecomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL CMISSDecomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
      & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
  ENDIF
!  IF(LastNodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
!      & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
!!      & CMISS_BOUNDARY_CONDITION_FIXED,1.0_CMISSDP,Err)
!  ENDIF
  
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Export results
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"ExtracellularBidomain","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"ExtracellularBidomain","FORTRAN",Err)
  
!--------------------------------------------------------------------------------------------------------------------------------    
  CALL CMISSFields_Finalise(Fields,Err)
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

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
END PROGRAM EXTRACELLULARBIDOMAINEXAMPLE
