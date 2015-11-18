!> \file
!> \author Mylena Mordhorst
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

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif
!  USE CONSTANTS   !for pi

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


!--------------------------------------------------------------------------------------------------------------------------------
  !Test program parameters

!  REAL(CMISSRP), PARAMETER :: WIDTH=6.0_CMISSRP   ! x-direction
!  REAL(CMISSRP), PARAMETER :: HEIGHT=1.9_CMISSRP  ! y-direction
!!  REAL(CMISSRP), PARAMETER :: LENGTH=1.8_CMISSRP  ! z-direction, muscle plus fat/skin
!  REAL(CMISSRP), PARAMETER :: LENGTH=0.8_CMISSRP  ! z-direction, only muscle
!   <-- not needed here (only for generated mesh)
  
  REAL(CMISSRP), PARAMETER :: PI=4.0_CMISSRP*DATAN(1.0_CMISSRP)

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=1
  INTEGER(CMISSIntg), PARAMETER :: MeshComponentNumber=1
      
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

!------------------------------------------------------------------------- 
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS_MUSCLE
  INTEGER(CMISSIntg) :: NumberOfMeshDimensions
  INTEGER(CMISSIntg) :: TotalNumberOfElements,TotalNumberOfNodes,TotalNumberOfElementsMuscle,TotalNumberOfNodesMuscle
  INTEGER(CMISSIntg) :: INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI
  INTEGER(CMISSIntg) :: node_idx,component_idx,elem_idx
  INTEGER(CMISSIntg) :: FibreFieldNumberOfComponents

  REAL(CMISSRP) :: FibreFieldAngle(3)
  INTEGER(CMISSIntg) :: time,dt
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT,Filename
  CHARACTER(LEN=255) :: name_part1,name_part3,numberstring,exnodefile,filename_results

  INTEGER(CMISSIntg),DIMENSION(:,:),ALLOCATABLE :: ElemTopology
  REAL(CMISSRP),DIMENSION(:,:),ALLOCATABLE :: NodeCoords,NodeCoordsMuscle
  
  INTEGER(CMISSIntg) :: clck_counts_beg, clck_counts_begloop, clck_counts_int1, clck_counts_int2, clck_counts_end, clck_rate

!-------------------------------------------------------------------------
  !CMISS variables

  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField,MaterialsField,FibreField,SourceField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_MeshElementsType) :: Elements  
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations

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
  
  CALL SYSTEM_CLOCK (clck_counts_beg, clck_rate)
  !WRITE(*,*) 'begin:', clck_counts_beg

!--------------------------------------------------------------------------------------------------------------------------------
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

!--------------------------------------------------------------------------------------------------------------------------------

  !Input data - either read in as arguments or default specifications
  
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
    NUMBER_GLOBAL_X_ELEMENTS=174
    NUMBER_GLOBAL_Y_ELEMENTS=19
    NUMBER_GLOBAL_Z_ELEMENTS=9

    NUMBER_GLOBAL_Z_ELEMENTS_MUSCLE=4
    TotalNumberOfNodesMuscle=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS_MUSCLE+1)
    TotalNumberOfElementsMuscle=NUMBER_GLOBAL_X_ELEMENTS*NUMBER_GLOBAL_Y_ELEMENTS*NUMBER_GLOBAL_Z_ELEMENTS_MUSCLE
    ALLOCATE(NodecoordsMuscle(TotalNumberOfNodesMuscle,4))           
    
    INTERPOLATION_TYPE=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
!    INTERPOLATION_TYPE=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION   
  ENDIF

!-------------------------------------------------------------------------  
  !Set the total number of nodes and elements 
  
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    NumberOfMeshDimensions=2
    TotalNumberOfElements=NUMBER_GLOBAL_X_ELEMENTS*NUMBER_GLOBAL_Y_ELEMENTS
    IF(INTERPOLATION_TYPE==CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION) THEN
      TotalNumberOfNodes=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
    ELSEIF(INTERPOLATION_TYPE==CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION) THEN
      TotalNumberOfNodes=(NUMBER_GLOBAL_X_ELEMENTS*2+1)*(NUMBER_GLOBAL_Y_ELEMENTS*2+1)    
    ENDIF
    ALLOCATE(Nodecoords(TotalNumberOfNodes,3))
  ELSE
    NumberOfMeshDimensions=3
    TotalNumberOfElements=NUMBER_GLOBAL_X_ELEMENTS*NUMBER_GLOBAL_Y_ELEMENTS*NUMBER_GLOBAL_Z_ELEMENTS
    IF(INTERPOLATION_TYPE==CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION) THEN    
      TotalNumberOfNodes=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
    ELSEIF(INTERPOLATION_TYPE==CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION) THEN
      TotalNumberOfNodes=(NUMBER_GLOBAL_X_ELEMENTS*2+1)*(NUMBER_GLOBAL_Y_ELEMENTS*2+1)*(NUMBER_GLOBAL_Z_ELEMENTS*2+1)   
    ENDIF
    ALLOCATE(Nodecoords(TotalNumberOfNodes,4))       
  ENDIF    

!--------------------------------------------------------------------------------------------------------------------------------
  !Intialise OpenCMISS
  
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  CALL cmfe_RandomSeedsSet(9999,Err)
  
  CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"],Err)
  
  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  WRITE(Filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "ExtracellularBidomain",NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS,INTERPOLATION_TYPE
    
  CALL cmfe_OutputSetOn(Filename,Err)

!--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of a new RC COORDINATE SYSTEM
  
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

!-------------------------------------------------------------------------
  !Start the creation of the REGION
  
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"ExtracellularBidomainRegion",Err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

!-------------------------------------------------------------------------
  !Start the creation of a BASIS (default is trilinear lagrange)
  
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_SIMPLEX_TYPE,Err)
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
    CALL cmfe_Basis_NumberOfXiSet(Basis,2,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ELSE
    !Set the basis to be a tri-interpolation basis
    CALL cmfe_Basis_NumberOfXiSet(Basis,3,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis,Err)

!-------------------------------------------------------------------------   
  !Start the creation of a MESH in the region
  
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err) 
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,TotalNumberOfElements,Err)    
  
  !Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,TotalNumberOfNodes,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)

  CALL cmfe_MeshElements_Initialise(Elements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,MeshComponentNumber,Basis,Elements,Err)
  
!  !read in information from existing exelem file for mesh generation 
!  OPEN(UNIT=2,FILE="./data/ELEMENTNODES.TXT")
!  CALL READ_ELEMENTNODES(2,TotalNumberOfElements,ElemTopology)
!  CLOSE(2)

  !create the mesh by calling CREATE_MESH and thus storing the information in ElemTopology
  CALL CREATE_MESH(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,TotalNumberOfElements,ElemTopology)
  
  !Set the nodes belonging to each element
  DO elem_idx=1,TotalNumberOfElements
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      IF(INTERPOLATION_TYPE==CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION) THEN  !4 nodes per element
        CALL cmfe_MeshElements_NodesSet(Elements,ElemTopology(elem_idx,1),[ElemTopology(elem_idx,2),ElemTopology(elem_idx,3), &
          & ElemTopology(elem_idx,4),ElemTopology(elem_idx,5)],Err)
      ELSEIF(INTERPOLATION_TYPE==CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION) THEN !9 nodes per element
        CALL cmfe_MeshElements_NodesSet(Elements,ElemTopology(elem_idx,1),[ElemTopology(elem_idx,2),ElemTopology(elem_idx,3), &
          & ElemTopology(elem_idx,4),ElemTopology(elem_idx,5),ElemTopology(elem_idx,6),ElemTopology(elem_idx,7), &
          & ElemTopology(elem_idx,8),ElemTopology(elem_idx,9),ElemTopology(elem_idx,10)],Err)
      ENDIF !interpolation
    ELSE  ! 3D
      IF(INTERPOLATION_TYPE==CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION) THEN  !8 nodes per element
        CALL cmfe_MeshElements_NodesSet(Elements,ElemTopology(elem_idx,1),[ElemTopology(elem_idx,2),ElemTopology(elem_idx,3), &
          & ElemTopology(elem_idx,4),ElemTopology(elem_idx,5),ElemTopology(elem_idx,6),ElemTopology(elem_idx,7), &
          & ElemTopology(elem_idx,8),ElemTopology(elem_idx,9)],Err)
      ELSEIF(INTERPOLATION_TYPE==CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION) THEN !27 nodes per element
        CALL cmfe_MeshElements_NodesSet(Elements,ElemTopology(elem_idx,1),[ElemTopology(elem_idx,2),ElemTopology(elem_idx,3), &
          & ElemTopology(elem_idx,4),ElemTopology(elem_idx,5),ElemTopology(elem_idx,6),ElemTopology(elem_idx,7), &
          & ElemTopology(elem_idx,8),ElemTopology(elem_idx,9),ElemTopology(elem_idx,10),ElemTopology(elem_idx,11), &
          & ElemTopology(elem_idx,12),ElemTopology(elem_idx,13),ElemTopology(elem_idx,14),ElemTopology(elem_idx,15), &
          & ElemTopology(elem_idx,16),ElemTopology(elem_idx,17),ElemTopology(elem_idx,18),ElemTopology(elem_idx,19), &
          & ElemTopology(elem_idx,20),ElemTopology(elem_idx,21),ElemTopology(elem_idx,22),ElemTopology(elem_idx,23), &
          & ElemTopology(elem_idx,24),ElemTopology(elem_idx,25),ElemTopology(elem_idx,26),ElemTopology(elem_idx,27), &
          & ElemTopology(elem_idx,28)],Err)
      ENDIF !interpolation          
    ENDIF !dimension
  ENDDO
  
  CALL cmfe_MeshElements_CreateFinish(Elements,Err)

  CALL cmfe_Mesh_CreateFinish(Mesh,Err)   
  
!--------------------------------------------------------------------------------------------------------------------------------
  !Create a DECOMPOSITION
  
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)
 
  !Destory the mesh now that we have decomposed it
!  CALL cmfe_Mesh_Destroy(Mesh,Err)

!--------------------------------------------------------------------------------------------------------------------------------
! FIELDS
!--------------------------------------------------------------------------------------------------------------------------------
!
! GEOMETRIC 
 
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)
  
  !the coordinates for each nodes will be specified in the time-loop around cmfe_Problem_Solve at the bottom
  !of this file as for Thomas case, the nodal positions will also change with time and should be updated
  
  !due to more elements for fat/skin-layer read exnode file for all nodal coordinates (x,y,z)
  OPEN(UNIT=4,FILE="./testinputdata/simple_geometryExample_M.part0.exnode")
!  OPEN(UNIT=4,FILE="../inputdata/data_Vm_2/simple_geometryExample_M.part0.exnode")
  CALL READ_EXNODE(4,TotalNumberOfNodes,NodeCoords)
  CLOSE(4)
    
  !initialise nodal coordinates, GEOMETRIC field 
  !for each node, specify the position in space
  DO node_idx=1,TotalNumberOfNodes
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
!      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!        & NodeCoords(node_idx,1),1,NodeCoords(node_idx,2),Err) !x-positon
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & node_idx,1,NodeCoords(node_idx,1),Err) !x-positon
!      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!        & NodeCoords(node_idx,1),2,NodeCoords(node_idx,3),Err) !y-position
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & node_idx,2,NodeCoords(node_idx,2),Err) !y-position
      IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
!        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!        & NodeCoords(node_idx,1),3,NodeCoords(node_idx,4),Err) !z-position
        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
          & node_idx,3,NodeCoords(node_idx,3),Err) !z-position
      ENDIF
    ENDIF
  ENDDO      
  
  
!-------------------------------------------------------------------------
!
! FIBRE  
  
  !Start to create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FibreFieldUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  !Set the decomposition to use  
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FibreField,FibreFieldNumberOfVariables,Err)
  CALL cmfe_Field_VariableTypesSet(FibreField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)

  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    FibreFieldNumberOfComponents=2  
  ELSE
    FibreFieldNumberOfComponents=3
  ENDIF    
  
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,FibreFieldNumberOfComponents,Err) 
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)   ! as 2D = 1, after if
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field  
  CALL cmfe_Field_CreateFinish(FibreField,Err)

!  CALL cmfe_Nodes_Initialise(Nodes,Err)
!  CALL cmfe_Region_NodesGet(Region,Nodes,Err)
!  CALL cmfe_Nodes_NumberOfNodesGet(Nodes,TotalNumberOfNodes,Err)
  
  !Rotation Angles (in radiant!!)
  ! in 2D an entry in Angle(1) means rotated x-axis, 
  !          entry in Angle(2) doesn't make sense, as rotates out of surface ...
  ! in 3D an entry in Angle(1) means rotated around z-axis, entry in Angle(2) means rotated around y-axis
  !          entry in Angle(3) means rotated around x-axis => no change
  ! 45° equivalent to pi/4, 90° equivalent to pi/2
  
!  FibreFieldAngle=[PI/4.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP]
  FibreFieldAngle=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP]

  DO component_idx=1,FibreFieldNumberOfComponents
    CALL cmfe_Field_ComponentValuesInitialise(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,component_idx, &
      & FibreFieldAngle(component_idx),Err)
  ENDDO

!-------------------------------------------------------------------------
!
! MATERIAL  
  
  !Start to create a material field on the region
  CALL cmfe_Field_Initialise(MaterialsField,Err)
  CALL cmfe_Field_CreateStart(MaterialsFieldUserNumber,Region,MaterialsField,Err)
  CALL cmfe_Field_TypeSet(MaterialsField,CMFE_FIELD_MATERIAL_TYPE,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(MaterialsField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialsField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialsField,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    ! need 2 tensors, sigma_i and sigma_e
    ! symmetric 2x2 tensor => 3 different entries
    ! 1 - 11, 2 - 22, 3 - 12=21
    ! first i then e
    CALL cmfe_Field_NumberOfComponentsSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,6,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,2, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,3, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,4, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,5, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,6, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ELSE
    ! need 2 tensors, sigma_i and sigma_e
    ! symmetric 3x3 tensor => 6 different entries
    ! 1 - 11, 2 - 22, 3 - 33, 4 - 12=21, 5 - 23=32, 6 - 13=31
    ! first i then e 
    CALL cmfe_Field_NumberOfComponentsSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,12,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,2, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,3, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,4, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,5, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,6, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,7, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,8, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,9, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,10, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,11, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,12, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !default is CMFE_FIELD_NODE_BASED_INTERPOLATION
  CALL cmfe_Field_VariableLabelSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(MaterialsField,Err)

  !Set material parameters (sigma_i = 8.93, sigma_e = 6.7)
  !Initialisation values stay as they are on muscle elements
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    ! sigma_i
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 1,8.93E-6_CMISSRP,Err)  ! 11
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
!      & 2,8.93E-6_CMISSRP,Err)  ! 22
      & 2,0.893E-6_CMISSRP,Err)  ! 22  
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 3,0.0_CMISSRP,Err)  ! 12=21
    ! sigma_e
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 4,6.7E-6_CMISSRP,Err)  ! 11
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 5,6.7E-6_CMISSRP,Err)  ! 22      
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 6,0.0_CMISSRP,Err)  ! 12=21      
  ELSE
    ! sigma_i 
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 1,8.93E-6_CMISSRP,Err)  ! 11
!      & 1,0.3E-6_CMISSRP,Err)  ! 11
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
!      & 2,8.93E-6_CMISSRP,Err)  ! 22
      & 2,0.893E-6_CMISSRP,Err)  ! 22
!      & 2,0.0_CMISSRP,Err)  ! 22      
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
!      & 3,8.93E-6_CMISSRP,Err)  ! 33
      & 3,0.893E-6_CMISSRP,Err)  ! 33
!      & 3,0.0E-6_CMISSRP,Err)  ! 33
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 4,0.0_CMISSRP,Err)  ! 12=21
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 5,0.0_CMISSRP,Err)  ! 23=32
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 6,0.0_CMISSRP,Err)  !13=31
    ! sigma_e 
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 7,6.7E-6_CMISSRP,Err)  ! 11
!      & 7,2.2E-6_CMISSRP,Err)  ! 11
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 8,6.7E-6_CMISSRP,Err)  ! 22
!      & 8,2.2E-6_CMISSRP,Err)  ! 22
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 9,6.7E-6_CMISSRP,Err)  ! 33
!      & 9,2.2E-6_CMISSRP,Err)  ! 33
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 10,0.0_CMISSRP,Err)  ! 12=21
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 11,0.0_CMISSRP,Err)  ! 23=32
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 12,0.0_CMISSRP,Err)  !13=31    
  ENDIF
  
  !On fat/skin-layer elements, material parameters are changed 
  !"sigma=sigma_i+sigma_e" => set sigma_i default to 0, thus sigma=sigma_e
  DO elem_idx=TotalNumberOfElementsMuscle+1,TotalNumberOfElements
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !sigma_i --> set to 0          
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,1,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,2,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,3,0.0_CMISSRP,Err)
      !sigma_e=sigma, literature value for fat/skin:
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,4,0.4E-6_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,5,0.4E-6_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,6,0.0_CMISSRP,Err)
    ELSE
      !sigma_i --> set to 0          
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,1,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,2,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,3,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,4,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,5,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,6,0.0_CMISSRP,Err)
      !sigma_e=sigma, literature value for fat/skin:
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,7,0.4E-6_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,8,0.4E-6_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,9,0.4E-6_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,10,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,11,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & elem_idx,12,0.0_CMISSRP,Err)            
    ENDIF
  ENDDO !elem_idx

!-------------------------------------------------------------------------
!
! SOURCE

  !Start to create a source field on the region
  CALL cmfe_Field_Initialise(SourceField,Err)
  CALL cmfe_Field_CreateStart(SourceFieldUserNumber,Region,SourceField,Err)
  CALL cmfe_Field_TypeSet(SourceField,CMFE_FIELD_GENERAL_TYPE,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(SourceField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(SourceField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(SourceField,1,Err)
  CALL cmfe_Field_VariableTypesSet(SourceField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
  CALL cmfe_Field_NumberOfComponentsSet(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)   
  CALL cmfe_Field_VariableLabelSet(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,"Vm",Err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(SourceField,Err)

  !initialise source values to -82.0  (leave like that for fat/skin, then RHS = 0 as divergence of a constant)
  CALL cmfe_Field_ComponentValuesInitialise(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 1,-82.0_CMISSRP,Err)

  !the value of Vm for each node (in muscle) will be specified in the time-loop around cmfe_Problem_Solve at the bottom
  !of this file as Vm will change with time and should be updated

!-------------------------------------------------------------------------
!
! DEPENDENT  

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_Field_CreateStart(DependentFieldUserNumber,Region,DependentField,Err) 
  CALL cmfe_Field_TypeSet(DependentField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_DependentTypeSet(DependentField,CMFE_FIELD_DEPENDENT_TYPE,Err)  
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentField,2,Err)  
  CALL cmfe_Field_VariableTypesSet(DependentField,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE],Err)
    
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Phi",Err)        
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"del Phi/del n",Err)        
    
  CALL cmfe_Field_DimensionSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_DimensionSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
                              
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,Err)
    
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,Err)   
    
  !Set the DOFs to be contiguous across components
!  CALL cmfe_Field_DOFOrderTypeSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
!  CALL cmfe_Field_DOFOrderTypeSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  
  CALL cmfe_Field_CreateFinish(DependentField,Err)
  
  !Initialise the field with an initial guess
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.5_CMISSRP, &
    & Err)  

!-------------------------------------------------------------------------
!
! EQUATIONSSET  

  !Create the extracellular bidomain Poisson Equations set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_POISSON_EQUATION_TYPE,CMFE_EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE], &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

!-------------------------------------------------------------------------

  !Create the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

!-------------------------------------------------------------------------

  !Create the equations set material field variables
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set material field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

!-------------------------------------------------------------------------

  !Create the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  !Finish the equations set material field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(EquationsSet,Err)

!--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set EQUATIONS
  
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

!-------------------------------------------------------------------------  
  !Start the creation of a PROBLEM
  
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS,CMFE_PROBLEM_POISSON_EQUATION_TYPE, &
    & CMFE_PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

!-------------------------------------------------------------------------   
  !Start the creation of the problem SOLVERS
  
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  
!  CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
!  CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-12_CMISSRP,Err)
!  CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(Solver,1.0E-12_CMISSRP,Err)
!  CALL cmfe_Solver_LinearIterativeTypeSet(Solver,CMFE_SOLVER_ITERATIVE_BiCGSTAB,Err)
!!  CALL cmfe_Solver_LinearIterativeTypeSet(Solver,CMFE_SOLVER_ITERATIVE_CONJUGATE_GRADIENT,Err)
!!  !CALL cmfe_Solver_LinearIterativeTypeSet(Solver,CMFE_SOLVER_ITERATIVE_GMRES,Err)
!  CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(Solver,1000000000,Err)
  
  CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  
  !CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  !CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_LAPACK_LIBRARY,Err)
  !CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_SUPERLU_LIBRARY,Err)
  !CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_PASTIX_LIBRARY,Err)
  
!  CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_PETSC_LIBRARY,Err)  
!  CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_HYPRE_LIBRARY,Err)   !<Hypre solver library.
  
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

!-------------------------------------------------------------------------   
  !Start the creation of the problem solver equations
  
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

!-------------------------------------------------------------------------
  !Start the creation of the equations set BOUNDARY CONDITIONS
  
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  
  !Set the (first node and the) last node to 0.0
  FirstNodeNumber=1
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Region_NodesGet(Region,Nodes,Err)
  CALL cmfe_Nodes_NumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
!  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
!    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
!      & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
!  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
  ENDIF
  
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

!-------------------------------------------------------------------------

  time=24
  dt=2
  DO WHILE(time<=26)
    
    CALL SYSTEM_CLOCK (clck_counts_begloop, clck_rate)
    WRITE (*,*) 'begloop-beg',  (clck_counts_begloop - clck_counts_beg) / REAL(clck_rate)
        
    name_part1="./testinputdata/MainTime_M_2_"
!    name_part1="../inputdata/data_Vm_2/MainTime_M_2_"
    name_part3=".part0.exnode"
    WRITE(numberstring,*) time
    numberstring=ADJUSTL(numberstring)
    exnodefile=TRIM(name_part1)//TRIM(numberstring)//TRIM(name_part3)
    
    WRITE(*,*) exnodefile
    
    !in each time step, read in information from existing exnode file for nodal information (node,x,y,z,Vm)
    OPEN(UNIT=3,FILE=exnodefile)    
!    CALL READ_EXNODE(3,TotalNumberOfNodes,NodeCoords)
    CALL READ_EXNODE(3,TotalNumberOfNodesMuscle,NodeCoordsMuscle)    
    CLOSE(3)
    
!    !update nodal coordinates, GEOMETRIC field (won't change in my case but later important for Thomas)
!    !for each node, specify the position in space
!    DO node_idx=1,TotalNumberOfNodes
!      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
!      IF(NodeDomain==ComputationalNodeNumber) THEN
!!        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!!          & NodeCoords(node_idx,1),1,NodeCoords(node_idx,2),Err) !x-positon
!        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!          & node_idx,1,NodeCoords(node_idx,1),Err) !x-positon
!!        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!!          & NodeCoords(node_idx,1),2,NodeCoords(node_idx,3),Err) !y-position
!        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!          & node_idx,2,NodeCoords(node_idx,2),Err) !y-position
!        IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
!!          CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!!          & NodeCoords(node_idx,1),3,NodeCoords(node_idx,4),Err) !z-position
!          CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!          & node_idx,3,NodeCoords(node_idx,3),Err) !z-position
!        ENDIF
!      ENDIF
!    ENDDO    
    
    !update the SOURCE field Vm
    !set source value on specific nodes  
    DO node_idx=1,TotalNumberOfNodesMuscle
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
!          CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
!            & 1,NodeCoordsMuscle(node_idx,1),1,NodeCoordsMuscle(node_idx,4),Err)
          CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
            & 1,node_idx,1,NodeCoordsMuscle(node_idx,3),Err)
        ELSE
!          CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
!            & 1,NodeCoordsMuscle(node_idx,1),1,NodeCoordsMuscle(node_idx,5),Err)
          CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
            & 1,node_idx,1,NodeCoordsMuscle(node_idx,4),Err)
        ENDIF
      ENDIF
    ENDDO    
    
    CALL SYSTEM_CLOCK (clck_counts_int1, clck_rate)
    !WRITE(*,*) 'int1:', clck_counts_int1
    WRITE (*,*) 'int1-begloop',  (clck_counts_int1 - clck_counts_begloop) / REAL(clck_rate)
    
    !Solve the problem
    CALL cmfe_Problem_Solve(Problem,Err)
    
    CALL SYSTEM_CLOCK (clck_counts_int2, clck_rate)
    !WRITE(*,*) 'int2', clck_counts_int2
    WRITE (*,*) 'int2-int1', (clck_counts_int2 - clck_counts_int1) / REAL(clck_rate)

    filename_results="ExtracellularBidomain_"//TRIM(numberstring)
  
    WRITE(*,*) filename_results
      
    !Export results
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,filename_results,"FORTRAN",Err)
!    CALL cmfe_Fields_ElementsExport(Fields,filename_results,"FORTRAN",Err)    
    
!    CALL cmfe_Fields_NodesExport(Fields,"ExtracellularBidomain","FORTRAN",Err)
!    CALL cmfe_Fields_ElementsExport(Fields,"ExtracellularBidomain","FORTRAN",Err)

    time=time+dt
    
    CALL SYSTEM_CLOCK (clck_counts_end, clck_rate)
    !WRITE(*,*) 'end', clck_counts_end
    WRITE (*,*) 'end-int2', (clck_counts_end - clck_counts_int2) / REAL(clck_rate)

  ENDDO !time
  
!--------------------------------------------------------------------------------------------------------------------------------    
  CALL cmfe_Fields_Finalise(Fields,Err)
  
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  DEALLOCATE(ElemTopology) 
  DEALLOCATE(NodeCoords)
  DEALLOCATE(NodeCoordsMuscle)  

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
  
  !--------------------------------------------------------------------------------------
  !Subroutine to read in mesh information from ELEMENTNODES.TXT (created by gen_elem.f90)
  !-------------------------------------------------------------------------------------- 

!  SUBROUTINE READ_ELEMENTNODES(fp,totalnumberofelements,elemtopology)
!    INTEGER(CMISSIntg), INTENT(IN)    :: fp  !< file 'pointer'
!    INTEGER(CMISSIntg), INTENT(IN)    :: totalnumberofelements !< total number of elements, needed for array size    
!    INTEGER(CMISSIntg), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: elemtopology !< element topology, array with rows 
!                                                                      !that contain element number and corresponding nodes
!    INTEGER(CMISSIntg) :: elemnr,n1,n2,n3,n4,n5,n6,n7,n8
!    INTEGER(CMISSIntg) :: k,j
!    
!    !select correct size of array depending on 2d/3d and linear/quadratic s.f.
!    !ALLOCATE(elemtopology(totalnumberofelements,5)) !2d,linear - 4
!    !ALLOCATE(elemtopology(totalnumberofelements,10)) !2d,quadratic - 9
!    ALLOCATE(elemtopology(totalnumberofelements,9)) !3d,linear - 8
!    !ALLOCATE(elemtopology(totalnumberofelements,28)) !3d,quadratic - 27

!    DO k=1,totalnumberofelements
!      READ(fp,*) elemnr,n1,n2,n3,n4,n5,n6,n7,n8
!      !WRITE(*,*) elemnr,n1,n2,n3,n4,n5,n6,n7,n8
!      elemtopology(k,:)=[elemnr,n1,n2,n3,n4,n5,n6,n7,n8]
!    ENDDO

!    
!!    WRITE(*,*) "Values for Elements are:"
!!    DO k=1,totalnumberofelements
!!      DO j=1,9
!!        WRITE(*,*) elemtopology(k,j)
!!      ENDDO
!!    ENDDO

!  END SUBROUTINE READ_ELEMENTNODES

  !--------------------------------------------------------------------------------------
  !Subroutine to create the mesh by the given geometry
  !-------------------------------------------------------------------------------------- 

  SUBROUTINE CREATE_MESH(n_elem_x,n_elem_y,n_elem_z,totalnumberofelements,elemtopology)
    INTEGER(CMISSIntg), INTENT(IN) :: n_elem_x,n_elem_y,n_elem_z    !< total number of elements in each direction
    INTEGER(CMISSIntg), INTENT(IN) :: totalnumberofelements     !< total number of elements    
    INTEGER(CMISSIntg), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: elemtopology !< element topology, array with rows 
                                                                      !that contain element number and corresponding nodes
    INTEGER(CMISSIntg) :: n_nodes_x,n_fibres_y,n_fibres_z,idx_x,idx_y,idx_z                                                     
    INTEGER(CMISSIntg) :: elemnr=0,n1,n2,n3,n4,n5,n6,n7,n8

    ALLOCATE(elemtopology(totalnumberofelements,9)) !3d,linear - 8
  
    n_nodes_x=n_elem_x+1
    n_fibres_y=n_elem_y+1
    n_fibres_z=n_elem_z+1
  
    DO idx_z=1,n_elem_z
      DO idx_y=1,n_elem_y
        DO idx_x=1,n_elem_x
          elemnr=elemnr+1
          n1=idx_x+(idx_y-1)*n_nodes_x+(idx_z-1)*n_nodes_x*n_fibres_y
          n2=n1+1
          n3=n1+n_nodes_x
          n4=n3+1
          n5=n1+n_nodes_x*n_fibres_y
          n6=n5+1
          n7=n5+n_nodes_x
          n8=n7+1
          elemtopology(elemnr,:)=[elemnr,n1,n2,n3,n4,n5,n6,n7,n8]
        ENDDO !idx_x
      ENDDO !idx_y
    ENDDO !idx_z
    
  END SUBROUTINE CREATE_MESH
  
  !--------------------------------------------------------------------------------------
  !Subroutine to read in nodal information from Thomas . (created by gen_elem.f90)
  !--------------------------------------------------------------------------------------
  
  SUBROUTINE READ_EXNODE(fp,totalnumberofnodes,nodecoords)
    INTEGER(CMISSIntg), INTENT(IN)    :: fp  !< file 'pointer'
    INTEGER(CMISSIntg), INTENT(IN)    :: totalnumberofnodes !< total number of nodes, needed for array size
    REAL(CMISSRP), DIMENSION(totalnumberofnodes,4), INTENT(INOUT) :: nodecoords !< nodal coordinates and Vm, array with rows
                                                                      !that contain (node number), corresponding coordinates 
                                                                      !and source value Vm for that node
    
!    REAL(CMISSRP), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: nodecoords !< nodal coordinates and Vm, array with rows 
                                                                      !that contain (node number), corresponding coordinates 
                                                                      !and source value Vm for that node
    REAL(CMISSRP) :: geoM1Dx,mat1,mat2,mat3,Vm,dVmdt,geoM3Dx,geoM3Dy,geoM3Dz
    INTEGER(CMISSIntg) :: iost,k,j,nodenr
    CHARACTER(len=256) :: str1

!    !select correct size of array depending on 2d/3d
!    !ALLOCATE(nodecoords(totalnumberofnodes,3)) !2d
!    ALLOCATE(nodecoords(totalnumberofnodes,4)) !3d
    
    nodenr=0
    nodecoords=0.0_CMISSRP
    
    ReadLoopNode: DO
      READ(fp,*,IOSTAT=iost) str1 !read data line by line
      IF (iost .LT. 0) EXIT !end of file reached before Node was found
      IF (TRIM(str1)=='Node:') THEN !prepare to begin capturing
        nodenr=nodenr+1
        READ(fp,*) geoM1Dx
        READ(fp,*) mat1
        READ(fp,*) mat2
        READ(fp,*) mat3
        READ(fp,*) Vm
        READ(fp,*) dVmdt
        READ(fp,*) geoM3Dx
        READ(fp,*) geoM3Dy
        READ(fp,*) geoM3Dz
        !WRITE(*,*) nodenr,geoM3Dx,geoM3Dy,geoM3Dz,Vm
        nodecoords(nodenr,:)=[geoM3Dx,geoM3Dy,geoM3Dz,Vm]
      ENDIF
    ENDDO ReadLoopNode
    
    
!    WRITE(*,*) "Values for nodes are:"
!    DO k=1,totalnumberofnodes
!      DO j=1,4
!        WRITE(*,*) nodecoords(k,j)
!      ENDDO
!    ENDDO
  
  END SUBROUTINE READ_EXNODE

  !------------------------------------------------------------------------------------------------------------------------------  
END PROGRAM EXTRACELLULARBIDOMAINEXAMPLE
