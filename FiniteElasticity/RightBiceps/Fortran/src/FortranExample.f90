!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a finite elasticity equation using openCMISS calls.
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
!> The Original Code is openCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): Kumar Mithraratne, Adam Reeve
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

!> \example FiniteElasticity/LargeUniAxialExtension/src/LargeUniAxialExtensionExample.f90
!! Example program to solve a finite elasticity equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/LargeUniAxialExtension/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/LargeUniAxialExtension/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM RIGHTBICEPSACTIVECONTRACTIONEXAMPLE

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

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: NumberOfElements=39
  INTEGER(CMISSIntg), PARAMETER :: NumberOfNodes=567
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussXi=3

  INTEGER(CMISSIntg), PARAMETER :: OutputFrequency=1

  INTEGER(CMISSIntg), DIMENSION(8), PARAMETER :: Entries=[1,3,7,9,19,21,25,27]

  INTEGER(CMISSIntg), DIMENSION(21), PARAMETER :: TendonTopFixNodes=[35,36,37,38,41,42,43,44,119,120,121,122,139,140,179,180, &
    & 181,182,199,200,205]
  INTEGER(CMISSIntg), DIMENSION(21), PARAMETER :: TendonBottomFixNodes=[5,6,7,8,21,22,73,74,89,90,91,92,127,128,149,150,151, &
    & 152,187,188,208]


!  CHARACTER(len=8), PARAMETER :: element_str = "Element:"
  
  !Program types

  !Program variables

  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx

  INTEGER(CMISSIntg) :: FT_1,FT_2,FT_3,FT_4,FT_5,FT_6,FT_7,FT_8,FT_9, &
    & FT_10,FT_11,FT_12,FT_13,FT_14,FT_15,FT_16,FT_17,FT_18,FT_19, &
    & FT_20,FT_21,FT_22,FT_23,FT_24,FT_25,FT_26,FT_27

  REAL(CMISSDP) :: x1,x2,x3
  REAL(CMISSDP) :: VALUE,INIT_PRESSURE

  INTEGER(CMISSIntg) :: i,j,k,elem_idx
  INTEGER(CMISSIntg) :: Elem,Node
  INTEGER(CMISSIntg) :: stat
  character(len=256) :: filename,string


  REAL(CMISSDP), DIMENSION(NumberOfNodes,3) :: AllNodes
  INTEGER(CMISSIntg), DIMENSION(NumberOfElements,27) :: AllElements

!  REAL(CMISSDP), DIMENSION(6) :: MAT_FE
  REAL(CMISSDP), DIMENSION(11) :: MAT_FE

  !CMISS variables
  TYPE(CMISSBasisType) :: Basis, PressureBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,FibreField,MaterialField,DependentField,EquationsSetField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSMeshElementsType) :: QuadraticElements
  TYPE(CMISSMeshElementsType) :: LinearElements

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

  !C(1)=c1_m1...Mooney Rivlin parameter material 1
  !C(2)=c2_m1...Mooney Rivlin parameter material 1
  !C(3)=c4_m1...polynomial coefficient (Markert model) material 1
  !C(4)=c5_m1...power coefficient (Markert model) material 1
  !C(5)=c1_m2...Mooney Rivlin parameter material 2
  !C(6)=c2_m2...Mooney Rivlin parameter material 2
  !C(7)=c4_m2...polynomial coefficient (Markert model) material 2
  !C(8)=c5_m2...power coefficient (Markert model) material 2
  !C(9)=alpha...activation parameter [0,1]
  !C(10)=trans...transition parameter [0,1] for the portion between the two materials
  !C(11)=P_max...maximum isometric stress

  MAT_FE= &
    & [0.0356_CMISSDP,0.00386_CMISSDP, &
    & 0.0000000357_CMISSDP,42.6_CMISSDP, &
    & 2.31_CMISSDP,0.00000115_CMISSDP, &
    & 7.99_CMISSDP,16.6_CMISSDP, &
    & 0.0_CMISSDP,1.0_CMISSDP, &
    & 0.3_CMISSDP]  ! MPa = N/mm^2 = 100 N/cm^2


  WRITE(*,*) "Reading file: input/biceps_3D27.dat"
  OPEN(UNIT=1,FILE="input/biceps_3D27_new.dat",IOSTAT=stat)

  !read the coodinates of the nodes from file
  k=1
  READ(1,*,IOSTAT=stat) 
  READ(1,*,IOSTAT=stat)
  DO WHILE(.true.)
    READ(1,*,IOSTAT=stat) Node,x1,x2,x3
    IF(stat<0) EXIT !end of file
    IF(k/=Node) THEN
      WRITE(*,*) "Error in reading nodal coordinates"
    END IF
    
    AllNodes(k,:)=[x1,x2,x3]
    IF(k==NumberOfNodes) THEN
      WRITE(*,*) "There are ", NumberOfNodes, " nodes in the mesh"
      EXIT
    END IF
    k=k+1
  END DO

  !read the element nodes from file
  k=1
  READ(1,*,IOSTAT=stat) 
  READ(1,*,IOSTAT=stat) 
  READ(1,*,IOSTAT=stat)
  DO WHILE(.true.)
    READ(1,*,IOSTAT=stat) Elem,FT_1,FT_2,FT_3,FT_4,FT_5,FT_6,FT_7,FT_8,FT_9, &
      & FT_10,FT_11,FT_12,FT_13,FT_14,FT_15,FT_16,FT_17,FT_18,FT_19, &
      & FT_20,FT_21,FT_22,FT_23,FT_24,FT_25,FT_26,FT_27

    IF(stat<0) EXIT !end of file
    IF(k/=Elem) THEN
      WRITE(*,*) "Error in reading element nodes"
      EXIT
    END IF
    
    !change order such that Xi_1 rather than Xi_3 direction is aligned with the fibre direction
    AllElements(k,:)=[FT_19,FT_10,FT_1,FT_22,FT_13,FT_4,FT_25,FT_16,FT_7, &
      & FT_20,FT_11,FT_2,FT_23,FT_14,FT_5,FT_26,FT_17,FT_8, &
      & FT_21,FT_12,FT_3,FT_24,FT_15,FT_6,FT_27,FT_18,FT_9]
    IF(k==NumberOfElements) THEN
      WRITE(*,*) "There are ", NumberOfElements, " elements in the mesh"
      EXIT
    END IF
    k=k+1
  END DO

  CLOSE(UNIT=1)
  WRITE(*,*) "Finished reading file: input/biceps_3D27.dat"


  !Intialise cmiss
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  !Set all diganostic levels on for testing
  !CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_RESIDUAL_EVALUATE"],Err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains=NumberOfComputationalNodes


  !Create a 3D rectangular cartesian coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)


  !Create a region and assign the coordinate system to the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_LabelSet(Region,"Region",Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)


  !Define geometric basis
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
  CALL CMISSBasis_InterpolationXiSet(Basis,[CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  CALL CMISSBasis_CreateFinish(Basis,Err)

  !Define pressure basis
  CALL CMISSBasis_Initialise(PressureBasis,Err)
  CALL CMISSBasis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
  CALL CMISSBasis_TypeSet(PressureBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(PressureBasis,3,Err)
  CALL CMISSBasis_InterpolationXiSet(PressureBasis,[CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  CALL CMISSBasis_CreateFinish(PressureBasis,Err)


  !Create a mesh with three-dimensional elements
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,3,Mesh,Err)
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,2,Err) 
  CALL CMISSMesh_NumberOfElementsSet(Mesh,NumberOfElements,Err)  
  !Define nodes for the mesh
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(Region,NumberOfNodes,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)

  CALL CMISSMeshElements_Initialise(QuadraticElements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,1,Basis,QuadraticElements,Err)
  DO j=1,NumberOfElements
    CALL CMISSMeshElements_NodesSet(QuadraticElements,j,AllElements(j,:),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(QuadraticElements,Err)

  !for the linear elements, we need only specific entries from the elements above
  CALL CMISSMeshElements_Initialise(LinearElements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,2,PressureBasis,LinearElements,Err)
  DO j=1,NumberOfElements
    CALL CMISSMeshElements_NodesSet(LinearElements,j,AllElements(j,Entries),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(LinearElements,Err)

  CALL CMISSMesh_CreateFinish(Mesh,Err) 


  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)


  !Create a field to put the geometry (default is geometry)
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSField_VariableLabelSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
!  CALL CMISSField_ScalingTypeSet(GeometricField,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !initialise the Geometry
  DO node_idx=1,SIZE(AllNodes,1)  
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
       & node_idx,1,AllNodes(node_idx,1),Err)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
       & node_idx,2,AllNodes(node_idx,2),Err)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
       & node_idx,3,AllNodes(node_idx,3),Err)
    END IF
  END DO


  !Create a fibre field and attach it to the geometric field
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_VariableLabelSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL CMISSField_CreateFinish(FibreField,Err)

  !Create the dependent field
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSField_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL CMISSField_TypeSet(DependentField,CMISS_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL CMISSField_DependentTypeSet(DependentField,CMISS_FIELD_DEPENDENT_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentField,2,Err)
  CALL CMISSField_VariableTypesSet(DependentField,[CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE],Err)
  CALL CMISSField_VariableLabelSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,Err)
  !Set the pressure to be nodally based and use the second mesh component if required
  CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
    & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,2,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
  CALL CMISSField_CreateFinish(DependentField,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  INIT_PRESSURE=-2.0_CMISSDP*MAT_FE(2)-MAT_FE(1)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4, &
    & INIT_PRESSURE,Err)


  !Create the material field
  CALL CMISSField_Initialise(MaterialField,Err)
  CALL CMISSField_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL CMISSField_TypeSet(MaterialField,CMISS_FIELD_MATERIAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(MaterialField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(MaterialField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialField,1,Err)
  CALL CMISSField_VariableLabelSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,"Material",Err)
!  CALL CMISSField_NumberOfComponentsSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,6,Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,11,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,3,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,5,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,6,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,7,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,8,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,9,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,10,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION, &
    & Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,11,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION, &
    & Err)
  CALL CMISSField_CreateFinish(MaterialField,Err)

  !C(1)=c1_m1...Mooney Rivlin parameter material 1
  !C(2)=c2_m1...Mooney Rivlin parameter material 1
  !C(3)=c4_m1...polynomial coefficient (Markert model) material 1
  !C(4)=c5_m1...power coefficient (Markert model) material 1
  !C(5)=c1_m2...Mooney Rivlin parameter material 2
  !C(6)=c2_m2...Mooney Rivlin parameter material 2
  !C(7)=c4_m2...polynomial coefficient (Markert model) material 2
  !C(8)=c5_m2...power coefficient (Markert model) material 2
  !C(9)=alpha...activation parameter [0,1]
  !C(10)=trans...transition parameter [0,1] for the portion between the two materials
  !C(11)=P_max...maximum isometric stress

  !Set material parameters
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,MAT_FE(1),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,MAT_FE(2),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,MAT_FE(3),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,MAT_FE(4),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,5,MAT_FE(5),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,6,MAT_FE(6),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,7,MAT_FE(7),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,8,MAT_FE(8),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,9,MAT_FE(9),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,10,MAT_FE(10),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,11,MAT_FE(11),Err)

  ! fem group elem  5..8,21..22,28..29,32..36,38..39      --> 100% muscle
  ! default -- do nothing

  ! fem group elem  3..4,13..18,20                        --> 50% muscle
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,10,0.5_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 4,10,0.5_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 13,10,0.5_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 14,10,0.5_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 15,10,0.5_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 16,10,0.5_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 17,10,0.5_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 18,10,0.5_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 20,10,0.5_CMISSDP,Err)

!  INIT_PRESSURE=0.5_CMISSDP*(-2.0_CMISSDP*MAT_FE(2)-MAT_FE(1))+0.5_CMISSDP*(-2.0_CMISSDP*MAT_FE(6)-MAT_FE(5))
!  CALL CMISSField_ParameterSetUpdateNode(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
!    & 3,4,INIT_PRESSURE,Err)

  ! fem group elem  11..12                                --> 80% muscle
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 11,10,0.8_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 12,10,0.8_CMISSDP,Err)

  ! fem group elem  30..31,37                             --> 99% muscle
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 30,10,0.99_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 31,10,0.99_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 37,10,0.99_CMISSDP,Err)

  ! fem group elem  1..2,9..10,19,23                      --> 99.9% muscle
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,10,0.999_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,10,0.999_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 9,10,0.999_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 10,10,0.999_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 19,10,0.999_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 23,10,0.999_CMISSDP,Err)

  ! fem group elem  24..27                                --> 100% soft tissue (all anisotropic contributions = 0)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 24,3,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 24,7,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 25,3,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 25,7,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 26,3,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 26,7,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 27,3,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 27,7,0.0_CMISSDP,Err)


  !Create the equations_set
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
!    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field 
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)


  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)


  !Define the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_ELASTICITY_CLASS,CMISS_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMISS_PROBLEM_NO_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL CMISSControlLoop_TypeSet(ControlLoop,CMISS_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,Err)
  CALL CMISSControlLoop_MaximumIterationsSet(ControlLoop,30,Err)
  CALL CMISSControlLoop_LoadOutputSet(ControlLoop,OutputFrequency,Err)
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(Solver,1.E-6_CMISSDP,Err)
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(Solver,1.E-6_CMISSDP,Err)
  CALL CMISSSolver_NewtonMaximumIterationsSet(Solver,300,Err)
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)


  !Fix the top tendon nodes in all directions
  DO node_idx=1,SIZE(TendonTopFixNodes,1)
    NodeNumber=TendonTopFixNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN

      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 1,VALUE,Err)
      VALUE=VALUE-15.8_CMISSDP
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED,VALUE,Err)

      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 2,VALUE,Err)
      VALUE=VALUE+5.2_CMISSDP
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED,VALUE,Err)

      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 3,VALUE,Err)
      VALUE=VALUE+31.8_CMISSDP
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED,VALUE,Err)
    ENDIF
  ENDDO

  !Fix the bottom tendon nodes in all directions
  DO node_idx=1,SIZE(TendonBottomFixNodes,1)
    NodeNumber=TendonBottomFixNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN

      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 1,VALUE,Err)
      VALUE=VALUE+10.0_CMISSDP
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED,VALUE,Err)

      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 2,VALUE,Err)
      VALUE=VALUE+1.0_CMISSDP
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED,VALUE,Err)

      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 3,VALUE,Err)
      VALUE=VALUE-17.7_CMISSDP
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED,VALUE,Err)
    ENDIF
  ENDDO

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Output solution
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"RightBicepsActiveContraction","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"RightBicepsActiveContraction","FORTRAN",Err)
!  CALL CMISSFields_Finalise(Fields,Err)



  DO i=1,10
  
    VALUE=REAL(i,CMISSDP)/10.0_CMISSDP
    WRITE(*,*) "activation: ", VALUE

    !set the activation level
    DO elem_idx=1,23
      CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,elem_idx,9, &
        & VALUE,ERR)
    END DO
    DO elem_idx=28,NumberOfElements
      CALL CMISSField_ParameterSetUpdateElement(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,elem_idx,9, &
        & VALUE,ERR)
    END DO

    !Solve problem
    CALL CMISSProblem_Solve(Problem,Err)
    
    WRITE(string,'(I2.2)') i
    filename="activation_"//trim(string)
    CALL CMISSFields_NodesExport(Fields,TRIM(filename),"FORTRAN",Err)

  END DO



  !Output solution
!  CALL CMISSFields_Initialise(Fields,Err)
!  CALL CMISSFields_Create(Region,Fields,Err)
!  CALL CMISSFields_NodesExport(Fields,"RightBicepsActiveContraction","FORTRAN",Err)
!  CALL CMISSFields_ElementsExport(Fields,"RightBicepsActiveContraction","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM RIGHTBICEPSACTIVECONTRACTIONEXAMPLE

