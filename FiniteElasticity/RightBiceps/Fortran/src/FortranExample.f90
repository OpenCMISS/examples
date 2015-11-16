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
!! Example program to solve a finite elasticity equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/LargeUniAxialExtension/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/LargeUniAxialExtension/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM RIGHTBICEPSACTIVECONTRACTIONEXAMPLE

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

  REAL(CMISSRP) :: x1,x2,x3
  REAL(CMISSRP) :: VALUE,INIT_PRESSURE

  INTEGER(CMISSIntg) :: i,j,k,elem_idx
  INTEGER(CMISSIntg) :: Elem,Node
  INTEGER(CMISSIntg) :: stat
  character(len=256) :: filename,string


  REAL(CMISSRP), DIMENSION(NumberOfNodes,3) :: AllNodes
  INTEGER(CMISSIntg), DIMENSION(NumberOfElements,27) :: AllElements

!  REAL(CMISSRP), DIMENSION(6) :: MAT_FE
  REAL(CMISSRP), DIMENSION(11) :: MAT_FE

  !CMISS variables
  TYPE(cmfe_BasisType) :: Basis, PressureBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,FibreField,MaterialField,DependentField,EquationsSetField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_MeshElementsType) :: QuadraticElements
  TYPE(cmfe_MeshElementsType) :: LinearElements

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
    & [0.0356_CMISSRP,0.00386_CMISSRP, &
    & 0.0000000357_CMISSRP,42.6_CMISSRP, &
    & 2.31_CMISSRP,0.00000115_CMISSRP, &
    & 7.99_CMISSRP,16.6_CMISSRP, &
    & 0.0_CMISSRP,1.0_CMISSRP, &
    & 0.3_CMISSRP]  ! MPa = N/mm^2 = 100 N/cm^2


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
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  !Set all diganostic levels on for testing
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_RESIDUAL_EVALUATE"],Err)

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains=NumberOfComputationalNodes


  !Create a 3D rectangular cartesian coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)


  !Create a region and assign the coordinate system to the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"Region",Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)


  !Define geometric basis
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(Basis,3,Err)
  CALL cmfe_Basis_InterpolationXiSet(Basis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  CALL cmfe_Basis_CreateFinish(Basis,Err)

  !Define pressure basis
  CALL cmfe_Basis_Initialise(PressureBasis,Err)
  CALL cmfe_Basis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
  CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(PressureBasis,3,Err)
  CALL cmfe_Basis_InterpolationXiSet(PressureBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  CALL cmfe_Basis_CreateFinish(PressureBasis,Err)


  !Create a mesh with three-dimensional elements
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,3,Mesh,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,2,Err) 
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NumberOfElements,Err)  
  !Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,NumberOfNodes,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)

  CALL cmfe_MeshElements_Initialise(QuadraticElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,1,Basis,QuadraticElements,Err)
  DO j=1,NumberOfElements
    CALL cmfe_MeshElements_NodesSet(QuadraticElements,j,AllElements(j,:),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(QuadraticElements,Err)

  !for the linear elements, we need only specific entries from the elements above
  CALL cmfe_MeshElements_Initialise(LinearElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,2,PressureBasis,LinearElements,Err)
  DO j=1,NumberOfElements
    CALL cmfe_MeshElements_NodesSet(LinearElements,j,AllElements(j,Entries),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(LinearElements,Err)

  CALL cmfe_Mesh_CreateFinish(Mesh,Err) 


  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)


  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
!  CALL cmfe_Field_ScalingTypeSet(GeometricField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !initialise the Geometry
  DO node_idx=1,SIZE(AllNodes,1)  
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
       & node_idx,1,AllNodes(node_idx,1),Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
       & node_idx,2,AllNodes(node_idx,2),Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
       & node_idx,3,AllNodes(node_idx,3),Err)
    END IF
  END DO


  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !Create the dependent field
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL cmfe_Field_TypeSet(DependentField,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL cmfe_Field_DependentTypeSet(DependentField,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentField,2,Err)
  CALL cmfe_Field_VariableTypesSet(DependentField,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE],Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,Err)
  !Set the pressure to be nodally based and use the second mesh component if required
  CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,2,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
  CALL cmfe_Field_CreateFinish(DependentField,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
  INIT_PRESSURE=-2.0_CMISSRP*MAT_FE(2)-MAT_FE(1)
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
    & INIT_PRESSURE,Err)


  !Create the material field
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL cmfe_Field_TypeSet(MaterialField,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialField,1,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
!  CALL cmfe_Field_NumberOfComponentsSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,6,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,11,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,5,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,6,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,7,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,8,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,9,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,10,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION, &
    & Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,11,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION, &
    & Err)
  CALL cmfe_Field_CreateFinish(MaterialField,Err)

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
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,MAT_FE(1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,MAT_FE(2),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,MAT_FE(3),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,MAT_FE(4),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,MAT_FE(5),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,MAT_FE(6),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,7,MAT_FE(7),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,8,MAT_FE(8),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,9,MAT_FE(9),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,10,MAT_FE(10),Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,11,MAT_FE(11),Err)

  ! fem group elem  5..8,21..22,28..29,32..36,38..39      --> 100% muscle
  ! default -- do nothing

  ! fem group elem  3..4,13..18,20                        --> 50% muscle
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 4,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 13,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 14,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 15,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 16,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 17,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 18,10,0.5_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 20,10,0.5_CMISSRP,Err)

!  INIT_PRESSURE=0.5_CMISSRP*(-2.0_CMISSRP*MAT_FE(2)-MAT_FE(1))+0.5_CMISSRP*(-2.0_CMISSRP*MAT_FE(6)-MAT_FE(5))
!  CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!    & 3,4,INIT_PRESSURE,Err)

  ! fem group elem  11..12                                --> 80% muscle
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 11,10,0.8_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 12,10,0.8_CMISSRP,Err)

  ! fem group elem  30..31,37                             --> 99% muscle
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 30,10,0.99_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 31,10,0.99_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 37,10,0.99_CMISSRP,Err)

  ! fem group elem  1..2,9..10,19,23                      --> 99.9% muscle
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,10,0.999_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,10,0.999_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 9,10,0.999_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 10,10,0.999_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 19,10,0.999_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 23,10,0.999_CMISSRP,Err)

  ! fem group elem  24..27                                --> 100% soft tissue (all anisotropic contributions = 0)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 24,3,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 24,7,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 25,3,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 25,7,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 26,3,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 26,7,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 27,3,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 27,7,0.0_CMISSRP,Err)


  !Create the equations_set
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE], &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field 
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)


  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)


  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_NO_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL cmfe_ControlLoop_TypeSet(ControlLoop,CMFE_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoop,30,Err)
  CALL cmfe_ControlLoop_LoadOutputSet(ControlLoop,OutputFrequency,Err)
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
  CALL cmfe_Solver_NewtonRelativeToleranceSet(Solver,1.E-6_CMISSRP,Err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(Solver,1.E-6_CMISSRP,Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(Solver,300,Err)
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)


  !Fix the top tendon nodes in all directions
  DO node_idx=1,SIZE(TendonTopFixNodes,1)
    NodeNumber=TendonTopFixNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN

      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 1,VALUE,Err)
      VALUE=VALUE-15.8_CMISSRP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,VALUE,Err)

      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 2,VALUE,Err)
      VALUE=VALUE+5.2_CMISSRP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,VALUE,Err)

      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 3,VALUE,Err)
      VALUE=VALUE+31.8_CMISSRP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,VALUE,Err)
    ENDIF
  ENDDO

  !Fix the bottom tendon nodes in all directions
  DO node_idx=1,SIZE(TendonBottomFixNodes,1)
    NodeNumber=TendonBottomFixNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN

      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 1,VALUE,Err)
      VALUE=VALUE+10.0_CMISSRP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,VALUE,Err)

      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 2,VALUE,Err)
      VALUE=VALUE+1.0_CMISSRP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,VALUE,Err)

      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
        & 3,VALUE,Err)
      VALUE=VALUE-17.7_CMISSRP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,VALUE,Err)
    ENDIF
  ENDDO

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !Output solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"RightBicepsActiveContraction","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"RightBicepsActiveContraction","FORTRAN",Err)
!  CALL cmfe_Fields_Finalise(Fields,Err)

  DO i=1,10
  
    VALUE=REAL(i,CMISSRP)/10.0_CMISSRP
    WRITE(*,*) "activation: ", VALUE

    !set the activation level
    DO elem_idx=1,23
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,elem_idx,9, &
        & VALUE,ERR)
    END DO
    DO elem_idx=28,NumberOfElements
      CALL cmfe_Field_ParameterSetUpdateElement(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,elem_idx,9, &
        & VALUE,ERR)
    END DO

    !Solve problem
    CALL cmfe_Problem_Solve(Problem,Err)
    
    WRITE(string,'(I2.2)') i
    filename="activation_"//trim(string)
    CALL cmfe_Fields_NodesExport(Fields,TRIM(filename),"FORTRAN",Err)

  END DO



  !Output solution
!  CALL cmfe_Fields_Initialise(Fields,Err)
!  CALL cmfe_Fields_Create(Region,Fields,Err)
!  CALL cmfe_Fields_NodesExport(Fields,"RightBicepsActiveContraction","FORTRAN",Err)
!  CALL cmfe_Fields_ElementsExport(Fields,"RightBicepsActiveContraction","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM RIGHTBICEPSACTIVECONTRACTIONEXAMPLE

