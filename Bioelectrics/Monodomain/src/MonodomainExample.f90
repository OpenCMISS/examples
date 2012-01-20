!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a Monodomain equation using OpenCMISS calls.
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

!> \example Bioelectrics/Monodomain/src/MonodomainExample.f90
!! Example program to solve a Monodomain equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/Bioelectrics/Monodomain/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/Bioelectrics/Monodomain/build-gnu'>Linux GNU Build</a>
!!
!<

!> Main program
PROGRAM MONODOMAINEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField


  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=3.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=16

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  
  INTEGER(CMISSIntg) :: MPI_IERROR

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: N,CELL_TYPE

  INTEGER(CMISSIntg) :: n98ModelIndex,JRWModelIndex,LRdModelIndex

  INTEGER(CMISSIntg) :: gK1component,gNacomponent,stimcomponent,node_idx

  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_ELEMENTS=25

  REAL(CMISSDP) :: X,Y,DISTANCE,gK1_VALUE,gNa_VALUE
  
  REAL(CMISSDP), PARAMETER :: STIM_VALUE = 100.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: STIM_STOP = 0.10_CMISSDP
  REAL(CMISSDP), PARAMETER :: TIME_STOP = 1.50_CMISSDP
  REAL(CMISSDP), PARAMETER :: ODE_TIME_STEP = 0.00001_CMISSDP
  REAL(CMISSDP), PARAMETER :: PDE_TIME_STEP = 0.001_CMISSDP
  REAL(CMISSDP), PARAMETER :: CONDUCTIVITY = 0.1_CMISSDP

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCellMLType) :: CellML
  TYPE(CMISSCellMLEquationsType) :: CellMLEquations
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField
  TYPE(CMISSFieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
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
  INTEGER(CMISSIntg) :: EquationsSetIndex,CellMLIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain
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

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Trap errors
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
  
  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  CALL CMISSOutputSetOn("Monodomain",Err)
    
  NUMBER_GLOBAL_X_ELEMENTS=NUMBER_OF_ELEMENTS
  NUMBER_GLOBAL_Y_ELEMENTS=NUMBER_OF_ELEMENTS
  NUMBER_GLOBAL_Z_ELEMENTS=0
  NUMBER_OF_DOMAINS=NumberOfComputationalNodes
  
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  !Read in the number of elements in the X & Y directions, and the number of partitions on the master node (number 0)

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
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Set the region label
  CALL CMISSRegion_LabelSet(Region,"Region",Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bilinear Lagrange basis
    CALL CMISSBasis_NumberOfXiSet(Basis,2,Err)
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
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
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)
  
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

  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
        
  !Create the equations_set
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSField_Initialise(EquationsSetField,Err)
  !Set the equations set to be a Monodomain equations set
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,CMISS_EQUATIONS_SET_BIOELECTRICS_CLASS, &
    & CMISS_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMISS_EQUATIONS_SET_NO_SUBTYPE,EquationsSetFieldUserNumber,EquationsSetField, &
    & EquationsSet,Err)
  
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)
  
  !Create the equations set materials field variables
  CALL CMISSField_Initialise(MaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)
  
  !Set Am
  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & 193.6_CMISSDP, &
    & Err)
  !Set Cm
  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
    & 0.014651_CMISSDP,Err)
  !Set conductivity
  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,CONDUCTIVITY, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,CONDUCTIVITY, &
    & Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,5, &
      & CONDUCTIVITY,Err)
  ENDIF

  !Create the CellML environment
  CALL CMISSCellML_Initialise(CellML,Err)
  CALL CMISSCellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import a Noble 1998 model from a file
  CALL CMISSCellML_ModelImport(CellML,"n98.xml",n98ModelIndex,Err)
  ! and import JRW 1998 from a file
  !CALL CMISSCellML_ModelImport(CellML,"jrw-1998.xml",JRWModelIndex,Err)
  ! and import LRd from a file
  !CALL CMISSCellML_ModelImport(CellML,"LRd.xml",LRdModelIndex,Err)
!  CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,[1,2,3,4,5],"",["CELLML_CREATE_FIELD_TO_CELLML_MAP_C", &
!    & "CELLML_CREATE_CELLML_TO_FIELD_MAP_C"],Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  !CALL CMISSCellML_VariableSetAsKnown(CellML,n98ModelIndex,"time_independent_potassium_current/g_K1",Err)
  CALL CMISSCellML_VariableSetAsKnown(CellML,n98ModelIndex,"fast_sodium_current/g_Na ",Err)
  CALL CMISSCellML_VariableSetAsKnown(CellML,n98ModelIndex,"membrane/IStim",Err)
  !CALL CMISSCellML_VariableSetAsKnown(CellML,JRWModelIndex,"L_type_Ca_channel/Ko",Err) ! this one should fail
  !CALL CMISSCellML_VariableSetAsKnown(CellML,JRWModelIndex,"membrane/I_stim",Err)
  !   - to get from the CellML side
  CALL CMISSCellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K1",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_to",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K_ATP",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Ca_L_K",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_b_K",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_NaK",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Na",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_b_Na",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Ca_L_Na",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_NaCa",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/IStimC",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_K1",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_Na",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_Ca_L_Ca",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_Ca_L_K",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_K",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_NaCa",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_Kp",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_p_Ca",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_Na_b",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_Ca_b",Err)
  !CALL CMISSCellML_VariableSetAsWanted(CellML,JRWModelIndex,"membrane/IStimC",Err)
  !   - and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
!  CALL CMISSDiagnosticsSetOff(Err)
  !Finish the CellML environment
  CALL CMISSCellML_CreateFinish(CellML,Err)

!  CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,[1,2,3,4,5],"",["CELLML_CREATE_FIELD_TO_CELLML_MAP_C", &
!    & "CELLML_CREATE_CELLML_TO_FIELD_MAP_C"],Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellML_FieldMapsCreateStart(CellML,Err)
  !Now we can set up the field variable component <--> CellML model variable mappings.
  !Map Vm
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & n98ModelIndex,"membrane/V",CMISS_FIELD_VALUES_SET_TYPE,Err)
  !CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
  !  & JRWModelIndex,"membrane/V",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,n98ModelIndex,"membrane/V",CMISS_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !CALL CMISSCellML_CreateCellMLToFieldMap(CellML,JRWModelIndex,"membrane/V",CMISS_FIELD_VALUES_SET_TYPE, &
  !  & DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellML_FieldMapsCreateFinish(CellML,Err)

  !todo - get vm initialial value.
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & -92.5_CMISSDP, &
    & Err)
  
!  CALL CMISSDiagnosticsSetOff(Err)

  !Start the creation of the CellML models field
  CALL CMISSField_Initialise(CellMLModelsField,Err)
  CALL CMISSCellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL CMISSCellML_ModelsFieldCreateFinish(CellML,Err)
  !Set up the models field
  !DO N=1,(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  !  IF(N < 5) THEN
  !    CELL_TYPE = 1
  !  ELSE
  !    CELL_TYPE = 2
  !  ENDIF
  !  CALL CMISSField_ParameterSetUpdateNode(CellMLModelsField, CMISS_FIELD_U_VARIABLE_TYPE, CMISS_FIELD_VALUES_SET_TYPE,1,N,1,CELL_TYPE,Err)
  !END DO
  !CALL CMISSField_ParameterSetUpdateStart(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !CALL CMISSField_ParameterSetUpdateFinish(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !CALL CMISSField_ComponentValuesInitialise(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2_CMISSIntg,Err)

  !Start the creation of the CellML state field
  CALL CMISSField_Initialise(CellMLStateField,Err)
  CALL CMISSCellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL CMISSCellML_StateFieldCreateFinish(CellML,Err)

  !Start the creation of the CellML intermediate field
  CALL CMISSField_Initialise(CellMLIntermediateField,Err)
  CALL CMISSCellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  !Finish the creation of the CellML intermediate field
  CALL CMISSCellML_IntermediateFieldCreateFinish(CellML,Err)
  
  !Start the creation of CellML parameters field
  CALL CMISSField_Initialise(CellMLParametersField,Err)
  CALL CMISSCellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL CMISSCellML_ParametersFieldCreateFinish(CellML,Err)
  
  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)

  CALL CMISSCellML_FieldComponentGet(CellML,n98ModelIndex,CMISS_CELLML_PARAMETERS_FIELD,"membrane/IStim",stimcomponent,Err)
  !Set the Stimulus at half the bottom nodes
  DO node_idx=1,NUMBER_OF_ELEMENTS/2
    IF(FirstNodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
        & node_idx, &
        & stimcomponent,STIM_VALUE,Err)
    ENDIF
  ENDDO

  !!Set up the g_K1 gradient
  !CALL CMISSCellML_FieldComponentGet(CellML,n98ModelIndex,CMISS_CELLML_PARAMETERS_FIELD,"time_independent_potassium_current/g_K1", &
  !  & gK1component,Err)
  !!Loop over the nodes
  !DO node_idx=1,(NUMBER_OF_ELEMENTS+1)*(NUMBER_OF_ELEMENTS+1)
  !  CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
  !    & X,Err)
  !  CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,2, &
  !    & Y,Err)
  !  DISTANCE=SQRT(X**2+Y**2)/SQRT(2.0_CMISSDP)
  !  gK1_VALUE=2.0_CMISSDP*(DISTANCE+0.5_CMISSDP)*77.11e-3_CMISSDP
  !  CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx, &
  !    & gK1component,gK1_VALUE,Err)
  !ENDDO
  
  !Set up the g_Na gradient
  CALL CMISSCellML_FieldComponentGet(CellML,n98ModelIndex,CMISS_CELLML_PARAMETERS_FIELD,"fast_sodium_current/g_Na", &
    & gNacomponent,Err)
  !Loop over the nodes
  DO node_idx=1,(NUMBER_OF_ELEMENTS+1)*(NUMBER_OF_ELEMENTS+1)
    CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
      & X,Err)
    CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,node_idx,2, &
      & Y,Err)
    DISTANCE=SQRT(X**2+Y**2)/SQRT(2.0_CMISSDP)
    gNa_VALUE=2.0_CMISSDP*(DISTANCE+0.5_CMISSDP)*385.5e-3_CMISSDP
    CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
      & node_idx, &
      & gNacomponent,gNa_VALUE,Err)
  ENDDO
  
  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_BIOELECTRICS_CLASS,CMISS_PROBLEM_MONODOMAIN_EQUATION_TYPE, &
    & CMISS_PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  !Loop in time for STIM_STOP with the Stimulus applied.
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoop_TimesSet(ControlLoop,0.0_CMISSDP,STIM_STOP,PDE_TIME_STEP,Err)
  !Set the output
  CALL CMISSControlLoop_OutputTypeSet(ControlLoop,CMISS_CONTROL_LOOP_TIMING_OUTPUT,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)
 
  !Start the creation of the problem solvers
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  !Get the first (DAE) solver
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  !Set the DAE time step to by 10 us
  CALL CMISSSolver_DAETimeStepSet(Solver,ODE_TIME_STEP,Err)
  !CALL CMISSSolver_DAESolverTypeSet(Solver,CMISS_SOLVER_DAE_EXTERNAL,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  !Get the second (Parabolic) solver
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver CellML equations
  CALL CMISSProblem_CellMLEquationsCreateStart(Problem,Err)
  !Get the first solver  
  !Get the CellML equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSCellMLEquations_Initialise(CellMLEquations,Err)
  CALL CMISSSolver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL CMISSCellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
  !Finish the creation of the problem solver CellML equations
  CALL CMISSProblem_CellMLEquationsCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !Get the second solver  
  !Get the solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
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
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
  ELSE
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  ENDIF
  CALL CMISSDecomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL CMISSDecomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
    !  & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    !CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
    !  & CMISS_BOUNDARY_CONDITION_FIXED,1.0_CMISSDP,Err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem for the first STIM_STOP
  !CALL CMISSProblem_Solve(Problem,Err)

  !Now turn the stimulus off
  !Set the Stimulus at node 1
  DO node_idx=1,NUMBER_OF_ELEMENTS/2
    IF(FirstNodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
        & node_idx, &
        & stimcomponent,0.0_CMISSDP,Err)
    ENDIF
  ENDDO !node_idx

  !Set the time loop from STIM_STOP to TIME_STOP
  CALL CMISSControlLoop_TimesSet(ControlLoop,STIM_STOP,TIME_STOP,PDE_TIME_STEP,Err)
  
  !Solve the problem for the next 900 ms
  !CALL CMISSProblem_Solve(Problem,Err)
  
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"MonodomainExample","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"MonodomainExample","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM MONODOMAINEXAMPLE
