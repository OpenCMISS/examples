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
  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)
  
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
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Set the region label
  CALL CMISSRegionLabelSet(Region,"Region",Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasisTypeInitialise(Basis,Err)
  CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bilinear Lagrange basis
    CALL CMISSBasisNumberOfXiSet(Basis,2,Err)
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL CMISSBasisNumberOfXiSet(Basis,3,Err)
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(Basis,Err)

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)
        
  !Create the equations_set
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  !Set the equations set to be a Monodomain equations set
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetBioelectricsClass, &
    & CMISSEquationsSetMonodomainEquationType,CMISSEquationsSetNoSubtype,EquationsSetFieldUserNumber,EquationsSetField, &
    & EquationsSet,Err)
  
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)
  
  !Create the equations set materials field variables
  CALL CMISSFieldTypeInitialise(MaterialsField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)
  
  !Set Am
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,193.6_CMISSDP,Err)
  !Set Cm
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,0.014651_CMISSDP,Err)
  !Set conductivity
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,CONDUCTIVITY,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,CONDUCTIVITY,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,CONDUCTIVITY,Err)
  ENDIF

  !Create the CellML environment
  CALL CMISSCellMLTypeInitialise(CellML,Err)
  CALL CMISSCellMLCreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import a Noble 1998 model from a file
  CALL CMISSCellMLModelImport(CellML,"n98.xml",n98ModelIndex,Err)
  ! and import JRW 1998 from a file
  !CALL CMISSCellMLModelImport(CellML,"jrw-1998.xml",JRWModelIndex,Err)
  ! and import LRd from a file
  !CALL CMISSCellMLModelImport(CellML,"LRd.xml",LRdModelIndex,Err)
!  CALL CMISSDiagnosticsSetOn(CMISSInDiagType,[1,2,3,4,5],"",["CELLML_CREATE_FIELD_TO_CELLML_MAP_C", &
!    & "CELLML_CREATE_CELLML_TO_FIELD_MAP_C"],Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  !CALL CMISSCellMLVariableSetAsKnown(CellML,n98ModelIndex,"time_independent_potassium_current/g_K1",Err)
  CALL CMISSCellMLVariableSetAsKnown(CellML,n98ModelIndex,"fast_sodium_current/g_Na ",Err)
  CALL CMISSCellMLVariableSetAsKnown(CellML,n98ModelIndex,"membrane/IStim",Err)
  !CALL CMISSCellMLVariableSetAsKnown(CellML,JRWModelIndex,"L_type_Ca_channel/Ko",Err) ! this one should fail
  !CALL CMISSCellMLVariableSetAsKnown(CellML,JRWModelIndex,"membrane/I_stim",Err)
  !   - to get from the CellML side
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K1",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_to",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K_ATP",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Ca_L_K",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_b_K",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_NaK",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Na",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_b_Na",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Ca_L_Na",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_NaCa",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/IStimC",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_K1",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_Na",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_Ca_L_Ca",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_Ca_L_K",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_K",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_NaCa",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_Kp",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_p_Ca",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_Na_b",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"membrane/i_Ca_b",Err)
  !CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"membrane/IStimC",Err)
  !   - and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
!  CALL CMISSDiagnosticsSetOff(Err)
  !Finish the CellML environment
  CALL CMISSCellMLCreateFinish(CellML,Err)

!  CALL CMISSDiagnosticsSetOn(CMISSInDiagType,[1,2,3,4,5],"",["CELLML_CREATE_FIELD_TO_CELLML_MAP_C", &
!    & "CELLML_CREATE_CELLML_TO_FIELD_MAP_C"],Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellMLFieldMapsCreateStart(CellML,Err)
  !Now we can set up the field variable component <--> CellML model variable mappings.
  !Map Vm
  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType, &
    & n98ModelIndex,"membrane/V",CMISSFieldValuesSetType,Err)
  !CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType, &
  !  & JRWModelIndex,"membrane/V",CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,n98ModelIndex,"membrane/V",CMISSFieldValuesSetType, &
    & DependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType,Err)
  !CALL CMISSCellMLCreateCellMLToFieldMap(CellML,JRWModelIndex,"membrane/V",CMISSFieldValuesSetType, &
  !  & DependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType,Err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellMLFieldMapsCreateFinish(CellML,Err)

  !todo - get vm initialial value.
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,-92.5_CMISSDP,Err)
  
!  CALL CMISSDiagnosticsSetOff(Err)

  !Start the creation of the CellML models field
  CALL CMISSFieldTypeInitialise(CellMLModelsField,Err)
  CALL CMISSCellMLModelsFieldCreateStart(CellMLModelsFieldUserNumber,CellML,CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL CMISSCellMLModelsFieldCreateFinish(CellML,Err)
  !Set up the models field
  !DO N=1,(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  !  IF(N < 5) THEN
  !    CELL_TYPE = 1
  !  ELSE
  !    CELL_TYPE = 2
  !  ENDIF
  !  CALL CMISSFieldParameterSetUpdateNode(CellMLModelsField, CMISSFieldUVariableType, CMISSFieldValuesSetType,1,N,1,CELL_TYPE,Err)
  !END DO
  !CALL CMISSFieldParameterSetUpdateStart(CellMLModelsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  !CALL CMISSFieldParameterSetUpdateFinish(CellMLModelsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  !CALL CMISSFieldComponentValuesInitialise(CellMLModelsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2_CMISSIntg,Err)

  !Start the creation of the CellML state field
  CALL CMISSFieldTypeInitialise(CellMLStateField,Err)
  CALL CMISSCellMLStateFieldCreateStart(CellMLStateFieldUserNumber,CellML,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL CMISSCellMLStateFieldCreateFinish(CellML,Err)

  !Start the creation of the CellML intermediate field
  CALL CMISSFieldTypeInitialise(CellMLIntermediateField,Err)
  CALL CMISSCellMLIntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber,CellML,CellMLIntermediateField,Err)
  !Finish the creation of the CellML intermediate field
  CALL CMISSCellMLIntermediateFieldCreateFinish(CellML,Err)
  
  !Start the creation of CellML parameters field
  CALL CMISSFieldTypeInitialise(CellMLParametersField,Err)
  CALL CMISSCellMLParametersFieldCreateStart(CellMLParametersFieldUserNumber,CellML,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL CMISSCellMLParametersFieldCreateFinish(CellML,Err)
  
  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)

  CALL CMISSCellMLFieldComponentGet(CellML,n98ModelIndex,CMISSCellMLParametersFieldType,"membrane/IStim",stimcomponent,Err)
  !Set the Stimulus at half the bottom nodes
  DO node_idx=1,NUMBER_OF_ELEMENTS/2
    IF(FirstNodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetUpdateNode(CellMLParametersField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx, &
        & stimcomponent,STIM_VALUE,Err)
    ENDIF
  ENDDO

  !!Set up the g_K1 gradient
  !CALL CMISSCellMLFieldComponentGet(CellML,n98ModelIndex,CMISSCellMLParametersFieldType,"time_independent_potassium_current/g_K1", &
  !  & gK1component,Err)
  !!Loop over the nodes
  !DO node_idx=1,(NUMBER_OF_ELEMENTS+1)*(NUMBER_OF_ELEMENTS+1)
  !  CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx,1, &
  !    & X,Err)
  !  CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx,2, &
  !    & Y,Err)
  !  DISTANCE=SQRT(X**2+Y**2)/SQRT(2.0_CMISSDP)
  !  gK1_VALUE=2.0_CMISSDP*(DISTANCE+0.5_CMISSDP)*77.11e-3_CMISSDP
  !  CALL CMISSFieldParameterSetUpdateNode(CellMLParametersField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx, &
  !    & gK1component,gK1_VALUE,Err)
  !ENDDO
  
  !Set up the g_Na gradient
  CALL CMISSCellMLFieldComponentGet(CellML,n98ModelIndex,CMISSCellMLParametersFieldType,"fast_sodium_current/g_Na", &
    & gNacomponent,Err)
  !Loop over the nodes
  DO node_idx=1,(NUMBER_OF_ELEMENTS+1)*(NUMBER_OF_ELEMENTS+1)
    CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx,1, &
      & X,Err)
    CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx,2, &
      & Y,Err)
    DISTANCE=SQRT(X**2+Y**2)/SQRT(2.0_CMISSDP)
    gNa_VALUE=2.0_CMISSDP*(DISTANCE+0.5_CMISSDP)*385.5e-3_CMISSDP
    CALL CMISSFieldParameterSetUpdateNode(CellMLParametersField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx, &
      & gNacomponent,gNa_VALUE,Err)
  ENDDO
  
  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemBioelectricsClass,CMISSProblemMonodomainEquationType, &
    & CMISSProblemMonodomainGudunovSplitSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  !Loop in time for STIM_STOP with the Stimulus applied.
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,0.0_CMISSDP,STIM_STOP,PDE_TIME_STEP,Err)
  !Set the output
  CALL CMISSControlLoopOutputTypeSet(ControlLoop,CMISSControlLoopTimingOutput,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)
 
  !Start the creation of the problem solvers
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  !Get the first (DAE) solver
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  !Set the DAE time step to by 10 us
  CALL CMISSSolverDAETimeStepSet(Solver,ODE_TIME_STEP,Err)
  !CALL CMISSSolverDAESolverTypeSet(Solver,CMISSSolverDAEExternal,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
  !Get the second (Parabolic) solver
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,2,Solver,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver CellML equations
  CALL CMISSProblemCellMLEquationsCreateStart(Problem,Err)
  !Get the first solver  
  !Get the CellML equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSCellMLEquationsTypeInitialise(CellMLEquations,Err)
  CALL CMISSSolverCellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL CMISSCellMLEquationsCellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
  !Finish the creation of the problem solver CellML equations
  CALL CMISSProblemCellMLEquationsCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the second solver  
  !Get the solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,2,Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)  
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Start the creation of the equations set boundary conditions
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
  ELSE
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  ENDIF
  CALL CMISSDecompositionNodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL CMISSDecompositionNodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,FirstNodeNumber,1, &
    !  & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,LastNodeNumber,1, &
    !  & CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem for the first STIM_STOP
  !CALL CMISSProblemSolve(Problem,Err)

  !Now turn the stimulus off
  !Set the Stimulus at node 1
  DO node_idx=1,NUMBER_OF_ELEMENTS/2
    IF(FirstNodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetUpdateNode(CellMLParametersField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx, &
        & stimcomponent,0.0_CMISSDP,Err)
    ENDIF
  ENDDO !node_idx

  !Set the time loop from STIM_STOP to TIME_STOP
  CALL CMISSControlLoopTimesSet(ControlLoop,STIM_STOP,TIME_STOP,PDE_TIME_STEP,Err)
  
  !Solve the problem for the next 900 ms
  !CALL CMISSProblemSolve(Problem,Err)
  
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"MonodomainExample","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"MonodomainExample","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
  ENDIF
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM MONODOMAINEXAMPLE
