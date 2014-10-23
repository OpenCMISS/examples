!> \file
!> \author David Nickerson
!> \brief This is an example program to demonstrate the use of CellML models in OpenCMISS using the Fortran API.
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

!> \example cellml/Fortran/src/FortranExample.f90
!! Example program to demonstrate the use of CellML models in OpenCMISS using the Fortran API.
!! \htmlinclude cellml/Fortran/history.html
!<

!> Main program
PROGRAM CELLMLFORTRANEXAMPLE

  USE OPENCMISS
  USE MPI

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

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  CHARACTER(LEN=255) :: ModelUrl
  
  INTEGER(CMISSIntg) :: modelIndex,MPI_IERROR
  INTEGER(CMISSIntg) :: gNacomponent,stimcomponent,node_idx

  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_ELEMENTS=2

  REAL(CMISSDP) :: X,Y,DISTANCE,gNa_VALUE
  
  REAL(CMISSDP), PARAMETER :: STIM_VALUE = 100.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: STIM_STOP = 0.10_CMISSDP
  REAL(CMISSDP), PARAMETER :: TIME_STOP = 1.50_CMISSDP
  REAL(CMISSDP), PARAMETER :: CONDUCTIVITY = 0.1_CMISSDP

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField
  !DOC-START cellml type declarations
  TYPE(CMISSCellMLType) :: CellML
  TYPE(CMISSFieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  !DOC-END cellml type declarations
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSRegionType) :: Region,WorldRegion

  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: nodeDomain
  INTEGER(CMISSIntg) :: Err

  ! we want to get the CellML model from the command line
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,ModelUrl,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) WRITE(*,'(">>ERROR: Error for command argument 1.")')
  ELSE
    ModelUrl="n98.xml"
  ENDIF

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Trap errors
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
  
  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  CALL CMISSOutputSetOn("CellMLFortran",Err)
    
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
  !Import the specified model into the CellML environment
  CALL CMISSCellML_ModelImport(CellML,ModelUrl,modelIndex,Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL CMISSCellML_VariableSetAsKnown(CellML,modelIndex,"fast_sodium_current/g_Na ",Err)
  CALL CMISSCellML_VariableSetAsKnown(CellML,modelIndex,"membrane/IStim",Err)
  !   - to get from the CellML side (state variables are wanted by default and can not be changed)
  CALL CMISSCellML_VariableSetAsWanted(CellML,modelIndex,"membrane/i_K1",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,modelIndex,"membrane/i_to",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,modelIndex,"membrane/i_K",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,modelIndex,"membrane/i_K_ATP",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,modelIndex,"membrane/i_Ca_L_K",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,modelIndex,"membrane/i_b_K",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,modelIndex,"membrane/i_NaK",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,modelIndex,"membrane/i_Na",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,modelIndex,"membrane/i_b_Na",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,modelIndex,"membrane/i_Ca_L_Na",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,modelIndex,"membrane/i_NaCa",Err)
  !   - and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
!  CALL CMISSDiagnosticsSetOff(Err)
  !Finish the CellML environment
  CALL CMISSCellML_CreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellML_FieldMapsCreateStart(CellML,Err)
  !Now we can set up the field variable component <--> CellML model variable mappings.
  !Map Vm - field to CellML
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & modelIndex,"membrane/V",CMISS_FIELD_VALUES_SET_TYPE,Err)
  ! - and CellML to field
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,modelIndex,"membrane/V",CMISS_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellML_FieldMapsCreateFinish(CellML,Err)

  !todo - get vm initial value.
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & -92.5_CMISSDP, &
    & Err)
  
  !Start the creation of the CellML models field
  CALL CMISSField_Initialise(CellMLModelsField,Err)
  CALL CMISSCellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL CMISSCellML_ModelsFieldCreateFinish(CellML,Err)

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
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)

  CALL CMISSCellML_FieldComponentGet(CellML,modelIndex,CMISS_CELLML_PARAMETERS_FIELD,"membrane/IStim",stimcomponent,Err)
  !Set the Stimulus at half the bottom nodes
  DO node_idx=1,NUMBER_OF_ELEMENTS/2
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node_idx,1,nodeDomain,Err)
    IF(nodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1, &
        & node_idx, &
        & stimcomponent,STIM_VALUE,Err)
    ENDIF
  ENDDO

  !Set up the g_Na gradient
  CALL CMISSCellML_FieldComponentGet(CellML,modelIndex,CMISS_CELLML_PARAMETERS_FIELD,"fast_sodium_current/g_Na", &
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
    
  !Finalise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM CELLMLFORTRANEXAMPLE
