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

  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=3.0_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=17

  !Program types
  
  !Program variables
  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT,CellmlFile
  LOGICAL :: fileExist

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: n98ModelIndex

  INTEGER(CMISSIntg) :: gNacomponent,stimcomponent,node_idx

  REAL(CMISSRP) :: X,Y,DISTANCE,gNa_VALUE
  
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_ELEMENTS=25
  INTEGER(CMISSIntg) :: OUTPUT_FREQUENCY = 1
  REAL(CMISSRP), PARAMETER :: STIM_VALUE = 100.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: STIM_STOP = 0.10_CMISSRP
  REAL(CMISSRP) :: TIME_STOP = 1.50_CMISSRP
  REAL(CMISSRP), PARAMETER :: ODE_TIME_STEP = 0.00001_CMISSRP
  REAL(CMISSRP) :: PDE_TIME_STEP = 0.001_CMISSRP
  REAL(CMISSRP), PARAMETER :: CONDUCTIVITY = 0.1_CMISSRP

  !CMISS variables
  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField,MaterialsField
  TYPE(cmfe_FieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh  
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations

  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex,CellMLIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain,NodeDomain
  INTEGER(CMISSIntg) :: Err


  ! process command line arguments before getting started.
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS >= 3) THEN
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    !IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) PDE_TIME_STEP
    WRITE(*, '("PDE Step Size: ", E14.7)') PDE_TIME_STEP
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) TIME_STOP
    WRITE(*, '("Stop Time: ", E14.7)') TIME_STOP
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) OUTPUT_FREQUENCY
    WRITE(*, '("Output Frequency: ", I10)') OUTPUT_FREQUENCY
    CALL GET_COMMAND_ARGUMENT(4,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    !IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 4.")
    CellmlFile = adjustl(COMMAND_ARGUMENT)
    WRITE(*, '("CellML File: ", A)') CellmlFile
    inquire(file=CellmlFile, exist=fileExist)
    if (.not. fileExist) then
      write(*, '(">>ERROR: File does not exist")')
      stop
    endif
  ELSE
    !If there are not enough arguments die horribly
    WRITE(*,'(">>USAGE: ",A)') "MonodomainExample <PDE step size> <stop time> <output frequency> <CellML Model URL>"
    STOP
  ENDIF

  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  !Trap errors
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  
  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !CALL cmfe_OutputSetOn("Monodomain",Err)
    
  NUMBER_GLOBAL_X_ELEMENTS=NUMBER_OF_ELEMENTS
  NUMBER_GLOBAL_Y_ELEMENTS=NUMBER_OF_ELEMENTS
  NUMBER_GLOBAL_Z_ELEMENTS=0
  
  !Start the creation of a new RC coordinate system
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

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Set the region label
  CALL cmfe_Region_LabelSet(Region,"Region",Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis,2,Err)
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(Basis,3,Err)
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis,Err)

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)
  
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
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
        
  !Create the equations_set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  !Set the equations set to be a Monodomain equations set
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
    & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_NO_SUBTYPE],EquationsSetFieldUserNumber,EquationsSetField, &
    & EquationsSet,Err)
  
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)
  
  !Create the equations set materials field variables
  CALL cmfe_Field_Initialise(MaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)
  
  !Set Am
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & 193.6_CMISSRP,Err)
  !Set Cm
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
    & 0.014651_CMISSRP,Err)
  !Set conductivity
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
    & CONDUCTIVITY,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
    & CONDUCTIVITY,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5, &
      & CONDUCTIVITY,Err)
  ENDIF

  !Create the CellML environment
  CALL cmfe_CellML_Initialise(CellML,Err)
  CALL cmfe_CellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import a Noble 1998 model from a file
  CALL cmfe_CellML_ModelImport(CellML,CellmlFile,n98ModelIndex,Err)
  CALL cmfe_CellML_VariableSetAsKnown(CellML,n98ModelIndex,"fast_sodium_current/g_Na ",Err)
  CALL cmfe_CellML_VariableSetAsKnown(CellML,n98ModelIndex,"membrane/IStim",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K1",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_to",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K_ATP",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Ca_L_K",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_b_K",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_NaK",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Na",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_b_Na",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Ca_L_Na",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_NaCa",Err)
  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(CellML,Err)
  !Now we can set up the field variable component <--> CellML model variable mappings.
  !Map Vm
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & n98ModelIndex,"membrane/V",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,n98ModelIndex,"membrane/V",CMFE_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateFinish(CellML,Err)

  !todo - get vm initialial value.
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & -92.5_CMISSRP,Err)
  
  !Start the creation of the CellML models field
  CALL cmfe_Field_Initialise(CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL cmfe_CellML_ModelsFieldCreateFinish(CellML,Err)

  !Start the creation of the CellML state field
  CALL cmfe_Field_Initialise(CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL cmfe_CellML_StateFieldCreateFinish(CellML,Err)

  !Start the creation of the CellML intermediate field
  CALL cmfe_Field_Initialise(CellMLIntermediateField,Err)
  CALL cmfe_CellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  !Finish the creation of the CellML intermediate field
  CALL cmfe_CellML_IntermediateFieldCreateFinish(CellML,Err)
  
  !Start the creation of CellML parameters field
  CALL cmfe_Field_Initialise(CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL cmfe_CellML_ParametersFieldCreateFinish(CellML,Err)
  
  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Find the domains of the first and last nodes
  FirstNodeNumber=1
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
  ELSE
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  ENDIF
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  
  CALL cmfe_CellML_FieldComponentGet(CellML,n98ModelIndex,CMFE_CELLML_PARAMETERS_FIELD,"membrane/IStim",stimcomponent,Err)
  !Set the Stimulus at half the bottom nodes
  DO node_idx=1,NUMBER_OF_ELEMENTS/2
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & node_idx,stimcomponent,STIM_VALUE,Err)
    ENDIF
  ENDDO
  
  !Set up the g_Na gradient
  CALL cmfe_CellML_FieldComponentGet(CellML,n98ModelIndex,CMFE_CELLML_PARAMETERS_FIELD,"fast_sodium_current/g_Na", &
    & gNacomponent,Err)
  !Loop over the nodes
  DO node_idx=1,LastNodeNumber
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
        & X,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,2, &
        & Y,Err)
      DISTANCE=SQRT(X**2+Y**2)/SQRT(2.0_CMISSRP)
      gNa_VALUE=2.0_CMISSRP*(DISTANCE+0.5_CMISSRP)*385.5e-3_CMISSRP
      CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & node_idx,gNacomponent,gNa_VALUE,Err)
    ENDIF
  ENDDO
  
  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_BIOELECTRICS_CLASS,CMFE_PROBLEM_MONODOMAIN_EQUATION_TYPE, &
    & CMFE_PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  !Loop in time for STIM_STOP with the Stimulus applied.
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,0.0_CMISSRP,STIM_STOP,PDE_TIME_STEP,Err)
  !Set the output
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoop,CMFE_CONTROL_LOOP_TIMING_OUTPUT,Err)
  !Set the output frequency (0 for no output, n for output every n time steps)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)
 
  !Start the creation of the problem solvers
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  !Get the first (DAE) solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  !Set the DAE time step to by 10 us
  CALL cmfe_Solver_DAETimeStepSet(Solver,ODE_TIME_STEP,Err)
  !CALL cmfe_Solver_DAESolverTypeSet(Solver,CMFE_SOLVER_DAE_EXTERNAL,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  !Get the second (Parabolic) solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)
  !Get the first solver  
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
  !Finish the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the second solver  
  !Get the solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Start the creation of the equations set boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the first node to 0.0 and the last node to 1.0
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
    !  & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
    !  & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,Err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem for the first STIM_STOP
  CALL cmfe_Problem_Solve(Problem,Err)

  !Now turn the stimulus off
  !Set the Stimulus at node 1
  DO node_idx=1,NUMBER_OF_ELEMENTS/2
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & node_idx,stimcomponent,0.0_CMISSRP,Err)
    ENDIF
  ENDDO !node_idx

  !Set the time loop from STIM_STOP to TIME_STOP
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,STIM_STOP,TIME_STOP,PDE_TIME_STEP,Err)
  
  !Solve the problem for the next 900 ms
  CALL cmfe_Problem_Solve(Problem,Err)
  
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"MonodomainExample","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"MonodomainExample","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF
  
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM MONODOMAINEXAMPLE
