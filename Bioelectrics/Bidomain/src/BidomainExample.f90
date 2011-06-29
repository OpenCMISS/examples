!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a Bidomain equation using OpenCMISS calls.
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

!> \example Bioelectrics/Bidomain/src/BidomainExample.f90
!! Example program to solve a Bidomain equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/Bioelectrics/Bidomain/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/Bioelectrics/Bidomain/build-gnu'>Linux GNU Build</a>
!!
!<

!> Main program
PROGRAM BIDOMAINEXAMPLE

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
  REAL(CMISSDP), PARAMETER :: WIDTH=2.0_CMISSDP
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
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: ParabolicEquationsSetUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: EllipticEquationsSetUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=18

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  
  INTEGER(CMISSIntg) :: MPI_IERROR

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: N,CELL_TYPE

  INTEGER(CMISSIntg) :: n98ModelIndex,JRWModelIndex

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCellMLType) :: CellML
  TYPE(CMISSCellMLEquationsType) :: CellMLEquations
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: ParabolicEquationsSet,EllipticEquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField,SourceField
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
  INTEGER(CMISSIntg) :: ParabolicEquationsSetIndex,EllipticEquationsSetIndex,CellMLIndex
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

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)
    
  NUMBER_GLOBAL_X_ELEMENTS=2
  NUMBER_GLOBAL_Y_ELEMENTS=2
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
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/WIDTH,HEIGHT/),Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/),Err)
  ELSE
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/WIDTH,HEIGHT,LENGTH/),Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS/),Err)
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
        
  !Create the first (parabolic) equations set for the bidomain equations
  CALL CMISSEquationsSetTypeInitialise(ParabolicEquationsSet,Err)
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(ParabolicEquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetBioelectricsClass, &
    & CMISSEquationsSetBidomainEquationType,CMISSEquationsSetFirstBidomainSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,ParabolicEquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(ParabolicEquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(ParabolicEquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(ParabolicEquationsSet,Err)
  
  !Create the equations set materials field variables
  CALL CMISSFieldTypeInitialise(MaterialsField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(ParabolicEquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(ParabolicEquationsSet,Err)
  
  !Create the equations set source field variables
  CALL CMISSFieldTypeInitialise(SourceField,Err)
  CALL CMISSEquationsSetSourceCreateStart(ParabolicEquationsSet,SourceFieldUserNumber,SourceField,Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSetSourceCreateFinish(ParabolicEquationsSet,Err)

  !Create the parabolic equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(ParabolicEquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(ParabolicEquationsSet,Err)

  !Create the second (elliptic) equations set for the bidomain equations
  CALL CMISSEquationsSetTypeInitialise(EllipticEquationsSet,Err)
  CALL CMISSEquationsSetCreateStart(EllipticEquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetBioelectricsClass, &
    & CMISSEquationsSetBidomainEquationType,CMISSEquationsSetSecondBidomainSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,EllipticEquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EllipticEquationsSet,Err)

  !Create the equations set dependent field variables using the field from the first (parabolic) equations set
  CALL CMISSEquationsSetDependentCreateStart(EllipticEquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EllipticEquationsSet,Err)
  
  !Create the equations set materials field variables using the field from the first (parabolic) equations set
  CALL CMISSEquationsSetMaterialsCreateStart(EllipticEquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EllipticEquationsSet,Err)
  
  !Create the equations set source field variables using the field from the first (parabolic) equations set
  CALL CMISSEquationsSetSourceCreateStart(EllipticEquationsSet,SourceFieldUserNumber,SourceField,Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSetSourceCreateFinish(EllipticEquationsSet,Err)

  !Create the elliptic equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EllipticEquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EllipticEquationsSet,Err)

  !Create the CellML environment for the source field
  CALL CMISSCellMLTypeInitialise(CellML,Err)
  CALL CMISSCellMLCreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import a Noble 1998 model from a file
  CALL CMISSCellMLModelImport(CellML,"n98.xml",n98ModelIndex,Err)
  ! and import JRW 1998 from a file
  CALL CMISSCellMLModelImport(CellML,"jrw-1998.xml",JRWModelIndex,Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL CMISSCellMLVariableSetAsKnown(CellML,n98ModelIndex,"ionic_concentrations/K_o",Err)
  CALL CMISSCellMLVariableSetAsKnown(CellML,n98ModelIndex,"environment/IStim",Err)
  CALL CMISSCellMLVariableSetAsKnown(CellML,JRWModelIndex,"L_type_Ca_channel/Ko",Err)
  CALL CMISSCellMLVariableSetAsKnown(CellML,JRWModelIndex,"environment/I_stim",Err)
  !   - to get from the CellML side
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"environment/i_K1",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"environment/IStimC",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"environment/i_K1",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,JRWModelIndex,"environment/IStimC",Err)
  !Finish the CellML environment
  CALL CMISSCellMLCreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellMLFieldMapsCreateStart(CellML,Err)
  !Now  set up the field variable component <--> CellML model variable mappings.
  !Map Vm
  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType, &
    & n98ModelIndex,"Vm",CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType, &
    & JRWModelIndex,"Vm",CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,n98ModelIndex,"Vm",CMISSFieldValuesSetType, &
    & DependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,JRWModelIndex,"Vm",CMISSFieldValuesSetType, &
    & DependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType,Err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellMLFieldMapsCreateFinish(CellML,Err)

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
  
  !Generate the CellML
  CALL CMISSCellMLGenerate(CellML,Err)

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemBioelectricsClass,CMISSProblemBidomainEquationType, &
    & CMISSProblemBidomainGudunovSplitSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,0.0_CMISSDP,1.01_CMISSDP,0.1_CMISSDP,Err)
  !Set the output type
  !CALL CMISSControlLoopOutputTypeSet(ControlLoop,CMISSControlLoopNoOutput,Err)
  CALL CMISSControlLoopOutputTypeSet(ControlLoop,CMISSControlLoopProgressOutput,Err)
  !CALL CMISSControlLoopOutputTypeSet(ControlLoop,CMISSControlLoopTimingOutput,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)
 
  !Start the creation of the problem solvers
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  !Get the first (DAE) solver
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !Get the second (Parabolic) solver
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,2,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
  !Get the third (Elliptic) solver
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,3,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
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
  !Get the second (parabolic) solver  
  !Get the solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,2,Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)  
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,ParabolicEquationsSet,ParabolicEquationsSetIndex,Err)
  !Get the third (elliptic) solver  
  !Get the solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,3,Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)  
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EllipticEquationsSet,EllipticEquationsSetIndex,Err)
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
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,FirstNodeNumber,1, &
      & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,LastNodeNumber,1, &
      & CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  ENDIF

  !Set boundary conditions for the elliptic equations
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
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldVVariableType,1,1,FirstNodeNumber,1, &
      & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldVVariableType,1,1,LastNodeNumber,1, &
      & CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  ENDIF
  !Finish the creation of the elliptic equations set boundary conditions
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)

  EXPORT_FIELD=.FALSE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"BidomainExample","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"BidomainExample","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
  ENDIF
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM BIDOMAINEXAMPLE
