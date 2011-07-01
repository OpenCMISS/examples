!> \file
!> \author Ishani Roy & Sander Land
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

!> \example MonodomainBuenoOrovioExample.f90
!! Example program to solve a Monodomain equation using OpenCMISS calls.
!! \htmlinclude Monodomain/MonodomainBuenoOrovio/history.html
!!
!<

!> Main program
PROGRAM MONODOMAINBUENOOROVIOEXAMPLE

  USE OPENCMISS
  USE MPI


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField


  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=7.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=20.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=3.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=11 
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumber=12

  !Program types
  
  !Program variables
  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, N, D, NA
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  
  INTEGER(CMISSIntg) :: MPI_IERROR

  LOGICAL :: EXPORT_FIELD


  REAL(CMISSDP) :: START_TIME = 0.0, END_TIME = 100.0, DT = 0.1, DX=0.2, ACTIV_R = 1.5+1e-6   ! ms ms ms mm mm
  REAL(CMISSDP), PARAMETER  :: FiberD = 0.095298372513562, TransverseD = 0.0125758411473; ! TODO CORRECT

  REAL(CMISSDP) :: x,y,z,activ, activmax

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField,IndependentField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver, LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSControlLoopType) :: ControlLoop

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
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
  IF(NUMBER_OF_ARGUMENTS >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*), DX
  END IF
  IF(NUMBER_OF_ARGUMENTS >= 2) THEN
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*), DT
  END IF
  IF(NUMBER_OF_ARGUMENTS >= 3) THEN
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*), END_TIME
  END IF
 

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  !CALL CMISSDiagnosticsSetOn(CMISSAllDiagType,(/1,2,3,4,5/),"Diagnostics",(/"FIELD_MAPPINGS_CALCULATE"/),Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  
  NUMBER_GLOBAL_X_ELEMENTS= ROUND(WIDTH / DX)
  NUMBER_GLOBAL_Y_ELEMENTS= ROUND(HEIGHT / DX)

  NUMBER_GLOBAL_Z_ELEMENTS= ROUND(LENGTH / DX)


  WRITE(*,*) '[',ComputationalNodeNumber,'] Solving ', NUMBER_GLOBAL_X_ELEMENTS*NUMBER_GLOBAL_Y_ELEMENTS*&
             & MAX(NUMBER_GLOBAL_Z_ELEMENTS,1),' elements on ', NumberOfComputationalNodes,' processors with DT=',DT,&
             & 'for T from 0 to ',END_TIME
  NUMBER_OF_DOMAINS=NumberOfComputationalNodes
    
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

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
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
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
  
  !Create the equations_set
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,CMISSProblemBioelectricsClass, &
    & CMISSEquationsSetMonodomainSSEquationType,CMISSEquationsSetMonodomainBuenoOrovioSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  CALL CMISSFieldTypeInitialise(MaterialsField,Err)

  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)
  ! activation params: cube near (0,0,0) of size 1.5 mm 
  DO N=1,(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
    CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,N,1,x,Err)
    CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,N,2,y,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,N,3,z,Err)
    ELSE
      z = 0
    END IF
    IF(x <= ACTIV_R .AND. y < ACTIV_R .AND. z < ACTIV_R) THEN 
      activ = 1.0 ! double current. default not enough?
    ELSE
      activ = 0.0
    END IF
    !WRITE(*,*) N,'@',x,y,z, '->',activ
    CALL CMISSFieldParameterSetUpdateNode(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,N,1,activ,Err)
  END DO
  ! diffusion coefficient, related to conductivity by D = sigma/Xi Cm
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,FiberD,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,TransverseD,Err)

  ! fiber direction unit vector
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,1.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,0.0_CMISSDP,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,6,0.0_CMISSDP,Err)
  END IF


  ! create/init independent field
  CALL CMISSFieldTypeInitialise(IndependentField,Err)
  CALL CMISSEquationsSetIndependentCreateStart(EquationsSet,IndependentFieldUserNumber,IndependentField,Err)
  CALL CMISSEquationsSetIndependentCreateFinish(EquationsSet,Err)

  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)
  
  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a standard Monodomain problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemBioelectricsClass, &
    & CMISSProblemMonodomainStrangSplittingEquationType,CMISSProblemMonodomainBuenoOrovioSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
   CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
   CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
   CALL CMISSControlLoopTimesSet(ControlLoop, START_TIME, END_TIME, DT, Err)! set begin/end timings dt = 0.1
   CALL CMISSControlLoopTimeOutputSet(ControlLoop,1,Err)    !Set the output timing

  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)
 
  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
 ! CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
  !!! In standard Monodomain
  !CALL CMISSSolverLinearTypeSet(Solver,CMISSSolverLinearDirectSolveType,Err)
  !CALL CMISSSolverLibraryTypeSet(Solver,CMISSSolverMUMPSLibrary,Err)
  ! In generalised
  CALL CMISSSolverDynamicLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolver,100,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
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
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)

  EXPORT_FIELD=.FALSE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"MonodomainBO","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"MonodomainBO","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
  ENDIF
  
  !Read activation times
  activmax = 0.0
  NA = 0
  DO N=1,(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
   CALL CMISSDecompositionNodeDomainGet(Decomposition,N,1,D,Err)
   IF(D==ComputationalNodeNumber) THEN
    CALL CMISSFieldParameterSetGetNode(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,N,7,activ,Err)
!    WRITE(*,*) 'Node ',N,' activation ',activ
    IF(activ > 0.0) THEN
      NA = NA+1
    ENDIF
    activmax = max(activ,activmax)
   ENDIF
  ENDDO

  WRITE(*,*) 100.0*NA/((NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)),&
             & '% complete, last activation=',activmax

  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

contains

INTEGER(CMISSIntg) FUNCTION ROUND(X)
  REAL(CMISSDP), INTENT(IN) :: X
  ROUND = INT(X+0.5)
  RETURN
END FUNCTION ROUND

END PROGRAM MONODOMAINBUENOOROVIOEXAMPLE
