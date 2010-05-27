!> \file
!> $Id: DiffusionExample.f90 20 2007-05-28 20:22:52Z cpb $
!> \author Chris Bradley
!> \brief This is an example program to solve a diffusion equation using openCMISS calls.
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

!> \example ClassicalField/AdvectionDiffusion/StaticAdvectionDiffusion/src/StaticAdvectionDiffusionExample.f90
!! Example program to solve a diffusion equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/Diffusion/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/Diffusion/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM STATICADVECTIONDIFFUSIONEXAMPLE

  USE OPENCMISS
  USE FIELDML_OUTPUT_ROUTINES
  USE FIELDML_UTIL_ROUTINES
  USE FIELDML_API
  USE MPI


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=3.0_CMISSDP 
  REAL(CMISSDP), POINTER :: GEOMETRIC_PARAMETERS(:)
  
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopNode=0
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=14

  INTEGER(CMISSIntg), PARAMETER :: MeshComponentNumber=1

  CHARACTER(C_CHAR), PARAMETER :: NUL=C_NULL_CHAR

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  
  INTEGER(CMISSIntg) :: MPI_IERROR
  
    !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField,IndependentField,AnalyticField,SourceField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver, LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

  LOGICAL :: EXPORT_FIELD,IMPORT_FIELD
 

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: Err
  
  !FieldML variables
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: filename = "StaticAdvectionDiffusion"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: basename = "static_advection_diffusion"
  TYPE(C_PTR) :: fmlHandle, writer
  INTEGER(C_INT) :: meshHandle, nodesHandle, elementsHandle, nodalDofsHandle, real1DHandle, importHandle
  INTEGER(C_INT) :: layoutHandle, connectivityHandle, femHandle, templateHandle, fieldHandle, tempHandle
  INTEGER(C_INT) :: geometricDofsHandle, dependentDofsHandle
  
  INTEGER(CMISSIntg) :: i, j, count1, count2, coordsType, basisNodeCount
  INTEGER(CMISSIntg) :: nodeCount, elementCount, dimensions, componentCount
  INTEGER(C_INT), ALLOCATABLE, TARGET :: iBuffer(:)
  INTEGER(C_INT), ALLOCATABLE, TARGET :: componentHandles(:)
  REAL(C_DOUBLE), ALLOCATABLE, TARGET :: dBuffer(:)
  REAL(C_DOUBLE) :: dValue
  INTEGER(C_INT), TARGET :: dummy(0)

  
  TYPE(CMISSMeshElementsType) :: meshElements
  
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

  NUMBER_GLOBAL_X_ELEMENTS=50
  NUMBER_GLOBAL_Y_ELEMENTS=100
  NUMBER_GLOBAL_Z_ELEMENTS=0
  NUMBER_OF_DOMAINS=1


  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  
  IF( NUMBER_GLOBAL_Z_ELEMENTS == 0 ) THEN
    dimensions = 2
  ELSE
    dimensions = 3
  ENDIF

    !Start the creation of a new RC coordinate system
    CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
    CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
    CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,dimensions,Err)
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
    CALL CMISSBasisNumberOfXiSet(Basis,dimensions,Err)
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(BASIS,Err)

    !Start the creation of a generated mesh in the region
    CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
    CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
    !Set up a regular x*y*z mesh
    CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
    !Set the default basis
    CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)   
    !Define the mesh on the region
    IF(dimensions == 2) THEN
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
  DO i = 1, dimensions
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,i,MeshComponentNumber,Err)
  ENDDO
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)

       
    !Update the geometric field parameters
    CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)
!  ENDIF

  !IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD)) GEOMETRIC_FIELD=>REGION%FIELDS%FIELDS(1)%PTR
  
  !Create the equations_set
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,EquationsSet,Err)
  !Set the equations set to be a standard Laplace problem
!   CALL CMISSEquationsSetSpecificationSet(EquationsSet,CMISSEquationsSetClassicalFieldClass, &
!     & CMISSEquationsSetAdvectionDiffusionEquationType,CMISSEquationsSetNoSourceStaticAdvecDiffSubtype,Err)
  CALL CMISSEquationsSetSpecificationSet(EquationsSet,CMISSEquationsSetClassicalFieldClass, &
    & CMISSEquationsSetAdvectionDiffusionEquationType,CMISSEquationsSetConstantSourceStaticAdvecDiffSubtype,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field variables
  CALL CMISSFieldTypeInitialise(MaterialsField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set source field variables
  !For comparison withe analytical solution used here, the source field must be set to the following:
  !f(x,y) = 2.0*tanh(-0.1E1+Alpha*(TanPhi*x-y))*(1.0-pow(tanh(-0.1E1+Alpha*(TanPhi*x-y)),2.0))*Alpha*Alpha*TanPhi*TanPhi
  !+2.0*tanh(-0.1E1+Alpha*(TanPhi*x-y))*(1.0-pow(tanh(-0.1E1+Alpha*(TanPhi*x-y)),2.0))*Alpha*Alpha
  !-Peclet*(-sin(6.0*y)*(1.0-pow(tanh(-0.1E1+Alpha*(TanPhi*x-y)),2.0))*Alpha*TanPhi+cos(6.0*x)*(1.0-pow(tanh(-0.1E1+Alpha*(TanPhi*x-y)),2.0))*Alpha)
  CALL CMISSFieldTypeInitialise(SourceField,Err)
  CALL CMISSEquationsSetSourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  CALL CMISSFieldComponentInterpolationSet(SourceField,CMISSFieldUVariableType,1,CMISSFieldNodeBasedInterpolation,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetSourceCreateFinish(EquationsSet,Err)

  CALL CMISSFieldParameterSetDataGet(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,GEOMETRIC_PARAMETERS,Err)
  !Create the equations set independent field variables
  CALL CMISSFieldTypeInitialise(IndependentField,Err)
  CALL CMISSEquationsSetIndependentCreateStart(EquationsSet,IndependentFieldUserNumber,IndependentField,Err)
  IF( dimensions == 2 ) THEN
  !For comparison withe analytical solution used here, the independent field must be set to the following:
  !w(x,y)=(sin 6y,cos 6x) FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION
!   CALL CMISSFieldComponentInterpolationSet(IndependentField,CMISSFieldUVariableType,1,CMISSFieldNodeBasedInterpolation,Err) 
!   CALL CMISSFieldComponentInterpolationSet(IndependentField,CMISSFieldUVariableType,2,CMISSFieldNodeBasedInterpolation,Err)
  !Loop over nodes to set the appropriate function value
!    DO

!    ENDDO
  ENDIF
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetIndependentCreateFinish(EquationsSet,Err)

  !Create the equations set analytic field variables
  CALL CMISSFieldTypeInitialise(AnalyticField,Err)
  IF( dimensions == 2) THEN  
    CALL CMISSEquationsSetAnalyticCreateStart(EquationsSet,CMISSEquationsSetAdvectionDiffusionTwoDim1,&
      & AnalyticFieldUserNumber,AnalyticField,Err)
  ELSE
    WRITE(*,'(A)') "Three dimensions is not implemented."
    STOP
  ENDIF
  !Finish the equations set analytic field variables
  CALL CMISSEquationsSetAnalyticCreateFinish(EquationsSet,Err)

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
  
!   !Create the equations set boundary conditions
!   CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
!   CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSet,BoundaryConditions,Err)
!   !Set the first node to 0.0 and the last node to 1.0
!   FirstNodeNumber=1
!   IF( dimensions == 2 ) THEN
!     LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
!   ELSE
!     LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
!   ENDIF
!   CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,FirstNodeNumber,1, &
!     & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
!   CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldDeludelnVariableType,1,LastNodeNumber,1, &
!     & CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
!   !Finish the creation of the equations set boundary conditions
!   CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSet,Err)
CALL CMISSEquationsSetBoundaryConditionsAnalytic(EquationsSet,Err)

!   EXPORT_FIELD=.TRUE.
!   IF(EXPORT_FIELD) THEN
!     CALL CMISSFieldsTypeInitialise(Fields,Err)
!     CALL CMISSFieldsTypeCreate(Region,Fields,Err)
!     CALL CMISSFieldIONodesExport(Fields,"StaticAdvectionDiffusionInitial","FORTRAN",Err)
!     CALL CMISSFieldIOElementsExport(Fields,"StaticAdvectionDiffusionInitial","FORTRAN",Err)
!     CALL CMISSFieldsTypeFinalise(Fields,Err)
! 
!   ENDIF


  !Create the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a no source static advection Diffusion problem
!   CALL CMISSProblemSpecificationSet(Problem,CMISSProblemClassicalFieldClass,CMISSProblemAdvectionDiffusionEquationType, &
!     & CMISSProblemNoSourceStaticAdvecDiffSubtype,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemClassicalFieldClass,CMISSProblemAdvectionDiffusionEquationType, &
    & CMISSProblemLinearSourceStaticAdvecDiffSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)


  !Create the problem control
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  !CALL CMISSProblemControlLoopGet(Problem,ControlLoopNode,ControlLoop,Err)
  !Set the times
  !CALL CMISSControlLoopTimesSet(ControlLoop,0.0_CMISSDP,1.0_CMISSDP,0.1_CMISSDP,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)


  !Start the creation of the problem solvers
!  
! !   !For the Direct Solver MUMPS, uncomment the below two lines and comment out the above five
! !   CALL SOLVER_LINEAR_TYPE_SET(LINEAR_SOLVER,SOLVER_LINEAR_DIRECT_SOLVE_TYPE,ERR,ERROR,*999)
! !   CALL SOLVER_LINEAR_DIRECT_TYPE_SET(LINEAR_SOLVER,SOLVER_DIRECT_MUMPS,ERR,ERROR,*999) 
! 

  CALL CMISSSolverTypeInitialise(Solver,Err)
  !CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverDynamicLinearSolverGet(Solver,LinearSolver,Err)
  !CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolver,300,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)


  !Create the problem solver equations
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

  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)

 !Output Analytic analysis
  Call CMISSAnalyticAnalysisOutput(DependentField,"StaticAdvectionDiffusionAnalytics",Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"StaticAdvectionDiffusion","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"StaticAdvectionDiffusion","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
    
    fmlHandle = Fieldml_Create()
    
    CALL CMISSNumberOfNodesGet( Region, nodeCount, err )
    
    nodesHandle = Fieldml_CreateEnsembleDomain( fmlHandle, "advection.nodes"//NUL, FML_INVALID_HANDLE )
    err = Fieldml_SetContiguousBoundsCount( fmlHandle, nodesHandle, nodeCount )
    err = Fieldml_SetMarkup( fmlHandle, nodesHandle, "geometric"//NUL, "point"//NUL )
    
    CALL CMISSMeshNumberOfElementsGet( Mesh, elementCount, err )

    CALL FieldmlUtil_GetXiEnsemble( fmlHandle, dimensions, tempHandle, err )
    meshHandle = Fieldml_CreateMeshDomain( fmlHandle, "advection.mesh"//NUL, tempHandle )
    err = Fieldml_SetContiguousBoundsCount( fmlHandle, meshHandle, elementCount )

    elementsHandle = Fieldml_GetMeshElementDomain( fmlHandle, meshHandle )
    
    !CALL FieldmlUtil_GetConnectivityEnsemble( fmlHandle, BasisUserNumber, layoutHandle, err )
    
    CALL FieldmlUtil_GetGenericDomain( fmlHandle, 1, real1DHandle, err )
    nodalDofsHandle = Fieldml_CreateContinuousVariable( fmlHandle, "advection.dofs.nodal"//NUL, real1DHandle )

    CALL CMISSMeshNumberOfComponentsGet( mesh, count1, err )
    ALLOCATE( componentHandles( count1 ) )
    DO i = 1, count1
      CALL CMISSMeshElementsGet( mesh, i, meshElements, err )
      CALL FieldmlOutput_AddConnectivity( fmlHandle, baseName, meshHandle, nodesHandle, nodalDofsHandle, i, meshElements, &
        & componentHandles( i ), err )
    ENDDO
    
    IF( dimensions == 2 ) THEN
      err = Fieldml_SetMeshDefaultShape( fmlHandle, meshHandle, "library.shape.square"//NUL )
    ELSE
    ENDIF
    
    CALL CMISSCoordinateSystemTypeGet( CoordinateSystem, coordsType, err )
    CALL CMISSFieldNumberOfComponentsGet( GeometricField, CMISSFieldUVariableType, count1, err )
    CALL FieldmlUtil_GetCoordinatesDomain( fmlHandle, coordsType, count1, tempHandle, err )

    IF( tempHandle /= FML_INVALID_HANDLE ) THEN

      fieldHandle = Fieldml_CreateContinuousAggregate( fmlHandle, "advection.geometric"//NUL, tempHandle )
      DO i = 1, count1
        CALL CMISSFieldComponentMeshComponentGet(GeometricField, CMISSFieldUVariableType,i,j,Err)
        err = Fieldml_SetEvaluator( fmlHandle, fieldHandle, i, componentHandles(j) )
      ENDDO
      geometricDofsHandle = Fieldml_CreateContinuousParameters( fmlHandle, "advection.dofs.geometric.nodal"//NUL, real1DHandle )
      tempHandle = Fieldml_GetDomainComponentEnsemble( fmlHandle, tempHandle )
      err = Fieldml_SetParameterDataDescription( fmlHandle, geometricDofsHandle, DESCRIPTION_SEMIDENSE )
      err = Fieldml_SetParameterDataLocation( fmlHandle, geometricDofsHandle, LOCATION_FILE )
      err = Fieldml_SetParameterFileData( fmlHandle, geometricDofsHandle, "advection.geometric"//NUL, TYPE_LINES, 0 )

      IF( tempHandle /= FML_INVALID_HANDLE ) THEN
        err = Fieldml_AddSemidenseIndex( fmlHandle, geometricDofsHandle, tempHandle, 0 )
      ENDIF
      err = Fieldml_AddSemidenseIndex( fmlHandle, geometricDofsHandle, nodesHandle, 0 )
      err = Fieldml_SetAlias( fmlHandle, fieldHandle, nodalDofsHandle, geometricDofsHandle )

      ALLOCATE( dBuffer( count1 ) )
      writer = Fieldml_OpenWriter( fmlHandle, geometricDofsHandle, 0 )
      DO i = 1, nodeCount
        DO j = 1, count1
          CALL CMISSFieldParameterSetGetNode( GeometricField, CMISSFieldUVariableType, CMISSFieldValuesSetType, & 
            & CMISSNoGlobalDerivative, i, j, dValue, err )
          dBuffer( j ) = dValue
        ENDDO
        err = Fieldml_WriteDoubleSlice( fmlHandle, writer, C_LOC(dummy), C_LOC(dBuffer) )
      ENDDO
      err = Fieldml_CloseWriter( fmlHandle, writer )
      DEALLOCATE( dBuffer )
    ENDIF


    CALL CMISSFieldNumberOfComponentsGet( DependentField, CMISSFieldUVariableType, count1, err )
    CALL FieldmlUtil_GetGenericDomain( fmlHandle, count1, tempHandle, err )

    IF( tempHandle /= FML_INVALID_HANDLE ) THEN
      fieldHandle = Fieldml_CreateContinuousAggregate( fmlHandle, "advection.dependent"//NUL, tempHandle )
      DO i = 1, count1
        CALL CMISSFieldComponentMeshComponentGet(DependentField, CMISSFieldUVariableType,i,j,Err)
        err = Fieldml_SetEvaluator( fmlHandle, fieldHandle, i, componentHandles(j) )
      ENDDO
      dependentDofsHandle = Fieldml_CreateContinuousParameters( fmlHandle, "advection.dofs.dependent.nodal"//NUL, real1DHandle )
      tempHandle = Fieldml_GetDomainComponentEnsemble( fmlHandle, tempHandle )
      err = Fieldml_SetParameterDataDescription( fmlHandle, dependentDofsHandle, DESCRIPTION_SEMIDENSE )
      err = Fieldml_SetParameterDataLocation( fmlHandle, dependentDofsHandle, LOCATION_FILE )
      err = Fieldml_SetParameterFileData( fmlHandle, dependentDofsHandle, "advection.dependent"//NUL, TYPE_LINES, 0 )

      IF( tempHandle /= FML_INVALID_HANDLE ) THEN
        err = Fieldml_AddSemidenseIndex( fmlHandle, dependentDofsHandle, tempHandle, 0 )
      ENDIF
      err = Fieldml_AddSemidenseIndex( fmlHandle, dependentDofsHandle, nodesHandle, 0 )
      err = Fieldml_SetAlias( fmlHandle, fieldHandle, nodalDofsHandle, dependentDofsHandle )

      ALLOCATE( dBuffer( count1 ) )
      writer = Fieldml_OpenWriter( fmlHandle, dependentDofsHandle, 0 )
      DO i = 1, nodeCount
        DO j = 1, count1
          CALL CMISSFieldParameterSetGetNode( DependentField, CMISSFieldUVariableType, CMISSFieldValuesSetType, & 
            & CMISSNoGlobalDerivative, i, j, dValue, err )
          dBuffer( j ) = dValue
        ENDDO
        err = Fieldml_WriteDoubleSlice( fmlHandle, writer, C_LOC(dummy), C_LOC(dBuffer) )
      ENDDO
      err = Fieldml_CloseWriter( fmlHandle, writer )
      DEALLOCATE( dBuffer )
    ENDIF

    err = Fieldml_WriteFile( fmlHandle, filename//".xml"//NUL )

    err = Fieldml_Destroy( fmlHandle )

  ENDIF

  !CALL CMISSFinalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM STATICADVECTIONDIFFUSIONEXAMPLE
