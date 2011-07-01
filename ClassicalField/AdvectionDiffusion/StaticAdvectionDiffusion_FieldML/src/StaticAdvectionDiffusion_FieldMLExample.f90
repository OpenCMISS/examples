!> \file
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

!> \example ClassicalField/AdvectionDiffusion/StaticAdvectionDiffusion_FieldML/src/StaticAdvectionDiffusion_FieldMLExample.f90
!! Example program to solve a diffusion equation using openCMISS calls.
!!
!! \htmlinclude ClassicalField/AdvectionDiffusion/StaticAdvectionDiffusion_FieldML/history.html
!<

!> Main program
PROGRAM STATICADVECTIONDIFFUSIONEXAMPLE

  USE OPENCMISS
  USE FIELDML_TYPES
  USE FIELDML_API
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

  INTEGER(CMISSIntg) :: dimensions, i
  
  !FieldML variables
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputDirectory = ""
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputFilename = "StaticAdvectionDiffusion.xml"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: basename = "static_advection_diffusion"

  TYPE(FieldmlInfoType) :: fieldmlInfo

  
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
    CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetClassicalFieldClass, &
    & CMISSEquationsSetAdvectionDiffusionEquationType,CMISSEquationsSetConstantSourceStaticAdvecDiffSubtype,&
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  !Set the equations set to be a standard Laplace problem
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
  !CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
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
CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
CALL CMISSProblemSolverEquationsBoundaryConditionsAnalytic(SolverEquations,Err)
CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)

 !Output Analytic analysis
  Call CMISSAnalyticAnalysisOutput(DependentField,"StaticAdvectionDiffusionAnalytics",Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    !CALL CMISSFieldsTypeInitialise(Fields,Err)
    !CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    !CALL CMISSFieldIONodesExport(Fields,"StaticAdvectionDiffusion","FORTRAN",Err)
    !CALL CMISSFieldIOElementsExport(Fields,"StaticAdvectionDiffusion","FORTRAN",Err)
    !CALL CMISSFieldsTypeFinalise(Fields,Err)
    
    CALL CMISSFieldmlOutput_InitialiseInfo( Region, Mesh, dimensions, outputDirectory, basename, fieldmlInfo, err )

    CALL CMISSFieldmlOutput_AddField( fieldmlInfo, baseName//".geometric", region, mesh, GeometricField, &
      & CMISSFieldUVariableType, err )

    CALL CMISSFieldmlOutput_AddField( fieldmlInfo, baseName//".dependent", region, mesh, DependentField, &
      & CMISSFieldUVariableType, err )

    CALL CMISSFieldmlOutput_AddField( fieldmlInfo, baseName//".independent", region, mesh, IndependentField, &
      & CMISSFieldUVariableType, err )

    CALL CMISSFieldmlOutput_AddField( fieldmlInfo, baseName//".source", region, mesh, SourceField, &
      & CMISSFieldUVariableType, err )

    CALL CMISSFieldmlOutput_AddField( fieldmlInfo, baseName//".materials", region, mesh, MaterialsField, &
      & CMISSFieldUVariableType, err )

    !CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".analytic", region, mesh, AnalyticField, &
    !  & CMISSFieldUVariableType, err )
    
    CALL CMISSFieldmlOutput_Write( fieldmlInfo, outputFilename, err )
    
    CALL CMISSFieldmlUtil_FinaliseInfo( fieldmlInfo, err )

  ENDIF

  !CALL CMISSFinalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM STATICADVECTIONDIFFUSIONEXAMPLE
