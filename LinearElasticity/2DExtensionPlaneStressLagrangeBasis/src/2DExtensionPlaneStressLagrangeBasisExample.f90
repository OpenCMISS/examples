! \file
!> \author Chris Bradley
!> \brief This is an example program to solve a linear elasticity equation using openCMISS calls.
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

!> \example LinearElasticity/src/LinearElasticityExample.f90
!! Example program to solve a linear elasticity equation using openCMISS calls.
!<

!> Main program
PROGRAM LinearElasticity2DExtensionPlaneStressLagrangeBasis

  USE MPI
  USE OPENCMISS

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: ORIGIN(3)=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
  REAL(CMISSDP), PARAMETER :: LENGTH=120.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=160.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=2
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Basis1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: Basis2UserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=2
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfNodes=4
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfElements=1
  INTEGER(CMISSIntg), PARAMETER :: MeshComponent1UserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshComponent2UserNumber=2

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=2

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariables=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponents=2

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=4

  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  REAL(CMISSDP), PARAMETER ::   ZERO = 0.0_CMISSDP

  !Program types


  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: node_idx,component_idx
  REAL(CMISSDP) :: Value1,Value2
  LOGICAL :: EXPORT_FIELD

  !CMISS variables

  TYPE(CMISSRegionType) :: WorldRegion
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  TYPE(CMISSBasisType) :: Basis(2)
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region
  TYPE(CMISSFieldType) :: EquationsSetField
  TYPE(CMISSSolverType) :: Solver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSMeshElementsType) :: Elements(2)

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

  !Intialise cmiss
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  CALL CMISSDiagnosticsSetOn(CMISSFromDiagType,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_FINITE_ELEMENT_CALCULATE"/),Err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberGlobalXElements=1
  NumberGlobalYElements=1
  NumberOfDomains=1

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemTypeSet(CoordinateSystem,CMISSCoordinateRectangularCartesianType,Err)
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystemOriginSet(CoordinateSystem,(ORIGIN),Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Create a region and assign the CS to the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !Define 2 sets of basis functions, one describing each independent coordinate specified by InterpolationType
  !NOTE if you change interpolation you need to change Boundary Conditions
  !NOTE:: Num of Gauss points must be the same across X,Y & Z coordinates and be sufficient to accurately integrate the hightest order interpolation being used

  CALL CMISSBasisTypeInitialise(Basis(1),Err)
  CALL CMISSBasisCreateStart(Basis1UserNumber,Basis(1),Err)
  CALL CMISSBasisTypeSet(Basis(1),CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(Basis(1),NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(Basis(1),(/CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis(1),(/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
  CALL CMISSBasisCreateFinish(Basis(1),Err)

  CALL CMISSBasisTypeInitialise(Basis(2),Err)
  CALL CMISSBasisCreateStart(Basis2UserNumber,Basis(2),Err)
  CALL CMISSBasisTypeSet(Basis(2),CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(Basis(2),NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(Basis(2),(/CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis(2),(/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
  CALL CMISSBasisCreateFinish(Basis(2),Err)

  !Create a mesh
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NumberOfMeshDimensions,Mesh,Err)

  CALL CMISSMeshNumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  CALL CMISSMeshNumberOfElementsSet(Mesh,TotalNumberOfElements,Err)

  !Define nodes for the mesh
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSNodesCreateStart(Region,TotalNumberOfNodes,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)

  !Create elements for the mesh
  !Mesh Component 1
  CALL CMISSMeshElementsTypeInitialise(Elements(1),Err)
  CALL CMISSMeshElementsCreateStart(Mesh,MeshComponent1UserNumber,Basis(1),Elements(1),Err)
  CALL CMISSMeshElementsNodesSet(Elements(1),1,(/1,2,3,4/),Err)
  CALL CMISSMeshElementsCreateFinish(Elements(1),Err)
  !Mesh Component 2
  CALL CMISSMeshElementsTypeInitialise(Elements(2),Err)
  CALL CMISSMeshElementsCreateStart(Mesh,MeshComponent2UserNumber,Basis(2),Elements(2),Err)
  CALL CMISSMeshElementsNodesSet(Elements(2),1,(/1,2,3,4/),Err)
  CALL CMISSMeshElementsCreateFinish(Elements(2),Err)

  CALL CMISSMeshCreateFinish(Mesh,Err) 

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)  
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldVariableLabelSet(GeometricField,CMISSFieldUVariableType,"Geometry",Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,2,Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Set geometric node coordinates (x)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,1,LENGTH,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,1,LENGTH,Err)

  !Set geometric node coordinates (y)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,3,2,WIDTH,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,4,2,WIDTH,Err)

  !Create a dependent field with two variables and three components
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSFieldCreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL CMISSFieldTypeSet(DependentField,CMISSFieldGeneralType,Err)  
  CALL CMISSFieldMeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(DependentField,GeometricField,Err) 
  CALL CMISSFieldDependentTypeSet(DependentField,CMISSFieldDependentType,Err) 
  CALL CMISSFieldNumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
  CALL CMISSFieldVariableLabelSet(DependentField,CMISSFieldUVariableType,"Dependent",Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldUVariableType,FieldDependentNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldDelUDelNVariableType,FieldDependentNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,2,2,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,2,2,Err)
  CALL CMISSFieldCreateFinish(DependentField,Err)

  !Create a material field and attach it to the geometric field  
  CALL CMISSFieldTypeInitialise(MaterialField,Err)
  CALL CMISSFieldCreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL CMISSFieldTypeSet(MaterialField,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialField,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(MaterialField,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL CMISSFieldVariableLabelSet(MaterialField,CMISSFieldUVariableType,"Material",Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialField,CMISSFieldUVariableType,FieldMaterialNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(MaterialField,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialField,CMISSFieldUVariableType,2,2,Err)
  CALL CMISSFieldCreateFinish(MaterialField,Err)

  !Set isotropic elasticity material parameters - Young's Modulus & Poisson's Ratio, thickness
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,0.036_CMISSDP,Err) !thickness
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,30E6_CMISSDP,Err) !E
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,0.25_CMISSDP,Err) !v


  !Create a Elasticity Class, Linear Elasticity type, no subtype, EquationsSet
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationSetUserNumber,Region,GeometricField,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetLinearElasticityType,CMISSEquationsSetPlaneStressSubtype,EquationsSetFieldUserNumber,EquationsSetField, &
    & EquationsSet,Err)
  
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquationsSparsityTypeSet(EQUATIONS,CMISSEquationsSparseMatrices,Err)
                                              !CMISSEquationsSparseMatrices=1 !<Use sparse matrices for the equations.
                                              !CMISSEquationsFullMatrices=2 !<Use fully populated matrices for the equations. 
  CALL CMISSEquationsOutputTypeSet(EQUATIONS,CMISSEquationsElementMatrixOutput,Err)
                                            !CMISSEquationsNoOutput !<No output from the equations.
                                            !CMISSEquationsTimingOutput !<Timing information output.
                                            !CMISSEquationsMatrixOutput !<All below and equation matrices output.
                                            !CMISSEquationsElementMatrixOutput !<All below and Element matrices output.
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)
  
  !Define the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemElasticityClass,CMISSProblemLinearElasticityType, &
    & CMISSProblemNoSubtype,Err)
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !Start the creation of the Problem solvers
  !Create the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverOutputTypeSet(SOLVER,CMISSSolverSolverMatrixOutput,Err)
                                      !CMISSSolverNoOutput !<No output from the solver routines. \see OPENCMISS_SolverOutputTypes,OPENCMISS
                                      !CMISSSolverProgressOutput !<Progress output from solver routines.
                                      !CMISSSolverTimingOutput !<Timing output from the solver routines plus below.
                                      !CMISSSolverSolverOutput !<Solver specific output from the solver routines plus below.
                                      !CMISSSolverSolverMatrixOutput !<Solver matrices output from the solver routines plus below.
  CALL CMISSSolverLibraryTypeSet(SOLVER,CMISSSolverPETScLibrary,Err)
  CALL CMISSSolverLinearTypeSet(SOLVER,CMISSSolverLinearDirectSolveType,Err)
                                      !CMISSSolverLinearDirectSolveType    !<Direct linear solver type.
                                      !CMISSSolverLinearIterativeSolveType !<Iterative linear solver type.
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)   
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
                                                          !CMISSSolverEquationsSparseMatrices !<Use sparse solver matrices.
                                                          !CMISSSolverEquationsFullMatrices !<Use fully populated solver matrices.
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Fix nodes 1,3 at x=0
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,1,1, &
    & CMISSBoundaryConditionFixed,ZERO,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,3,1, &
    & CMISSBoundaryConditionFixed,ZERO,Err)

  !Fix nodes 1,2 at y=0
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,1,2, &
    & CMISSBoundaryConditionFixed,ZERO,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,2,2, &
    & CMISSBoundaryConditionFixed,ZERO,Err)

  !Apply force at nodes 2,4 at x=l
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldDelUDelNVariableType,1,1,2,1, &
    & CMISSBoundaryConditionFixed, &
    & -800.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldDelUDelNVariableType,1,1,4,1, &
    & CMISSBoundaryConditionFixed, &
    & -800.0_CMISSDP,Err)

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !=SOLVE Problem==================================================================================================================
  !Solve the Problem
  CALL CMISSProblemSolve(Problem,Err)

  !Add dependent field to geometric field
  DO node_idx=1,TotalNumberOfNodes
    DO component_idx=1,FieldGeometryNumberOfComponents
      CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx, &
        & component_idx,Value1,Err)
      CALL CMISSFieldParameterSetGetNode(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx, &
        & component_idx,Value2,Err)
      CALL CMISSFieldParameterSetAddNode(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx, &
        & component_idx,Value1+Value2,Err)
    ENDDO
  ENDDO

  !=OUTPUT SOLUTION================================================================================================================
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"LinearElasticity2DExtensionPlaneStressLagrangeBasisExample","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"LinearElasticity2DExtensionPlaneStressLagrangeBasisExample","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
  ENDIF
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM LinearElasticity2DExtensionPlaneStressLagrangeBasis


