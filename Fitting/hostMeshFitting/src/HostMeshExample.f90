!> \file
!> \author Hashem Yousefi
!> \brief This is an example program to solve/simulate a test problem for
!> Host-Mesh fitting technique in 2D using OpenCMISS calls.
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

!> \OpenCMISS/examples/Fitting/GeometricFitting/Host-MeshExample.f90
!! Example program to update enslaved mesh inside Host Element using  
!! OpenCMISS calls.
!!
!<

!> Main program
PROGRAM HOSTMESH_FITTING_EXAMPLE

  USE OPENCMISS
  USE MPI

  IMPLICIT NONE
										
 !Parameters
 ! =====================================================================
 ! statring the control panel ... 
 ! =====================================================================

  REAL(CMISSDP), PARAMETER :: HEIGHT = 1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH = 1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH = 1.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber = 4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: DataPointFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: MaterialFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: DependentDataFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: DataProjectionUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=17
   
  !Program types
  !Program variables
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_Y_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: INTERPOLATION_TYPE

  !CMISS variables
  TYPE(CMISSBasisType) :: Basis
!  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType)    :: Decomposition
  TYPE(CMISSDataPointsType)       :: DataPoints
  TYPE(CMISSDataProjectionType)   :: DataProjection
  TYPE(CMISSEquationsType)        :: Equations
  TYPE(CMISSEquationsSetType)     :: EquationsSet
  TYPE(CMISSFieldType)            :: GeometricField,EquationsSetField,DependentField
  TYPE(CMISSFieldType)            :: IndependentField, MaterialField
  TYPE(CMISSFieldsType)           :: Fields
  TYPE(CMISSGeneratedMeshType)    :: GeneratedMesh  
  TYPE(CMISSMeshType)             :: Mesh
  TYPE(CMISSNodesType)            :: Nodes
  TYPE(CMISSProblemType)          :: Problem
  TYPE(CMISSRegionType)           :: Region, WorldRegion
  TYPE(CMISSSolverType)           :: Solver
  TYPE(CMISSSolverEquationsType)  :: SolverEquations
  
  !Generic CMISS variables
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NumberofDataPoints, NumberOfProjectedDataPoints 
  INTEGER(CMISSIntg) :: DataPointID, DataPointsNumber, DataPointsNumberIndex
  INTEGER(CMISSIntg) :: NumberOfDimensions, NumberOfGaussXi, NumberOfIterations
  INTEGER(CMISSIntg) :: Element_Domain, ElementID
  INTEGER(CMISSIntg) :: EquationsSetIndex, SolverIndex
  INTEGER(CMISSIntg) :: i, j , k  
  REAL(CMISSDP)      :: ValueSet, ZeroTolerance, Tau, Kappa,x,y,z
  REAL(CMISSDP)      :: absoluteTolerance, relativeTolerance
  
  INTEGER(CMISSIntg) :: Err
  REAL (CMISSDP), Allocatable :: A(:,:),B(:,:)
  REAL (CMISSDP), Allocatable :: DataPointLocations(:,:)
  REAL (CMISSDP), Allocatable :: location(:)
! based on the geometry number of Global X, Y and Z elements can be defined 
! and also the interpolation type can be defined ...
  NumberofDataPoints = 4
  NumberOfComputationalNodes = 4
  NumberOfDimensions = 2
  NumberOfGaussXi = 2
  NUMBER_OF_ELEMENTS = 1
  NUMBER_GLOBAL_X_ELEMENTS = 1
  NUMBER_GLOBAL_Y_ELEMENTS = 1
  NUMBER_GLOBAL_Z_ELEMENTS = 0
  INTERPOLATION_TYPE = 1
	
! default of the number of iterations would be one .. 	
  NumberOfIterations = 1
! at this step we would like to neglect smoothing terms in objective function...
  Tau = 0.0
  Kappa = 0.0
  ZeroTolerance = 0.00001
  Absolutetolerance = 1.0E-10 
  relativetolerance = 1.0E-5

  allocate (A(NumberofDataPoints,NumberofDimensions))
  allocate (B(NumberofDataPoints,NumberofDimensions))
  allocate (DataPointLocations(NumberOfDataPoints,NumberOfDimensions))
  allocate (location(NumberOfDimensions))

! Getting the datapoints from a file for the location of datapoints before and after deformation  
  
!Gib: do not use file numbers < 10
  OPEN (10,FILE='/people/hyou267/OpenCMISS/Examples/Fitting/HostMesh/datafile1.txt',STATUS = 'OLD')
  DO i = 1, NumberofDataPoints
    READ (10,*) A(i,1), A(i,2)
  END DO
  CLOSE (10)
	
  OPEN (11,FILE='/people/hyou267/OpenCMISS/Examples/Fitting/HostMesh/datafile2.txt',STATUS = 'OLD')
  DO i = 1, NumberofDataPoints
    READ (11,*) B(i,1), B(i,2)
  END DO
  CLOSE (11)

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  CALL CMISSRandomSeedsSet(9999,Err)

 
	!Gib: You do not write to a file name, you write to a file number.  You must open the file for writing first,
	! or you can do write(*,... which writes to the console.  In any case filename is neither declared nor defined.
  WRITE(*,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "HostMeshFitting","NUMBER_GLOBAL_X_ELEMENTS=", NUMBER_GLOBAL_X_ELEMENTS, & 
    & "NUMBER_GLOBAL_Y_ELEMENTS=",  NUMBER_GLOBAL_Y_ELEMENTS,"NUMBER_GLOBAL_Z_ELEMENTS=", NUMBER_GLOBAL_Z_ELEMENTS, &
    & "INTERPOLATION_TYPE=",INTERPOLATION_TYPE, "NumberofIterations=", NumberOfIterations
  write(*,*) 'A: '
  WRITE (*,*) A
  !Gib: I suggest you put STOP here to debug the above code, and make sure that you are writing out the data you read in correctly
  ! before attempting to go further.
  STOP
  
  
  CALL CMISSOutputSetOn("test.out",Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)
! ======================================================================    
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
  CALL CMISSRegion_LabelSet(Region,"HostRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

! ======================================================================
! Mesh 
! ======================================================================

  NUMBER_OF_ELEMENTS = 1
  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1,2,3,4)
    CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bi-interpolation basis
    CALL CMISSBasis_NumberOfXiSet(Basis,2,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
!    IF(NUMBER_OF_GAUSS_XI>0) THEN
!      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
!    ENDIF
  ELSE
    !Set the basis to be a tri-interpolation basis
    CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
!    IF(NUMBER_OF_GAUSS_XI>0) THEN
!      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
!    ENDIF
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
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

! ====================================================================== 
! Geometry Field 
! ======================================================================

  !Destory the mesh now that we have decomposed it
  !CALL CMISSMesh_Destroy(Mesh,Err)
 
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

! =================================

  CALL CMISSDataPoints_Initialise(DataPoints,Err)
  CALL CMISSDataPoints_CreateStartNumber(regionUserNumber,numberOfDataPoints,Err)
! calculation for the data points for the corners of the square as the Host element 
  DO i = 1, NumberOfDataPoints
     IF (NumberOfDimensions == 3) THEN 
         x = A(i,1)
         y = A(i,2)
         z = A(i,3)
         DataPointLocations(i,:) = [x,y,z]
     END IF
     IF (NumberOfDimensions == 2) THEN
        x = A(i,1)
        y = A(i,2)
        DataPointLocations(i,:) = [x,y]
     END IF
  END DO
  DO i = 1, NumberOfDataPoints
    location(:) = DataPointLocations(i,:)
    CALL CMISSDataPoints_ValuesSet(DataPoints,i,location,Err)
  END DO
!!! now the data points are inside the A matrix ... and it should be called inside datapoints PTR

  CALL CMISSDataPoints_CreateFinishNumber(regionUserNumber,numberOfDataPoints,Err)


  ! ====================================================================


  CALL CMISSDataProjection_Initialise(DataProjection,Err)
  CALL CMISSDataProjection_CreateStart(dataProjectionUserNumber,dataPoints,mesh,dataprojection,Err)
  CALL CMISSDataProjection_ProjectionTypeSetObj(dataProjection,CMISS_DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE,Err)
  CALL CMISSDataProjection_CreateFinishNumber(dataProjectionUserNumber,regionUserNumber,Err)


  CALL CMISSDataProjection_ProjectionEvaluate(dataProjection,geometricfield,Err)
  CALL CMISSMesh_TopologyDataPointsCalculateProjection(Mesh,DataProjection,Err)
  CALL CMISSDecomposition_TopologyDataProjectionCalculate(Decomposition,Err)
  Write(*,*) 'projection complete'

! ==================================
  
  !Create the Finite Element Calculation for Fitting and Equations set
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,CMISS_EQUATIONS_SET_FITTING_CLASS, & 
    CMISS_EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE,CMISS_EQUATIONS_SET_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE, &
    EquationsSetFieldUserNumber, EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

! ======================================================================
  !Create the equations set dependent field variables
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Set the DOFs to be contiguous across components
  CALL CMISSField_VariableLabelSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,'Dependent',Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)
  CALL CMISSField_ComponentValuesInitialiseDP(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.0)
  ! in order to initialise 
  DO i=1, NumberOfDimensions
    CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE, &
      CMISS_FIELD_VALUES_SET_TYPE,i, DependentField, CMISS_FIELD_U_VARIABLE_TYPE, CMISS_FIELD_VALUES_SET_TYPE,i,Err) 
  END DO
  
! ======================================================================
! Independent Field 
! Create DataPoint Field, 

  CALL CMISSField_Initialise(IndependentField,Err)
  CALL CMISSEquationsSet_IndependentCreateStart(EquationsSet,IndependentFieldUserNumber,IndependentField,Err)
  CALL CMISSField_VariableLabelSet(IndependentField,CMISS_FIELD_U_VARIABLE_TYPE,'data point vector',Err)
  CALL CMISSField_VariableLabelSet(IndependentField,CMISS_FIELD_V_VARIABLE_TYPE,'data point weight',Err)
  CALL CMISSField_DataProjectionSet(IndependentField,DataProjection,Err)
  CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSet,Err)
  ! initialise data point vector 
  CALL CMISSField_ComponentValuesInitialiseDP(IndependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.0)
  ! initialise data point weight
  CALL CMISSField_ComponentValuesInitialiseDP(IndependentField,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1.0)  
  !Initialise the field with an initial guess

! update independent field ... 
  Do i=1, NUMBER_OF_ELEMENTS
    ElementID = i
    CALL CMISSDecomposition_ElementDomainGet(Decomposition,ElementID,Element_Domain,Err)
    IF (Element_Domain == ComputationalNodeNumber) THEN 
        CALL CMISSDecomposition_TopologyNumberOfElementDataPointsGet(Decomposition,ElementID,NumberOfProjectedDataPoints,err)
        DO j = 1, NumberOfProjectedDataPoints
            DataPointID = j
            CALL CMISSDecomposition_TopologyElementDataPointUserNumberGet(Decomposition,ElementID,DataPointID,DataPointsNumber,err)
            CALL CMISSDataPoints_ValuesGet(DataPoints,DataPointsNumber,location,Err)
            Do k = 1, NumberOfDimensions
 !               DataPointNumberIndex = DataPointsNumber - 1
                
                CALL CMISSField_ParameterSetUpdateElementDataPoint(IndependentField,CMISS_FIELD_U_VARIABLE_TYPE, & 
                  CMISS_FIELD_VALUES_SET_TYPE,i,j,k,location(k),Err)
            END DO
        END DO
     END IF
  END DO
! ======================================================================
! Material Field
  CALL CMISSField_Initialise(MaterialField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,MaterialFieldUserNumber,MaterialField,Err)
  CALL CMISSField_VariableLabelSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,'Smoothing Parameters',Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(equationsSet,Err)
 
! set Kappa and Tau for Sobolev smoothing ... 
  CALL CMISSField_ComponentValuesInitialiseDP(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,tau)
  CALL CMISSField_ComponentValuesInitialiseDP(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Kappa)
  
! Equations 
  ! Create Equations 
  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquations_SparsityTypeSet(equations,CMISS_EQUATIONS_FULL_MATRICES,Err)   
  CALL CMISSEquations_OutputTypeSet(equations,CMISS_EQUATIONS_NO_OUTPUT,Err)         
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)
 
! problem Setup 
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)  
  !CALL CMISSControlLoop_Initialise(Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_FITTING_CLASS,CMISS_PROBLEM_DATA_FITTING_TYPE, & 
  & CMISS_PROBLEM_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(Problem,Err)
   
! creating control loops 
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)
  
 ! creating problem solver 
 !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT,Err)
  CALL CMISSSolver_LinearTypeSet(Solver,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)
  CALL CMISSSolver_LibraryTypeSet(Solver,CMISS_SOLVER_UMFPACK_LIBRARY,Err)
  CALL CMISSSolver_LinearIterativeAbsoluteToleranceSetObj(Solver,absoluteTolerance,err)
  CALL CMISSSolver_LinearIterativeRelativeToleranceSetObj(Solver,relativetolerance,err) 
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)
  
! now solver equations   
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(ProblemUserNumber,err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_EQUATIONS_FULL_MATRICES,Err) 
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,err)
  CALL CMISSProblem_SolverEquationsCreateFinish(ProblemUserNumber,err)
  
  
  ! ====================================================================
!   BoundaryConditions can be defined later ... 
! ======================================================================
! iteration for doing the optimisation process
! Solve the problem
  DO j=1, NumberOfIterations
    write(*,'(a,i3)') 'Solving HostMesh Fitting Problem, Iteration: ',j
    CALL CMISSProblem_Solve(Problem,Err)
    DO i = 1, NumberOfDimensions
        CALL CMISSField_ParametersToFieldParametersComponentCopy(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,  & 
           CMISS_FIELD_VALUES_SET_TYPE,i,DependentField, CMISS_FIELD_U_VARIABLE_TYPE, CMISS_FIELD_VALUES_SET_TYPE,i,Err) 
    END DO
!Export results
    CALL CMISSFields_Initialise(Fields,Err)  
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"Deformed Geometry","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"Deformed Geometry","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  END DO
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
    
END PROGRAM HOSTMESH_FITTING_EXAMPLE
