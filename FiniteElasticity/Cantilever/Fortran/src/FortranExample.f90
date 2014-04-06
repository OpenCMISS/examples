!> \file
!> \author Adam Reeve
!> \brief This is an example program to solve a finite elasticity equation using openCMISS calls.
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
!> Contributor(s): Adam Reeve
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

!> \example FiniteElasticity/Cantilever/src/CantileverExample.f90
!! Example program to solve a finite elasticity equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/Cantilever/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/Cantilever/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM CANTILEVEREXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP) :: Width
  REAL(CMISSDP) :: Length
  REAL(CMISSDP) :: Height
  REAL(CMISSDP) :: YoungsModulus,PoissonsRatio,ShearModulus,MooneyRivlin1,MooneyRivlin2
  INTEGER(CMISSIntg) :: DisplacementInterpolationType
  INTEGER(CMISSIntg) :: PressureInterpolationType
  INTEGER(CMISSIntg) :: PressureMeshComponent
  INTEGER(CMISSIntg) :: NumberOfGaussXi,TopNormalXi
  INTEGER(CMISSIntg) :: ScalingType,ForceNode
  REAL(CMISSDP) :: Density
  REAL(CMISSDP), PARAMETER :: Gravity(3)=[0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP]
  REAL(CMISSDP) :: ExternalForce,Load,EI,v,F_d,C_d,u_y
  INTEGER(CMISSIntg), PARAMETER :: NumberOfLoadIncrements=5
  INTEGER(CMISSIntg), ALLOCATABLE :: TopNodes(:)

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DisplacementBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldSourceUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx,component_idx,deriv_idx
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg) :: LeftNormalXi
  INTEGER(CMISSIntg) :: NumberOfArguments,ArgumentLength,ArgStatus
  CHARACTER(LEN=255) :: CommandArgument
  LOGICAL :: SetupFlag=.FALSE.

  !CMISS variables
  TYPE(CMISSBasisType) :: DisplacementBasis,PressureBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,FibreField,MaterialField,DependentField,SourceField,EquationsSetField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSControlLoopType) :: ControlLoop

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: Err

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
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
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
  CALL CMISSOutputSetOn("Cantilever",Err)

  !Read in arguments and overwrite default values
  !Usage: CantileverExample [Displacement Interpolation Type] [X elements] [Y elements] [Z elements] [Scaling Type]
  !Defaults:
  DisplacementInterpolationType=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  NumberGlobalXElements=3
  NumberGlobalYElements=2
  NumberGlobalZElements=1
  ScalingType=CMISS_FIELD_ARITHMETIC_MEAN_SCALING
  
  !=================================================================================================================================
  
  
  !Set solver parameters
!    RelativeTolerance=1.0E-4_CMISSDP !default: 1.0E-05_CMISSDP
!    AbsoluteTolerance=1.0E-4_CMISSDP !default: 1.0E-10_CMISSDP
!    DivergenceTolerance=1.0E5 !default: 1.0E5
!    MaximumIterations=100000000 !default: 100000
!    MaxFunctionEvaluations=100000
!    RestartValue=30 !default: 30
!    LinesearchAlpha=1.0_CMISSDP
    
    
  !=================================================================================================================================
  
  SetupFlag=.TRUE.
  IF(SetupFlag) THEN
    DisplacementInterpolationType=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
    NumberGlobalXElements=2
    NumberGlobalYElements=1
    NumberGlobalZElements=1
    ScalingType=CMISS_FIELD_ARITHMETIC_MEAN_SCALING
    !===========================================================================================================================================
    Width=2000.0_CMISSDP ! m
    Length=200.0_CMISSDP ! m
    Height=200.0_CMISSDP ! m
    !Solid
    Density=300.0E-6_CMISSDP ! kg / m^3
    !Young's modulus E
    YoungsModulus=2.0E1_CMISSDP ! N / m^2
    !Poisson's ratio
    PoissonsRatio=0.3 ! [.]
    !Homogenous, isotropic material G=E/(2*(1+poissonsRatio))
    ShearModulus=YoungsModulus/(2.0_CMISSDP*(1.0_CMISSDP+PoissonsRatio)) ! N / m^2
    !Neo-Hookean material: c1=G/2 c2=0
    MooneyRivlin1=0.5_CMISSDP*ShearModulus ! N / m^2
    MooneyRivlin2=0.0_CMISSDP !
    
    EI=(YoungsModulus*(Height)**3.0_CMISSDP)/(12.0_CMISSDP*(1-PoissonsRatio)) ! Nm
    C_d=1.05_CMISSDP
    C_d=9.46_CMISSDP
    C_d=3.26_CMISSDP
    v=0.1_CMISSDP ! m / s
    F_d=0.5_CMISSDP*C_d*(Width*Length/1000/1000)*100.0_CMISSDP*v*v ! N
    u_y=-(F_d*(Width)**4)/(8.0_CMISSDP*EI*(Width*Length)) ! m
    PRINT *,"Constant load, u_y = ",u_y," m"
    
    !+========
    !point force
    EI=YoungsModulus/12.0_CMISSDP*(Height)*(Length)**3.0_CMISSDP
    u_y=-(F_d*(Width)**3.0_CMISSDP)/(3.0_CMISSDP*EI)
    
    ExternalForce=-F_d/(DisplacementInterpolationType+1)*1000!N
    Load=ExternalForce/(Width*Length)
    !===========================================================================================================================================
    Width=2.0_CMISSDP ! m
    Length=0.2_CMISSDP ! m
    Height=0.2_CMISSDP ! m
    !Solid
    Density=300.0_CMISSDP ! kg / m^3
    !Young's modulus E
    YoungsModulus=2.0E4_CMISSDP ! N / m^2
    !Poisson's ratio
    PoissonsRatio=0.3 ! [.]
    !Homogenous, isotropic material G=E/(2*(1+poissonsRatio))
    ShearModulus=YoungsModulus/(2.0_CMISSDP*(1.0_CMISSDP+PoissonsRatio)) ! N / m^2
    !Neo-Hookean material: c1=G/2 c2=0
    MooneyRivlin1=0.5_CMISSDP*ShearModulus ! N / m^2
    MooneyRivlin2=0.0_CMISSDP !
    
    EI=(YoungsModulus*(Height)**3.0_CMISSDP)/(12.0_CMISSDP*(1-PoissonsRatio**2)) ! Nm
  !  C_d=1.05_CMISSDP
  !  C_d=9.46_CMISSDP
    C_d=3.26_CMISSDP
    v=0.1_CMISSDP ! m / s
    F_d=0.5_CMISSDP*C_d*(Width*Length)*100.0_CMISSDP*v*v ! N
    u_y=-(F_d*(Width)**4)/(8.0_CMISSDP*EI*(Width*Length)) ! m
    PRINT *,"Constant load, u_y = ",u_y," m"
    
    !+========
    !point force
    EI=YoungsModulus/12.0_CMISSDP*(Height)*(Length)**3.0_CMISSDP
    u_y=-(F_d*(Width)**3.0_CMISSDP)/(3.0_CMISSDP*EI)
    
    ExternalForce=-F_d/(DisplacementInterpolationType+1)!N
    Load=ExternalForce/(Width*Length)
  ENDIF
  
  
    

  
    PRINT *, "Width = ",Width," m"
    PRINT *, "Length = ",Length," m"
    PRINT *, "Height = ",Height," m"
    PRINT *, "Rho: ",Density," kg m^-3"
    PRINT *, "E = ",YoungsModulus," N/m^2"
    PRINT *, "nu = ",PoissonsRatio
    PRINT *, "c10 = ",MooneyRivlin1," N/m^2"
    PRINT *, "c01 = ",MooneyRivlin2
    PRINT *, "EI = ",EI," N"
    PRINT *, "C_d  =",C_d
    PRINT *, "v = ",v," m/s"
    PRINT *, "F_d = ",F_d," N"
    PRINT *, "ExternalForce = ",ExternalForce,' N (point force)'
    PRINT *, "Point force, u_y = ",u_y," m"
  
  
  
  
  !=================================================================================================================================

  NumberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(NumberOfArguments >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(CommandArgument(1:ArgumentLength),*) DisplacementInterpolationType
  ENDIF
  IF(NumberOfArguments >= 2) THEN
    CALL GET_COMMAND_ARGUMENT(2,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalXElements
    IF(NumberGlobalXElements<1) CALL HANDLE_ERROR("Invalid number of X elements.")
  ENDIF
  IF(NumberOfArguments >= 3) THEN
    CALL GET_COMMAND_ARGUMENT(3,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalYElements
    IF(NumberGlobalYElements<1) CALL HANDLE_ERROR("Invalid number of Y elements.")
  ENDIF
  IF(NumberOfArguments >= 4) THEN
    CALL GET_COMMAND_ARGUMENT(4,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 4.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalZElements
    IF(NumberGlobalZElements<1) CALL HANDLE_ERROR("Invalid number of Z elements.")
  ENDIF
  IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
    IF(NumberOfArguments >= 5) THEN
      CALL GET_COMMAND_ARGUMENT(5,CommandArgument,ArgumentLength,ArgStatus)
      IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 5.")
      READ(CommandArgument(1:ArgumentLength),*) ScalingType
      IF(ScalingType<0.OR.ScalingType>5) CALL HANDLE_ERROR("Invalid scaling type.")
    ENDIF
  ELSE
    ScalingType=CMISS_FIELD_NO_SCALING
  ENDIF
  SELECT CASE(DisplacementInterpolationType)
  CASE(CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION)
    NumberOfGaussXi=2
    PressureMeshComponent=1
  CASE(CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
    NumberOfGaussXi=3
    PressureMeshComponent=2
    PressureInterpolationType=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  CASE(CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)
    NumberOfGaussXi=4
    PressureMeshComponent=2
    !Should generally use quadratic interpolation but use linear to match CMISS example 5e
    PressureInterpolationType=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  CASE DEFAULT
    NumberOfGaussXi=0
    PressureMeshComponent=1
  END SELECT
  WRITE(*,'("Interpolation: ", i3)') DisplacementInterpolationType
  WRITE(*,'("Elements: ", 3 i3)') NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  WRITE(*,'("Scaling type: ", i3)') ScalingType

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains=NumberOfComputationalNodes

  !Create a 3D rectangular cartesian coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_LabelSet(Region,"Region",Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Define basis function for displacement
  CALL CMISSBasis_Initialise(DisplacementBasis,Err)
  CALL CMISSBasis_CreateStart(DisplacementBasisUserNumber,DisplacementBasis,Err)
  SELECT CASE(DisplacementInterpolationType)
  CASE(1,2,3,4)
    CALL CMISSBasis_TypeSet(DisplacementBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL CMISSBasis_TypeSet(DisplacementBasis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  CALL CMISSBasis_NumberOfXiSet(DisplacementBasis,3,Err)
  CALL CMISSBasis_InterpolationXiSet(DisplacementBasis,[DisplacementInterpolationType,DisplacementInterpolationType, &
      & DisplacementInterpolationType],Err)
  IF(NumberOfGaussXi>0) THEN
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(DisplacementBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL CMISSBasis_CreateFinish(DisplacementBasis,Err)

  IF(PressureMeshComponent/=1) THEN
    !Basis for pressure
    CALL CMISSBasis_Initialise(PressureBasis,Err)
    CALL CMISSBasis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
    SELECT CASE(PressureInterpolationType)
    CASE(1,2,3,4)
      CALL CMISSBasis_TypeSet(PressureBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CASE(7,8,9)
      CALL CMISSBasis_TypeSet(PressureBasis,CMISS_BASIS_SIMPLEX_TYPE,Err)
    END SELECT
    CALL CMISSBasis_NumberOfXiSet(PressureBasis,3,Err)
    CALL CMISSBasis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType, &
        & PressureInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
    CALL CMISSBasis_CreateFinish(PressureBasis,Err)
  ENDIF

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  IF(PressureMeshComponent==1) THEN
    CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,[DisplacementBasis],Err)
  ELSE
    CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,[DisplacementBasis,PressureBasis],Err)
  ENDIF
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[Width,Height],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[Width,Length,Height],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
      & NumberGlobalZElements],Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSField_VariableLabelSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL CMISSField_ScalingTypeSet(GeometricField,ScalingType,Err)
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Create a fibre field and attach it to the geometric field
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_VariableLabelSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL CMISSField_ScalingTypeSet(FibreField,ScalingType,Err)
  CALL CMISSField_CreateFinish(FibreField,Err)

  !Create the equations_set
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EquationsSetFieldUserNumber, &
      & EquationsSetField, &
    & EquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSField_VariableLabelSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  DO component_idx=1,3
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,component_idx,1,Err)
  ENDDO
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,PressureMeshComponent,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,PressureMeshComponent,Err)
  IF(PressureMeshComponent==1) THEN
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION, &
      & Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ELSE
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
  ENDIF
  CALL CMISSField_ScalingTypeSet(DependentField,ScalingType,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL CMISSField_Initialise(MaterialField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL CMISSField_VariableLabelSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL CMISSField_VariableLabelSet(MaterialField,CMISS_FIELD_V_VARIABLE_TYPE,"Density",Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & MooneyRivlin1,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
    & MooneyRivlin2,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Density,Err)

  !Create the source field with the gravity vector
  CALL CMISSField_Initialise(SourceField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(EquationsSet,FieldSourceUserNumber,SourceField,Err)
  CALL CMISSField_ScalingTypeSet(SourceField,ScalingType,Err)
  CALL CMISSEquationsSet_SourceCreateFinish(EquationsSet,Err)
  DO component_idx=1,3
    CALL CMISSField_ComponentValuesInitialise(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
        & component_idx,Gravity(component_idx),Err)
  ENDDO

  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4, &
    & -MooneyRivlin1, &
    & Err)
  CALL CMISSField_ParameterSetUpdateStart(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Define the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_ELASTICITY_CLASS,CMISS_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMISS_PROBLEM_NO_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL CMISSControlLoop_MaximumIterationsSet(ControlLoop,NumberOfLoadIncrements,Err)
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
!  CALL CMISSSolver_NewtonLineSearchTypeSet(Solver,CMISS_SOLVER_NEWTON_LINESEARCH_LINEAR,Err)
  CALL CMISSSolver_NewtonMaximumFunctionEvaluationsSet(Solver,100000,Err)
!  CALL CMISSSolver_OutputTypeSet(Solver,NonlinearSolver_OutputType,Err)
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(Solver,1.0E-4_CMISSDP,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(Solver,1.0E-5_CMISSDP,Err)
  CALL CMISSSolver_NewtonMaximumIterationsSet(Solver,100000,Err)
!  CALL CMISSSolver_NewtonLineSearchAlphaSet(Solver,LinesearchAlpha,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)

  !Fix x=0 nodes in x, y and z
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      DO component_idx=1,3
        CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
          & component_idx,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        IF(DisplacementInterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
          DO deriv_idx=3,8
            CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,deriv_idx, &
              & NodeNumber, &
              & component_idx,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_TOP_SURFACE,TopNodes,TopNormalXi,Err)
  !DO node_idx=1,SIZE(TopNodes,1)
  IF(.TRUE.) THEN
    IF(SetupFlag) THEN
      ForceNode=9
      IF(DisplacementInterpolationType==CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION) ForceNode=35
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1, &
        & ForceNode,3,CMISS_BOUNDARY_CONDITION_FIXED,ExternalForce,Err)
      ForceNode=12
      IF(DisplacementInterpolationType==CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION) ForceNode=40
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1, &
        & ForceNode,3,CMISS_BOUNDARY_CONDITION_FIXED,ExternalForce,Err)
      IF(DisplacementInterpolationType==CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION) THEN
        ForceNode=45
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1, &
        & ForceNode,3,CMISS_BOUNDARY_CONDITION_FIXED,ExternalForce,Err)
      ENDIF
    ELSE
      ForceNode=TopNodes(1)+NumberGlobalXElements+(NumberGlobalXElements+1)*NumberGlobalYElements/2
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1, &
        & ForceNode,3,CMISS_BOUNDARY_CONDITION_FIXED,ExternalForce,Err)
    ENDIF
  ELSE
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 8,1,CMISS_BOUNDARY_CONDITION_FIXED,0.00593_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 11,1,CMISS_BOUNDARY_CONDITION_FIXED,0.00593_CMISSDP,Err)
  !  CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
  !    & 2,1,CMISS_BOUNDARY_CONDITION_FIXED,-0.028_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 5,1,CMISS_BOUNDARY_CONDITION_FIXED,-0.028_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 9,1,CMISS_BOUNDARY_CONDITION_FIXED,-0.00065_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 12,1,CMISS_BOUNDARY_CONDITION_FIXED,-0.00065_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 3,1,CMISS_BOUNDARY_CONDITION_FIXED,-0.036_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 6,1,CMISS_BOUNDARY_CONDITION_FIXED,-0.036_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 8,2,CMISS_BOUNDARY_CONDITION_FIXED,0.016_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 11,2,CMISS_BOUNDARY_CONDITION_FIXED,-0.016_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 2,2,CMISS_BOUNDARY_CONDITION_FIXED,-0.00016_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 5,2,CMISS_BOUNDARY_CONDITION_FIXED,0.00016_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 9,2,CMISS_BOUNDARY_CONDITION_FIXED,0.000043_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 12,2,CMISS_BOUNDARY_CONDITION_FIXED,-0.000043_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 3,2,CMISS_BOUNDARY_CONDITION_FIXED,-0.000031_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 6,2,CMISS_BOUNDARY_CONDITION_FIXED,0.000031_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 8,3,CMISS_BOUNDARY_CONDITION_FIXED,-0.164_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 11,3,CMISS_BOUNDARY_CONDITION_FIXED,-0.164_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 2,3,CMISS_BOUNDARY_CONDITION_FIXED,-0.161_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 5,3,CMISS_BOUNDARY_CONDITION_FIXED,-0.161_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 9,3,CMISS_BOUNDARY_CONDITION_FIXED,-0.2496_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 12,3,CMISS_BOUNDARY_CONDITION_FIXED,-0.2496_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 3,3,CMISS_BOUNDARY_CONDITION_FIXED,-0.2465_CMISSDP,Err)
    CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
      & 6,3,CMISS_BOUNDARY_CONDITION_FIXED,-0.2465_CMISSDP,Err)
    CALL CMISSBoundaryConditions_SetElement(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,4, &
      & CMISS_BOUNDARY_CONDITION_FIXED,-3846.15_CMISSDP,err)
    CALL CMISSBoundaryConditions_SetElement(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,2,4, &
      & CMISS_BOUNDARY_CONDITION_FIXED,-3846.15_CMISSDP,err)
    ForceNode=12
  ENDIF
  !    & TopNodes(node_idx),3,CMISS_BOUNDARY_CONDITION_FIXED,ExternalForce,Err)
  !ENDDO
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
!  PRINT *,"Top surface nodes: ",TopNodes
!  PRINT *,"Point force applied to node number ",ForceNode
  CALL CMISSProblem_Solve(Problem,Err)
  CALL CMISSField_ParameterSetGetNode(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,ForceNode, &
    & 3,u_y,err)
  PRINT *,"Calculated u_y = ",u_y-Height,' m'

  !Output solution
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"Cantilever","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"Cantilever","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP
  END SUBROUTINE HANDLE_ERROR

END PROGRAM CANTILEVEREXAMPLE

