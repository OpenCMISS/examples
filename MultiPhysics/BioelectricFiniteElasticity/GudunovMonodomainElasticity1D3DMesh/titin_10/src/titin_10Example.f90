!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a coupled Monodomain equation Finite Elasticity problem using OpenCMISS calls.
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
!> Contributor(s): Thomas Heidlauf
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

!> \example MultiPhysics/BioelectricFiniteElasticity/GudunovMonodomainElasticity1D3DMesh/simple_geometry/src/STAExample.f90
!! Example program to solve a Monodomain equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/BioelectricFiniteElasticity/GudunovMonodomainElasticity1D3DMesh/simple_geometry/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/BioelectricFiniteElasticity/GudunovMonodomainElasticity1D3DMesh/simple_geometry/build-gnu'>Linux GNU Build</a>
!!
!<

!> Main program
PROGRAM simple_geometryEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !--------------------------------------------------------------------------------------------------------------------------------
  !Test program parameters
!  real etime          ! Declare the type of etime()
  real :: elapsed(2)     ! For receiving user and system time
  real :: total 

!  REAL(CMISSDP), PARAMETER :: P_max=2.7_CMISSDP ! N/cm^2
  REAL(CMISSDP), PARAMETER :: P_max=7.3_CMISSDP ! N/cm^2
!  REAL(CMISSDP), PARAMETER :: P_max=27.0_CMISSDP ! N/cm^2

  REAL(CMISSDP), PARAMETER :: TK_lin_param=0.0_CMISSDP


  REAL(CMISSDP), PARAMETER :: tol=1.0E-8_CMISSDP
  
  LOGICAL :: independent_field_auto_create=.FALSE.

  INTEGER(CMISSIntg), PARAMETER :: NumberOfElementsFE=8

  INTEGER(CMISSIntg) :: NumberOfElementsM
  INTEGER(CMISSIntg) :: NumberOfNodesM
  INTEGER(CMISSIntg) :: NumberOfNodesPerFibre
  INTEGER(CMISSIntg), PARAMETER :: NumberOfInSeriesFibres=1

  integer(CMISSIntg) :: stat
  character(len=256) :: filename,filename2

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: mu_nr,Ftype,fibre_nr,NearestGP,InElement
  
  logical :: less_info,fast_twitch

  
  !--------------------------------------------------------------------------------------------------------------------------------

  !all lengths in [cm]
!  REAL(CMISSDP), PARAMETER :: LENGTH=2.0_CMISSDP ! X-direction
  REAL(CMISSDP), PARAMETER :: LENGTH=1.0_CMISSDP ! X-direction
  REAL(CMISSDP), PARAMETER :: WIDTH= 1.0_CMISSDP ! Y-direction
  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP ! Z-direction
  
  !all times in [ms]
  REAL(CMISSDP) :: time !=10.00_CMISSDP 
!  REAL(CMISSDP), PARAMETER :: PERIODD=10.00_CMISSDP
  REAL(CMISSDP), PARAMETER :: PERIODD=20.00_CMISSDP
!  REAL(CMISSDP), PARAMETER :: TIME_STOP=300.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: TIME_STOP=500.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: TIME_STOP_2=750.0_CMISSDP !time the muscle is stimulated at fixed length + stretch



!  REAL(CMISSDP), PARAMETER :: ODE_TIME_STEP=0.00001_CMISSDP
  REAL(CMISSDP), PARAMETER :: ODE_TIME_STEP=0.0001_CMISSDP
  REAL(CMISSDP), PARAMETER :: PDE_TIME_STEP=0.0005_CMISSDP
  REAL(CMISSDP), PARAMETER :: ELASTICITY_TIME_STEP=0.10000000001_CMISSDP!0.5_CMISSDP!0.05_CMISSDP!0.8_CMISSDP
!tomo keep ELASTICITY_TIME_STEP and STIM_STOP at the same values
  REAL(CMISSDP), PARAMETER :: STIM_STOP=0.1_CMISSDP!ELASTICITY_TIME_STEP

  INTEGER(CMISSIntg), PARAMETER :: OUTPUT_FREQUENCY=1
  
  REAL(CMISSDP) :: stretch_sarcolength_ratio
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !stimulation current in [uA/cm^2]
  REAL(CMISSDP) :: STIM_VALUE

  !condctivity in [mS/cm]
  REAL(CMISSDP), PARAMETER :: CONDUCTIVITY=3.828_CMISSDP

  !surface area to volume ratio in [cm^-1]
  REAL(CMISSDP), PARAMETER :: Am=500.0_CMISSDP

  !membrane capacitance in [uF/cm^2]
  REAL(CMISSDP), PARAMETER :: Cm_fast=1.0_CMISSDP 
  REAL(CMISSDP), PARAMETER :: Cm_slow=0.58_CMISSDP
  

  !maximum contraction velocity in [cm/ms]
  REAL(CMISSDP), PARAMETER :: Vmax=-0.02_CMISSDP ! =0.2 m/s, rat GM
!  REAL(CMISSDP), PARAMETER :: Vmax=-0.2_CMISSDP !to stabilise...
  
  REAL(CMISSDP), PARAMETER :: alpha=5.0_CMISSDP  !1.0_CMISSDP Original Value  !5.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: beta=0.01_CMISSDP
  REAL(CMISSDP), PARAMETER :: gama=1.0_CMISSDP
  !CAUTION - what are the units???
  REAL(CMISSDP), PARAMETER, DIMENSION(4) :: MAT_FE= &
!    & [0.0000000000635201_CMISSDP,0.3626712895523322_CMISSDP,0.0000027562837093_CMISSDP,43.372873938671383_CMISSDP] !original values TOMO [N/cm^2]
!    & [alpha*0.0000000000635201_CMISSDP,alpha*0.3626712895523322_CMISSDP,beta*0.0000027562837093_CMISSDP, &
!      & gama*43.372873938671383_CMISSDP] !new values TOMO [N/cm^2]
    & [alpha*0.0000000000635201_CMISSDP,alpha*0.3626712895523322_CMISSDP,beta*1.074519943356914_CMISSDP, &
      & gama*9.173311371574769_CMISSDP] !new values TOMO [N/cm^2]


  REAL(CMISSDP) :: INIT_PRESSURE

  !--------------------------------------------------------------------------------------------------------------------------------
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumberM=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3

  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumberM=2

  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumberM=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=NumberOfSpatialCoordinates
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussPoints=3

  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensionsFE=3

  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumberM=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponentsFE=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MonodomainMeshComponentNumber=1

  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumberM=2

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumberM=2
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=NumberOfSpatialCoordinates

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumberM=4
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariablesM=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponentsM=3 !Am, Cm, Conductiity   !(scalar, since 1D)

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumberFE=5
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariablesFE=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponentsFE1=5
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponentsFE2=1

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumberM=6
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariablesM=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponentsM=1

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumberFE=7
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariablesFE=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponentsFE=NumberOfSpatialCoordinates+1

  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentUserNumberFE=8
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfVariablesFE=2
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfComponentsFE=2

  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentUserNumberM=9
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfVariablesM=4
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfComponentsM1=1
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfComponentsM2=5

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberM=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberFE=11

  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=15

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetsUserNumberM=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetsUserNumberFE=2

  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: SolverDAEIndex=1
  INTEGER(CMISSIntg), PARAMETER :: SolverParabolicIndex=2
  INTEGER(CMISSIntg), PARAMETER :: SolverFEIndex=1

  INTEGER(CMISSIntg), PARAMETER :: ControlLoopMonodomainNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopElasticityNumber=2
  
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: EquationsSetIndexM,EquationsSetIndexFE
  INTEGER(CMISSIntg) :: CellMLIndex
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx,ComponentNumber,domain_idx,ElementDomain
  INTEGER(CMISSIntg) :: NumberOfNodesInXi1,NumberOfNodesInXi2,NumberOfNodesInXi3
  INTEGER(CMISSIntg) :: i,j,k,my_node_idx,elem_idx,node1,node2,elem_idx2,NumberOfElementsPerElasticityElement

  INTEGER(CMISSIntg), ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg) :: BottomNormalXi,LeftNormalXi,RightNormalXi,FrontNormalXi

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: shortenModelIndex!,shortenModelIndex2
  INTEGER(CMISSIntg) :: stimcomponent

!  REAL(CMISSDP) :: YVALUE
  REAL(CMISSDP) :: VALUE

  INTEGER(CMISSIntg) :: Err


  !CMISS variables

  TYPE(CMISSBasisType) :: QuadraticBasis,LinearBasis,LinearBasisM
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsM,BoundaryConditionsFE
  TYPE(CMISSCellMLType) :: CellML
  TYPE(CMISSCellMLEquationsType) :: CellMLEquations
  TYPE(CMISSControlLoopType) :: ControlLoopMain
  TYPE(CMISSControlLoopType) :: ControlLoopM,ControlLoopFE
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystemFE,CoordinateSystemM,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: DecompositionFE,DecompositionM
  TYPE(CMISSEquationsType) :: EquationsM,EquationsFE
  TYPE(CMISSEquationsSetType) :: EquationsSetM,EquationsSetFE
  TYPE(CMISSFieldType) :: EquationsSetFieldM,EquationsSetFieldFE
  TYPE(CMISSFieldType) :: GeometricFieldM,GeometricFieldFE
  TYPE(CMISSFieldType) :: DependentFieldM,DependentFieldFE
  TYPE(CMISSFieldType) :: IndependentFieldFE,IndependentFieldM
  TYPE(CMISSFieldType) :: MaterialFieldM,MaterialFieldFE
  TYPE(CMISSFieldType) :: FibreField
  TYPE(CMISSFieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSMeshType) :: MeshFE,MeshM
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: RegionFE,RegionM,WorldRegion
  TYPE(CMISSSolverType) :: SolverDAE,SolverParabolic
  TYPE(CMISSSolverType) :: SolverFE,LinearSolverFE
  TYPE(CMISSSolverEquationsType) :: SolverEquationsM,SolverEquationsFE
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSMeshElementsType) :: QuadraticElements
  TYPE(CMISSMeshElementsType) :: LinearElements
  TYPE(CMISSMeshElementsType) :: ElementsM


#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables

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

!##################################################################################################################################
!##################################################################################################################################
!##################################################################################################################################

  NumberGlobalXElements=2
  NumberGlobalYElements=2
  NumberGlobalZElements=2

!##################################################################################################################################
  less_info = .false.!.true.!
  if(less_info) then
    !note that the NumberOfNodesInXi1 only applies to elements in which fibres begin. Otherwise it's NumberOfNodesInXi1-1
    NumberOfNodesInXi1=50!500!240
    NumberOfNodesInXi2=2
    NumberOfNodesInXi3=1
  else
    if(NumberGlobalXElements*NumberGlobalYElements*NumberGlobalZElements==1) then
      NumberOfNodesInXi1=51
    else
 !      NumberOfNodesInXi1=30!500
      NumberOfNodesInXi1=15
    endif
    NumberOfNodesInXi2=3
    NumberOfNodesInXi3=3
!    NumberOfNodesInXi2=4
!    NumberOfNodesInXi3=4
  endif
!  NumberOfNodesPerFibre=(NumberOfNodesInXi1-1)*NumberGlobalXElements+1
!  NumberOfNodesM=NumberOfNodesPerFibre*NumberGlobalYElements*NumberGlobalZElements*NumberOfNodesInXi2*NumberOfNodesInXi3
!  NumberOfElementsM=(NumberOfNodesPerFibre-1)*NumberGlobalYElements*NumberGlobalZElements*NumberOfNodesInXi2*NumberOfNodesInXi3

!  NumberOfNodesPerFibre=NumberOfNodesInXi1
!  NumberOfNodesM=NumberOfNodesPerFibre*NumberGlobalYElements*NumberGlobalZElements*NumberOfNodesInXi2*NumberOfNodesInXi3* &
!    & NumberGlobalXElements
!  NumberOfElementsM=(NumberOfNodesPerFibre-1)*NumberGlobalYElements*NumberGlobalZElements*NumberOfNodesInXi2*NumberOfNodesInXi3* &
!    & NumberGlobalXElements

  NumberOfNodesPerFibre=(NumberOfNodesInXi1-1)*NumberGlobalXElements/NumberOfInSeriesFibres+1
  NumberOfNodesM=NumberOfNodesPerFibre*NumberOfInSeriesFibres*NumberGlobalYElements*NumberGlobalZElements*NumberOfNodesInXi2* &
    & NumberOfNodesInXi3
  NumberOfElementsM=(NumberOfNodesPerFibre-1)*NumberOfInSeriesFibres*NumberGlobalYElements*NumberGlobalZElements* &
    & NumberOfNodesInXi2*NumberOfNodesInXi3

!##################################################################################################################################
!  fast_twitch=.true.
!  if(fast_twitch) then
    filename= &
!   &"data/homes/heidlauf/OpenCMISS/OpenCMISS/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/fast_2012_07_23.xml"
!    "/data/home/heidlauf/OpenCMISS/OpenCMISS/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/fast_2012_07_23.xml"
!    "/home/heidlauf/OpenCMISS/OpenCMISS/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/slow_2014_11_28.xml"
!     "/home/heidlauf/OpenCMISS/OpenCMISS/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/slow_TK_2014_12_08.xml"
     "/home/heidlauf/OpenCMISS/OpenCMISS/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/slow_TK_2015_02_13.xml"
!   &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/shorten_mod_2011_07_04.xml"
    STIM_VALUE=1200.0_CMISSDP
!  else !slow twitch
!    filename2= &
!  &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/fast_stim_2012_07_23.xml"
!  &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/slow_2012_07_23.xml_0.401"
!  &"/home/heidlauf/OpenCMISS/opencmiss/examples/MultiPhysics/BioelectricFiniteElasticity/cellModelFiles/slow_twitch_2012_01_27.xml"
!    STIM_VALUE=2000.0_CMISSDP
!  endif
!##################################################################################################################################


  !================================================================================================================================
  !  G E N E R A L   F E A T U R E S
  !================================================================================================================================

  !--------------------------------------------------------------------------------------------------------------------------------
  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Trap errors
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
  
!  CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,[1,2,3,4,5],"dignostics",["FIELD_MAPPINGS_CALCULATE"],err)
  
  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  CALL CMISSOutputSetOn("simple_geometry",Err)
  
  NumberOfDomains=NumberOfComputationalNodes
  
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystemFE,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumberFE,CoordinateSystemFE,Err)
  !Set the coordinate system to be 3D
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystemFE,3,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystemFE,Err)


  ! CREATE A 1D COORDINATE SYSTEM
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystemM,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumberM,CoordinateSystemM,Err)
  !Set the coordinate system to be 1D
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystemM,1,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystemM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Start the creation of the region
  CALL CMISSRegion_Initialise(RegionFE,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumberFE,WorldRegion,RegionFE,Err)
  CALL CMISSRegion_CoordinateSystemSet(RegionFE,CoordinateSystemFE,Err)
  CALL CMISSRegion_LabelSet(RegionFE,"Region3D",Err)
  CALL CMISSRegion_CreateFinish(RegionFE,Err)


  ! CREATE A SECOND REGION
  CALL CMISSRegion_Initialise(RegionM,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumberM,WorldRegion,RegionM,Err)
  CALL CMISSRegion_CoordinateSystemSet(RegionM,CoordinateSystemM,Err)
  CALL CMISSRegion_LabelSet(RegionM,"Region1D",Err)
  CALL CMISSRegion_CreateFinish(RegionM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the bases
  !Define basis functions - tri-Quadratic Lagrange 
  CALL CMISSBasis_Initialise(QuadraticBasis,Err)
  CALL CMISSBasis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasis_TypeSet(QuadraticBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(QuadraticBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(QuadraticBasis,[CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
   & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
   & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL CMISSBasis_CreateFinish(QuadraticBasis,Err)

  !Define basis functions - tri-Linear Lagrange
  CALL CMISSBasis_Initialise(LinearBasis,Err)
  CALL CMISSBasis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasis_TypeSet(LinearBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(LinearBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(LinearBasis,[CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
   & CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearBasis, &
   & [NumberOfGaussPoints,NumberOfGaussPoints,NumberOfGaussPoints],Err)
  CALL CMISSBasis_CreateFinish(LinearBasis,Err)


  ! CREATE A SECOND LINEAR BASIS FOR THE 1D GRID
  !Define basis functions - tri-Linear Lagrange
  CALL CMISSBasis_Initialise(LinearBasisM,Err)
  CALL CMISSBasis_CreateStart(LinearBasisUserNumberM,LinearBasisM,Err)
  CALL CMISSBasis_TypeSet(LinearBasisM,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(LinearBasisM,1,Err)
  CALL CMISSBasis_InterpolationXiSet(LinearBasisM,[CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearBasisM,[2],Err)
  CALL CMISSBasis_CreateFinish(LinearBasisM,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a mesh with 8 three-dimensional elements
  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,RegionFE,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,[QuadraticBasis,LinearBasis],Err)
  !Define the mesh on the region
  CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH,WIDTH,HEIGHT],Err)
  CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements],Err)
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(MeshFE,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumberFE,MeshFE,Err)
  


  ! CREATE A SECOND MESH
  !Create a mesh in the region
  CALL CMISSMesh_Initialise(MeshM,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumberM,RegionM,1,MeshM,Err) ! 1=NumberOfMeshDimensionsM
  CALL CMISSMesh_NumberOfComponentsSet(MeshM,1,Err) ! 1=NumberOfComponentsM
  CALL CMISSMesh_NumberOfElementsSet(MeshM,NumberOfElementsM,Err)

  CALL CMISSMeshElements_Initialise(ElementsM,Err)
  CALL CMISSMeshElements_CreateStart(MeshM,MonodomainMeshComponentNumber,LinearBasisM,ElementsM,Err)

  !Define nodes for the mesh
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(RegionM,NumberOfNodesM,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)

  elem_idx=0
  DO node1=1,NumberOfNodesM
    IF(mod(node1,NumberOfNodesPerFibre)==0) CYCLE
    elem_idx=elem_idx+1
    CALL CMISSMeshElements_NodesSet(ElementsM,elem_idx,[node1,node1+1],Err)
!    WRITE(*,*) elem_idx,node1,node1+1
  ENDDO    
  write(*,*) "Finished setting up 1D elements"


  CALL CMISSMeshElements_CreateFinish(ElementsM,Err)
  CALL CMISSMesh_CreateFinish(MeshM,Err) 



  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a decomposition
  CALL CMISSDecomposition_Initialise(DecompositionFE,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumberFE,MeshFE,DecompositionFE,Err)

  IF(NumberOfDomains>1) THEN
    CALL CMISSDecomposition_TypeSet(DecompositionFE,CMISS_DECOMPOSITION_USER_DEFINED_TYPE,Err)
    elem_idx2=0
    DO domain_idx=0,NumberOfDomains-1
      DO elem_idx=1,NumberOfElementsFE/NumberOfDomains
        elem_idx2=elem_idx2+1
        CALL CMISSDecomposition_ElementDomainSet(DecompositionFE,elem_idx2,domain_idx,Err)
      ENDDO
    ENDDO       
    CALL CMISSDecomposition_NumberOfDomainsSet(DecompositionFE,NumberOfDomains,Err)
  ELSE
    CALL CMISSDecomposition_TypeSet(DecompositionFE,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL CMISSDecomposition_NumberOfDomainsSet(DecompositionFE,NumberOfDomains,Err)
  ENDIF
  CALL CMISSDecomposition_CalculateFacesSet(DecompositionFE,.TRUE.,Err)
  CALL CMISSDecomposition_CreateFinish(DecompositionFE,Err)

  ! CREATE A SECOND DECOMPOSITION
  CALL CMISSDecomposition_Initialise(DecompositionM,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumberM,MeshM,DecompositionM,Err)
  IF(NumberOfDomains>1) THEN
    CALL CMISSDecomposition_TypeSet(DecompositionM,CMISS_DECOMPOSITION_USER_DEFINED_TYPE,Err)
    elem_idx2=0
    DO domain_idx=0,NumberOfDomains-1
      DO elem_idx=1,NumberOfElementsFE/NumberOfDomains/2
        DO i=1,NumberOfNodesInXi3
          DO j=1,NumberOfNodesInXi2
            DO k=1,NumberOfNodesPerFibre-1
              elem_idx2=elem_idx2+1
              CALL CMISSDecomposition_ElementDomainSet(DecompositionM,elem_idx2,domain_idx,Err)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    IF(elem_idx2/=NumberOfElementsM) THEN
      WRITE(*,*) "Error in setting up the decomposition for monodomain!"
      STOP
    ENDIF

    CALL CMISSDecomposition_NumberOfDomainsSet(DecompositionM,NumberOfDomains,Err)
  ELSE
    CALL CMISSDecomposition_TypeSet(DecompositionM,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL CMISSDecomposition_NumberOfDomainsSet(DecompositionM,NumberOfDomains,Err)
  ENDIF
  CALL CMISSDecomposition_CreateFinish(DecompositionM,Err)


  !================================================================================================================================
  !  F I N I T E   E L A S T C I T Y
  !================================================================================================================================
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for finite elasticity - quadratic interpolation
  CALL CMISSField_Initialise(GeometricFieldFE,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumberFE,RegionFE,GeometricFieldFE,Err)
  CALL CMISSField_TypeSet(GeometricFieldFE,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricFieldFE,DecompositionFE,Err)
  CALL CMISSField_NumberOfVariablesSet(GeometricFieldFE,FieldGeometryNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_VariableLabelSet(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL CMISSField_CreateFinish(GeometricFieldFE,Err)

!  CALL CMISSField_ParameterSetUpdateStart(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
!  CALL CMISSField_ParameterSetUpdateFinish(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricFieldFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a fibre field and attach it to the geometric field - quadratic interpolation
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FieldFibreUserNumber,RegionFE,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,DecompositionFE,Err)
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricFieldFE,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err) 
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_VariableLabelSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL CMISSField_CreateFinish(FibreField,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a material field for Finite Elasticity and attach it to the geometric field - quadratic interpolation
  CALL CMISSField_Initialise(MaterialFieldFE,Err)
  CALL CMISSField_CreateStart(FieldMaterialUserNumberFE,RegionFE,MaterialFieldFE,Err)
  CALL CMISSField_TypeSet(MaterialFieldFE,CMISS_FIELD_MATERIAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(MaterialFieldFE,DecompositionFE,Err)
  CALL CMISSField_GeometricFieldSet(MaterialFieldFE,GeometricFieldFE,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialFieldFE,FieldMaterialNumberOfVariablesFE,Err)
  CALL CMISSField_VariableTypesSet(MaterialFieldFE,[CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_V_VARIABLE_TYPE],Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,6,Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,FieldMaterialNumberOfComponentsFE2,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,3,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,5,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,6,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,1,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_VariableLabelSet(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,"MaterialFE",Err)
  CALL CMISSField_VariableLabelSet(MaterialFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,"Gravity",Err)
  CALL CMISSField_CreateFinish(MaterialFieldFE,Err)
  !Set Mooney-Rivlin constants c10 and c01.
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,MAT_FE(1),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,MAT_FE(2),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,MAT_FE(3),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,MAT_FE(4),Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,5,0.0_CMISSDP, &
   & Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,6,TK_lin_param,&
   & Err)
!  CALL CMISSField_ComponentValuesInitialise(MaterialFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)


  !Create the dependent field for FE with 2 variables and * components 
  !  3-d: 3 displacement (quad interpol), 1 pressure (lin interpol) --> * = 4
  !  2-d: 2 displacement (quad interpol), 1 pressure (lin interpol) --> * = 3
  !  1-d: 1 displacement (quad interpol), 1 pressure (lin interpol) --> * = 2
  CALL CMISSField_Initialise(DependentFieldFE,Err)
  CALL CMISSField_CreateStart(FieldDependentUserNumberFE,RegionFE,DependentFieldFE,Err)
  CALL CMISSField_TypeSet(DependentFieldFE,CMISS_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentFieldFE,DecompositionFE,Err)
  CALL CMISSField_GeometricFieldSet(DependentFieldFE,GeometricFieldFE,Err)
  CALL CMISSField_DependentTypeSet(DependentFieldFE,CMISS_FIELD_DEPENDENT_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentFieldFE,FieldDependentNumberOfVariablesFE,Err)
  CALL CMISSField_VariableTypesSet(DependentFieldFE,[CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE],Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponentsFE,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponentsFE,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
!  CALL CMISSField_ScalingTypeSet(DependentFieldFE,CMISS_FIELD_UNIT_SCALING,Err)
  CALL CMISSField_VariableLabelSet(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,"DependentFE",Err)
  CALL CMISSField_VariableLabelSet(DependentFieldFE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,"Reaction_Force",Err)
  CALL CMISSField_CreateFinish(DependentFieldFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the independent field for the active stress in the mechanics mesh
!  independent_field_auto_create = .TRUE.
  CALL CMISSField_Initialise(IndependentFieldFE,Err)
  IF(.NOT. independent_field_auto_create) THEN
    CALL CMISSField_CreateStart(FieldIndependentUserNumberFE,RegionFE,IndependentFieldFE,Err)
    CALL CMISSField_TypeSet(IndependentFieldFE,CMISS_FIELD_GENERAL_TYPE,Err)
    CALL CMISSField_MeshDecompositionSet(IndependentFieldFE,DecompositionFE,Err)
    CALL CMISSField_GeometricFieldSet(IndependentFieldFE,GeometricFieldFE,Err)
    CALL CMISSField_DependentTypeSet(IndependentFieldFE,CMISS_FIELD_INDEPENDENT_TYPE,Err)
    CALL CMISSField_NumberOfVariablesSet(IndependentFieldFE,FieldIndependentNumberOfVariablesFE,Err)
    CALL CMISSField_VariableTypesSet(IndependentFieldFE,[CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_V_VARIABLE_TYPE],Err)
    CALL CMISSField_DimensionSet(IndependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VECTOR_DIMENSION_TYPE,Err)
    CALL CMISSField_NumberOfComponentsSet(IndependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,5,Err)
    CALL CMISSField_NumberOfComponentsSet(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,4,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1, &
     & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,2, &
     & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! titin force (unbound)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,3, &
     & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! titin force (bound)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,4, &
     & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! titin force in XF-direction (unbound)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,5, &
     & CMISS_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! titin force in XF-direction (bound)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,1, &
     & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,2, &
     & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,3, &
     & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,4, &
     & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentMeshComponentSet(IndependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(IndependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
    CALL CMISSField_DataTypeSet(IndependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DP_TYPE,Err)
    CALL CMISSField_DataTypeSet(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_INTG_TYPE,Err)
    CALL CMISSField_VariableLabelSet(IndependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,"Active_Stress_FE",Err)
    CALL CMISSField_VariableLabelSet(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,"subgrid_info",Err)
    CALL CMISSField_CreateFinish(IndependentFieldFE,Err)
  ENDIF


  !================================================================================================================================
  !  M O N O D O M A I N
  !================================================================================================================================
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for monodomain - quadratic interpolation
  CALL CMISSField_Initialise(GeometricFieldM,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumberM,RegionM,GeometricFieldM,Err)
  CALL CMISSField_TypeSet(GeometricFieldM,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricFieldM,DecompositionM,Err)
  CALL CMISSField_TypeSet(GeometricFieldM,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(GeometricFieldM,1,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMISSField_VariableLabelSet(GeometricFieldM,CMISS_FIELD_U_VARIABLE_TYPE,"GeometryM",Err)
  CALL CMISSField_CreateFinish(GeometricFieldM,Err)



  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a materials field for monodomain and attach it to the geometric field - constant interpolation
  CALL CMISSField_Initialise(MaterialFieldM,Err)
  CALL CMISSField_CreateStart(FieldMaterialUserNumberM,RegionM,MaterialFieldM,Err)
  CALL CMISSField_TypeSet(MaterialFieldM,CMISS_FIELD_MATERIAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(MaterialFieldM,DecompositionM,Err)
  CALL CMISSField_GeometricFieldSet(MaterialFieldM,GeometricFieldM,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialFieldM,FieldMaterialNumberOfVariablesM,Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponentsM,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,3,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_VariableLabelSet(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,"MaterialM",Err)
  CALL CMISSField_CreateFinish(MaterialFieldM,Err)
  !Set Am
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Am,Err)
  !Set Cm
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Cm_fast,Err)
  !Set Conductivity
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3, &
   & CONDUCTIVITY,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the dependent field for monodomain with 2 variables and 1 components 
  CALL CMISSField_Initialise(DependentFieldM,Err)
  CALL CMISSField_CreateStart(FieldDependentUserNumberM,RegionM,DependentFieldM,Err)
  CALL CMISSField_TypeSet(DependentFieldM,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentFieldM,DecompositionM,Err)
  CALL CMISSField_GeometricFieldSet(DependentFieldM,GeometricFieldM,Err)
  CALL CMISSField_DependentTypeSet(DependentFieldM,CMISS_FIELD_DEPENDENT_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentFieldM,FieldDependentNumberOfVariablesM,Err)
  CALL CMISSField_VariableTypesSet(DependentFieldM,[CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE, &
   & CMISS_FIELD_V_VARIABLE_TYPE],Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponentsM,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldM,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponentsM,Err)
  !additional the V_Var_Type with 3 components for the 3-d position of the geometry
  CALL CMISSField_NumberOfComponentsSet(DependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,3,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,MonodomainMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldM,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,MonodomainMeshComponentNumber,Err)
!  CALL CMISSField_ComponentMeshComponentSet(DependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,1,MonodomainMeshComponentNumber,Err)
  CALL CMISSField_DimensionSet(DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMISSField_DimensionSet(DependentFieldM,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMISSField_VariableLabelSet(DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,"Vm",Err)
  CALL CMISSField_VariableLabelSet(DependentFieldM,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,"dVm/dt",Err)
  CALL CMISSField_VariableLabelSet(DependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,"GeometryM3D",Err)
  CALL CMISSField_CreateFinish(DependentFieldM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the independent field for the active stress in the electrics mesh
!  independent_field_auto_create = .TRUE.
  CALL CMISSField_Initialise(IndependentFieldM,Err)
  IF(.NOT. independent_field_auto_create) THEN
    CALL CMISSField_CreateStart(FieldIndependentUserNumberM,RegionM,IndependentFieldM,Err)
    CALL CMISSField_TypeSet(IndependentFieldM,CMISS_FIELD_GENERAL_TYPE,Err)
    CALL CMISSField_MeshDecompositionSet(IndependentFieldM,DecompositionM,Err)
    CALL CMISSField_GeometricFieldSet(IndependentFieldM,GeometricFieldM,Err)
    CALL CMISSField_DependentTypeSet(IndependentFieldM,CMISS_FIELD_INDEPENDENT_TYPE,Err)
    CALL CMISSField_NumberOfVariablesSet(IndependentFieldM,FieldIndependentNumberOfVariablesM,Err) !4
    CALL CMISSField_VariableTypesSet(IndependentFieldM,[CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_V_VARIABLE_TYPE, &
     & CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_U2_VARIABLE_TYPE],Err)
    
    !first variable:   CMISS_FIELD_U_VARIABLE_TYPE -- 1) active stress   2) initial titin stress P_0 (or passive sarcomere stress when non activated)    3) titin stress delta P    4) titin stress in cross-fibre direction
    CALL CMISSField_DataTypeSet(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DP_TYPE,Err)
    CALL CMISSField_DimensionSet(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VECTOR_DIMENSION_TYPE,Err)
!    CALL CMISSField_NumberOfComponentsSet(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,Err)
    CALL CMISSField_NumberOfComponentsSet(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,5,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,2, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the initial titin stress P_0 (or passive sarcomere stress when non activated)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,3, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the additional titin stress delta P
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,4, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the titin stress in the XF-direction
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,5, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err) !this component is for the titin stress in the XF-direction
    CALL CMISSField_VariableLabelSet(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,"Active_Stress_M",Err)

    !second variable:   CMISS_FIELD_V_VARIABLE_TYPE -- 1) motor unit number   2) fibre type   3) fibre number   4) nearest Gauss point   5) in element number (LOCAL NODE NUMBERING!!!)
    CALL CMISSField_DataTypeSet(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_INTG_TYPE,Err)
    CALL CMISSField_NumberOfComponentsSet(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,FieldIndependentNumberOfComponentsM2,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,1, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,2, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,3, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,4, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,5, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_VariableLabelSet(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,"fibre_info",Err)

    !third variable:   FIELD_U1_VARIABLE_TYPE -- 1) sarcomere half length   2) inital sarcomere half length   3) initial node distance   4) sarcomere half length when activation starts
    CALL CMISSField_DataTypeSet(IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_DP_TYPE,Err)
    CALL CMISSField_NumberOfComponentsSet(IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,4,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,1, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,2, &
     & CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,3, &
     & CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,4, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_VariableLabelSet(IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,"sarcomere_half_length",Err)

    !fourth variable:   FIELD_U2_VARIABLE_TYPE -- 1) old node distance   2) maximum contraction velocity   3) relative contraction velocity
    CALL CMISSField_DataTypeSet(IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_DP_TYPE,Err)
    CALL CMISSField_NumberOfComponentsSet(IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,3,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,1, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,2, &
     & CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,3, &
     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_VariableLabelSet(IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,"contraction velocity",Err)

    CALL CMISSField_CreateFinish(IndependentFieldM,Err)
  ENDIF
  

  
  !================================================================================================================================
  !  EQUATIONS SET

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for Finite Elasticity
  CALL CMISSField_Initialise(EquationsSetFieldFE,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetFE,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetsUserNumberFE,RegionFE,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
   & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, & 
   & EquationsSetFieldUserNumberFE,EquationsSetFieldFE,EquationsSetFE,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSetFE,Err)

  !Create the equations set dependent field variables for Finite Elasticity
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetFE,FieldDependentUserNumberFE,DependentFieldFE,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetFE,Err)

  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL CMISSEquationsSet_IndependentCreateStart(EquationsSetFE,FieldIndependentUserNumberFE,IndependentFieldFE,Err)
  CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSetFE,Err)

  !Create the equations set materials field variables for Finite Elasticity
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetFE,FieldMaterialUserNumberFE,MaterialFieldFE,Err)  
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for monodomain
  CALL CMISSField_Initialise(EquationsSetFieldM,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetM,Err)
  !Set the equations set to be a Monodomain equations set
  !> \todo solve the monodomain problem on the fibre field rather than on the geometric field: GeometricField <--> FibreField
  CALL CMISSEquationsSet_CreateStart(EquationsSetsUserNumberM,RegionM,GeometricFieldM,CMISS_EQUATIONS_SET_BIOELECTRICS_CLASS, &
   & CMISS_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMISS_EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
   & EquationsSetFieldUserNumberM,EquationsSetFieldM,EquationsSetM,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSetM,Err)

  !Create the equations set dependent field variables for monodomain
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetM,FieldDependentUserNumberM,DependentFieldM,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetM,Err)

  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL CMISSEquationsSet_IndependentCreateStart(EquationsSetM,FieldIndependentUserNumberM,IndependentFieldM,Err)
  CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSetM,Err)

  !Create the equations set materials field variables for monodomain
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetM,FieldMaterialUserNumberM,MaterialFieldM,Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetM,Err)

 
 
  !UPDATE THE INDEPENDENT FIELD IndependentFieldM
  !first variable
  !  components: 
  !    1) active stress
!  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
!    & 0.0_CMISSDP,Err)
!  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3, &
!    & 0.0_CMISSDP,Err)
!  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4, &
!    & 0.0_CMISSDP,Err)
  !
  !second variable
  !  components: 
  !    1) motor unit number
  !    2) fibre type
  !    3) fibre number
  !    4) nearest Gauss point
  !    5) in element number (LOCAL NODE NUMBERING!!!)
  !
  !init the motor unit number and fibre type to 1
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,Err) !mu_nr=1
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,1,Err) !Ftype=1
  !init the fibre number, the nearest Gauss point info and the inElem info to 0
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,0,Err)
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,0,Err)
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,5,0,Err) !(LOCAL NODE NUMBERING!!!)
  !third variable:
  !  components:
  !    1) sarcomere half length
  !    2) initial sarcomere half length
  !    3) initial node distance
  !    4) sarcomere half length when activation starts
  stretch_sarcolength_ratio=1.0_CMISSDP/1.0_CMISSDP
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
!   & 1.0_CMISSDP,Err)
   & stretch_sarcolength_ratio,Err)
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3, &
   & LENGTH/NumberOfInSeriesFibres/(NumberOfNodesPerFibre-1),Err)
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4, &
   & 1.0_CMISSDP,Err)
!   & stretch_sarcolength_ratio,Err)
!   & 5.0_CMISSDP,Err)
  !fourth variable:
  !  components:
  !    1) old node distance
  !    2) maximum contraction velocity
  !    3) relative contraction velocity
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
   & LENGTH/NumberOfInSeriesFibres/(NumberOfNodesPerFibre-1),Err)
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
   & Vmax/NumberOfInSeriesFibres/(NumberOfNodesPerFibre-1),Err)
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3, &
   & 0.0_CMISSDP,Err)


  !UPDATE THE INDEPENDENT FIELD IndependentFieldFE
  !second variable of IndependentFieldFE
  !  components:
  !    1) number of nodes in Xi(1) direction per element
  !    2) number of nodes in Xi(2) direction per element
  !    3) number of nodes in Xi(3) direction per element
  !    4) beginning of fibres in this FE element??? 1=yes, 0=no
  !
  !initialise as if the fibres would not start in any element, and adjust below
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
   & NumberOfNodesInXi1-1,Err)
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
   & NumberOfNodesInXi2,Err)
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3, &
   & NumberOfNodesInXi3,Err)
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4, &
   & 0,Err)

  !fibres are starting in elements 1,4,7,10,...
  DO elem_idx=1,NumberOfElementsFE,NumberGlobalXElements/NumberOfInSeriesFibres
    CALL CMISSDecomposition_ElementDomainGet(DecompositionFE,elem_idx,ElementDomain,Err)
    IF(ElementDomain==ComputationalNodeNumber) THEN
      !fibres begin in this element
      CALL CMISSField_ParameterSetUpdateElement(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
       & elem_idx,4,1,Err) 
      CALL CMISSField_ParameterSetUpdateElement(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
       & elem_idx,1,NumberOfNodesInXi1,Err)
    ENDIF
  ENDDO

!!  CALL CMISSField_ComponentValuesInitialise(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,0,Err)
!  !fibres start in all elements
!  CALL CMISSField_ComponentValuesInitialise(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,1,Err)
!  !adjust numer of nodes in each element in Xi1 direction
!  CALL CMISSField_ComponentValuesInitialise(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
!   & NumberOfNodesInXi1,Err)

!  !fibres are starting only in element 1
!  elem_idx=1
!  CALL CMISSDecomposition_ElementDomainGet(DecompositionFE,elem_idx,ElementDomain,Err)
!  IF(ElementDomain==ComputationalNodeNumber) THEN
!    !fibres begin in this element
!    CALL CMISSField_ParameterSetUpdateElement(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
!     & elem_idx,4,1,Err) 
!    CALL CMISSField_ParameterSetUpdateElement(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
!     & elem_idx,1,NumberOfNodesInXi1,Err)
!  ENDIF

!  !fibres are starting in elements 1,3,5, and 7
!  DO elem_idx=1,NumberOfElementsFE,2
!    CALL CMISSDecomposition_ElementDomainGet(DecompositionFE,elem_idx,ElementDomain,Err)
!    IF(ElementDomain==ComputationalNodeNumber) THEN
!      !fibres begin in this element
!      CALL CMISSField_ParameterSetUpdateElement(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
!       & elem_idx,4,1,Err) 
!      CALL CMISSField_ParameterSetUpdateElement(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
!       & elem_idx,1,NumberOfNodesInXi1,Err)
!    ENDIF
!  ENDDO

!  !fibres are not starting in elements 2,4,6, and 8
!  DO elem_idx=2,NumberOfElementsFE,2
!    CALL CMISSDecomposition_ElementDomainGet(DecompositionFE,elem_idx,ElementDomain,Err)
!    IF(ElementDomain==ComputationalNodeNumber) THEN
!      CALL CMISSField_ParameterSetUpdateElement(IndependentFieldFE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
!       & elem_idx,4,0,Err) !fibres do not begin in this element
!    ENDIF
!  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set equations for monodomain
  CALL CMISSEquations_Initialise(EquationsM,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetM,EquationsM,Err)
  CALL CMISSEquations_SparsityTypeSet(EquationsM,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(EquationsM,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetM,Err)

  !Create the equations set equations for Finite Elasticity
  CALL CMISSEquations_Initialise(EquationsFE,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetFE,EquationsFE,Err)
  CALL CMISSEquations_SparsityTypeSet(EquationsFE,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(EquationsFE,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetFE,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML environment
  CALL CMISSCellML_Initialise(CellML,Err)
  CALL CMISSCellML_CreateStart(CellMLUserNumber,RegionM,CellML,Err)
  !Import the Shorten et al. 2007 model from a file
  CALL CMISSCellML_ModelImport(CellML,filename,shortenModelIndex,Err)
!  CALL CMISSCellML_ModelImport(CellML,filename2,shortenModelIndex2,Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !,- to set from this side
  CALL CMISSCellML_VariableSetAsKnown(CellML,shortenModelIndex,"wal_environment/I_HH",Err)
  CALL CMISSCellML_VariableSetAsKnown(CellML,shortenModelIndex,"razumova/L_S",Err)
!  CALL CMISSCellML_VariableSetAsKnown(CellML,shortenModelIndex,"razumova/rel_velo",Err)
!
!  CALL CMISSCellML_VariableSetAsKnown(CellML,shortenModelIndex2,"wal_environment/I_HH",Err)
!  CALL CMISSCellML_VariableSetAsKnown(CellML,shortenModelIndex2,"razumova/L_S",Err)
!  CALL CMISSCellML_VariableSetAsKnown(CellML,shortenModelIndex2,"razumova/rel_velo",Err)
  !,- to get from the CellML side
!  CALL CMISSCellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_T",Err)
!  CALL CMISSCellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_ionic_s",Err)
!  CALL CMISSCellML_VariableSetAsWanted(CellML,shortenModelIndex,"wal_environment/I_ionic_t",Err)
  !
  !NOTE: If an INTERMEDIATE (or ALGEBRAIC) variable should be used in a mapping, it has to be set as known or wanted first!
  !,  --> set "razumova/stress" as wanted!
  !,  --> no need to set "wal_environment/vS" since all STATE variables are automatically set as wanted! 
  CALL CMISSCellML_VariableSetAsWanted(CellML,shortenModelIndex,"razumova/stress",Err)
!  CALL CMISSCellML_VariableSetAsWanted(CellML,shortenModelIndex2,"razumova/stress",Err)
  !,- and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
  !Finish the CellML environment
  CALL CMISSCellML_CreateFinish(CellML,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML <--> OpenCMISS field maps
  CALL CMISSCellML_FieldMapsCreateStart(CellML,Err)
  !Map the sarcomere half length L_S
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
   & shortenModelIndex,"razumova/L_S",CMISS_FIELD_VALUES_SET_TYPE,Err)
  !Map the sarcomere relative contraction velocity
!  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,3,CMISS_FIELD_VALUES_SET_TYPE, &
!   & shortenModelIndex,"razumova/rel_velo",CMISS_FIELD_VALUES_SET_TYPE,Err)
  !Map the transmembrane voltage V_m
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
   & shortenModelIndex,"wal_environment/vS",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"wal_environment/vS",CMISS_FIELD_VALUES_SET_TYPE, &
   & DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !Map the active stress
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,shortenModelIndex,"razumova/stress",CMISS_FIELD_VALUES_SET_TYPE, &
   & IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

!  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
!   & shortenModelIndex2,"razumova/L_S",CMISS_FIELD_VALUES_SET_TYPE,Err)
!  !Map the sarcomere relative contraction velocity
!  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,3,CMISS_FIELD_VALUES_SET_TYPE, &
!   & shortenModelIndex2,"razumova/rel_velo",CMISS_FIELD_VALUES_SET_TYPE,Err)
!  !Map the transmembrane voltage V_m
!  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
!   & shortenModelIndex2,"wal_environment/vS",CMISS_FIELD_VALUES_SET_TYPE,Err)
!  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,shortenModelIndex2,"wal_environment/vS",CMISS_FIELD_VALUES_SET_TYPE, &
!   & DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)
!  !Map the active stress
!  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,shortenModelIndex2,"razumova/stress",CMISS_FIELD_VALUES_SET_TYPE, &
!   & IndependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

  CALL CMISSCellML_FieldMapsCreateFinish(CellML,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Initialise dependent field for monodomain
  !> \todo - get V_m initialial value.
  CALL CMISSField_ComponentValuesInitialise(DependentFieldM,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
   & -79.974_CMISSDP,Err)
  
  !Initialise dependent field for Finite Elasticity from undeformed geometry and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE, &
   & CMISS_FIELD_VALUES_SET_TYPE,1,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE, & 
   & CMISS_FIELD_VALUES_SET_TYPE,2,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE, &
   & CMISS_FIELD_VALUES_SET_TYPE,3,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  INIT_PRESSURE=-2.0_CMISSDP*MAT_FE(2)-MAT_FE(1)
  CALL CMISSField_ComponentValuesInitialise(DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4, & 
   & INIT_PRESSURE,Err)
!   & -8.0_CMISSDP,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML models field
  CALL CMISSField_Initialise(CellMLModelsField,Err)
  CALL CMISSCellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  CALL CMISSCellML_ModelsFieldCreateFinish(CellML,Err)

  !Set up the models field
  CALL CMISSField_ComponentValuesInitialise(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
   & shortenModelIndex,Err)

!  DO NodeNumber=NumberOfNodesPerFibre/2,NumberOfNodesM,NumberOfNodesPerFibre
!    CALL CMISSDecomposition_NodeDomainGet(DecompositionM,NodeNumber,1,NodeDomain,Err)
!    IF(NodeDomain==ComputationalNodeNumber) THEN
!      CALL CMISSField_ParameterSetUpdateNode(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE, &
!        & CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,shortenModelIndex2,Err)
!    ENDIF
!  ENDDO

!  CALL CMISSField_ParameterSetUpdateStart(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
!  CALL CMISSField_ParameterSetUpdateFinish(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)







  !Create the CellML state field
  CALL CMISSField_Initialise(CellMLStateField,Err)
  CALL CMISSCellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,Err)
  CALL CMISSCellML_StateFieldCreateFinish(CellML,Err)

  !Create the CellML intermediate field
  CALL CMISSField_Initialise(CellMLIntermediateField,Err)
  CALL CMISSCellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  CALL CMISSCellML_IntermediateFieldCreateFinish(CellML,Err)
  
  !Create the CellML parameters field
  CALL CMISSField_Initialise(CellMLParametersField,Err)
  CALL CMISSCellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  CALL CMISSCellML_ParametersFieldCreateFinish(CellML,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Define the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_MULTI_PHYSICS_CLASS,CMISS_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE, &
   & CMISS_PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)

  !set the main control loop (time loop type)
  CALL CMISSControlLoop_Initialise(ControlLoopMain,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoopMain,Err)
  CALL CMISSControlLoop_LabelSet(ControlLoopMain,'MAIN_TIME_LOOP',Err)
  !Loop in time for STIM_STOP with the Stimulus applied.
  CALL CMISSControlLoop_TimesSet(ControlLoopMain,0.0_CMISSDP,ELASTICITY_TIME_STEP,ELASTICITY_TIME_STEP,Err)
  CALL CMISSControlLoop_TimeOutputSet(ControlLoopMain,OUTPUT_FREQUENCY,Err)
  CALL CMISSControlLoop_OutputTypeSet(ControlLoopMain,CMISS_CONTROL_LOOP_TIMING_OUTPUT,Err) !DO NOT CHANGE!!!

  !set the monodomain loop (time loop type)
  CALL CMISSControlLoop_Initialise(ControlLoopM,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,[ControlLoopMonodomainNumber,CMISS_CONTROL_LOOP_NODE],ControlLoopM,Err)
  CALL CMISSControlLoop_LabelSet(ControlLoopM,'MONODOMAIN_TIME_LOOP',Err)
  CALL CMISSControlLoop_TimesSet(ControlLoopM,0.0_CMISSDP,ELASTICITY_TIME_STEP,PDE_TIME_STEP,Err)
  CALL CMISSControlLoop_OutputTypeSet(ControlLoopM,CMISS_CONTROL_LOOP_NO_OUTPUT,Err)

  !set the finite elasticity loop (simple type)
  CALL CMISSControlLoop_Initialise(ControlLoopFE,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,[ControlLoopElasticityNumber,CMISS_CONTROL_LOOP_NODE],ControlLoopFE,Err)
  CALL CMISSControlLoop_TypeSet(ControlLoopFE,CMISS_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,Err)
  CALL CMISSControlLoop_MaximumIterationsSet(ControlLoopFE,20,Err) ! tomo
!  CALL CMISSControlLoop_MaximumIterationsSet(ControlLoopFE,1,Err)
  CALL CMISSControlLoop_LabelSet(ControlLoopFE,'ELASTICITY_LOOP',Err)

  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solvers
  CALL CMISSProblem_SolversCreateStart(Problem,Err)

  !Create the DAE solver
  CALL CMISSSolver_Initialise(SolverDAE,Err)
  CALL CMISSProblem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMISS_CONTROL_LOOP_NODE], &
   & SolverDAEIndex,SolverDAE,Err)
  CALL CMISSSolver_DAETimeStepSet(SolverDAE,ODE_TIME_STEP,Err)
  !> \todo - solve the CellML equations on the GPU for efficiency (later)
  !CALL CMISSSolver_DAESolverTypeSet(SolverDAE,CMISS_SOLVER_DAE_EXTERNAL,Err) 
  CALL CMISSSolver_OutputTypeSet(SolverDAE,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverDAE,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverDAE,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverDAE,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverDAE,CMISS_SOLVER_MATRIX_OUTPUT,Err)

  !Create the parabolic solver
  CALL CMISSSolver_Initialise(SolverParabolic,Err)
  CALL CMISSProblem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMISS_CONTROL_LOOP_NODE], &
   & SolverParabolicIndex,SolverParabolic,Err)
  CALL CMISSSolver_DynamicSchemeSet(SolverParabolic,CMISS_SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME,Err)
  CALL CMISSSolver_OutputTypeSet(SolverParabolic,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverParabolic,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverParabolic,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverParabolic,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverParabolic,CMISS_SOLVER_MATRIX_OUTPUT,Err)

  !Create the Finte Elasticity solver
  CALL CMISSSolver_Initialise(SolverFE,Err)
  CALL CMISSSolver_Initialise(LinearSolverFE,Err)
  CALL CMISSProblem_SolverGet(Problem,[ControlLoopElasticityNumber,CMISS_CONTROL_LOOP_NODE], &
   & SolverFEIndex,SolverFE,Err)
  CALL CMISSSolver_OutputTypeSet(SolverFE,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverFE,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverFE,CMISS_SOLVER_TIMING_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverFE,CMISS_SOLVER_SOLVER_OUTPUT,Err)
  !CALL CMISSSolver_OutputTypeSet(SolverFE,CMISS_SOLVER_MATRIX_OUTPUT,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(SolverFE,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL CMISSSolver_NewtonMaximumIterationsSet(SolverFE,500,Err)
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(SolverFE,1.E-6_CMISSDP,Err) !6
  CALL CMISSSolver_NewtonSolutionToleranceSet(SolverFE,2.E-6_CMISSDP,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(SolverFE,LinearSolverFE,Err)
!  CALL CMISSSolver_LinearTypeSet(LinearSolverFE,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolverFE,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)

  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver CellML equations
  CALL CMISSProblem_CellMLEquationsCreateStart(Problem,Err)

  CALL CMISSSolver_Initialise(SolverDAE,Err)
  CALL CMISSCellMLEquations_Initialise(CellMLEquations,Err)
  CALL CMISSProblem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMISS_CONTROL_LOOP_NODE],SolverDAEIndex,SolverDAE,Err)
  CALL CMISSSolver_CellMLEquationsGet(SolverDAE,CellMLEquations,Err)
  CALL CMISSCellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  CALL CMISSProblem_CellMLEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)

  !Create the problem solver parabolic equations (Monodomain)
  CALL CMISSSolver_Initialise(SolverParabolic,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsM,Err)
  CALL CMISSProblem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMISS_CONTROL_LOOP_NODE], &
   & SolverParabolicIndex,SolverParabolic,Err)
  CALL CMISSSolver_SolverEquationsGet(SolverParabolic,SolverEquationsM,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsM,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsM,CMISS_SOLVER_FULL_MATRICES,Err)  
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsM,EquationsSetM,EquationsSetIndexM,Err)

  !Create the problem solver Finite Elasticity equations
  CALL CMISSSolver_Initialise(SolverFE,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsFE,Err)
  CALL CMISSProblem_SolverGet(Problem,[ControlLoopElasticityNumber,CMISS_CONTROL_LOOP_NODE],SolverFEIndex,SolverFE,Err)
  CALL CMISSSolver_SolverEquationsGet(SolverFE,SolverEquationsFE,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsFE,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsFE,CMISS_SOLVER_FULL_MATRICES,Err)  
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsFE,EquationsSetFE,EquationsSetIndexFE,Err)

  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !boundary conditions

  !Prescribe boundary conditions for monodomain
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsM,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsM,BoundaryConditionsM,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsM,Err)

  !Prescribe boundary conditions for Finite Elasticity (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsFE,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsFE,BoundaryConditionsFE,Err)

  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi, &
    & Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)

  !Set x=0 nodes to no x displacment in x. Set x=WIDTH nodes to 100% x displacement
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO

  DO node_idx=1,SIZE(RightSurfaceNodes,1)
    NodeNumber=RightSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
!        & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED,LENGTH*1.0_CMISSDP,Err) ! To change the initial half-sarcomere length. To create cases similar to Rode. 
        & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED,LENGTH*1.0_CMISSDP/stretch_sarcolength_ratio,Err) ! To change the initial half-sarcomere length. To create cases similar to Rode. 
    ENDIF
  ENDDO

  !Set y=0 nodes to no y displacement
  DO node_idx=1,SIZE(FrontSurfaceNodes,1)
    NodeNumber=FrontSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO

  !Set z=0 nodes to no z displacement
  DO node_idx=1,SIZE(BottomSurfaceNodes,1)
    NodeNumber=BottomSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO

!  !fix node 1 in all directions
!  NodeNumber=1
!  CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
!  IF(NodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
!     & 2,VALUE,Err)
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
!     & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!    CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,&
!     & 3,VALUE,Err)
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
!     & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!  ENDIF

!  !fix nodes 26 and 51 in direction z
!  NodeNumber=(2*NumberGlobalXElements+1)+1
!  CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
!  IF(NodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
!     & 3,VALUE,Err)
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
!     & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!  ENDIF
!  NodeNumber=2*(2*NumberGlobalXElements+1)+1
!  CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
!  IF(NodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
!     & 3,VALUE,Err)
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
!     & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!  ENDIF

!  !fix nodes 76 and 151 in direction y
!  NodeNumber=3*(2*NumberGlobalXElements+1)+1
!  CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
!  IF(NodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
!     & 2,VALUE,Err)
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
!     & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!  ENDIF
!  NodeNumber=6*(2*NumberGlobalXElements+1)+1
!  CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
!  IF(NodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
!     & 2,VALUE,Err)
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
!     & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!  ENDIF

!  !fix node 25 in all directions
!  NodeNumber=(2*NumberGlobalXElements+1)
!  CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
!  IF(NodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
!     & 2,VALUE,Err)
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
!     & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!    CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,&
!     & 3,VALUE,Err)
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
!     & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!  ENDIF

!  !fix nodes 50 and 75 in direction z
!  NodeNumber=2*(2*NumberGlobalXElements+1)
!  CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
!  IF(NodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
!     & 3,VALUE,Err)
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
!     & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!  ENDIF
!  NodeNumber=3*(2*NumberGlobalXElements+1)
!  CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
!  IF(NodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
!     & 3,VALUE,Err)
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
!     & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!  ENDIF

!  !fix nodes 100 and 175 in direction y
!  NodeNumber=4*(2*NumberGlobalXElements+1)
!  CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
!  IF(NodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
!     & 2,VALUE,Err)
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
!     & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!  ENDIF
!  NodeNumber=7*(2*NumberGlobalXElements+1)
!  CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
!  IF(NodeDomain==ComputationalNodeNumber) THEN
!    CALL CMISSField_ParameterSetGetNode(GeometricFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,& 
!     & 2,VALUE,Err)
!    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
!     & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!  ENDIF

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Calculate the bioelectrics geometric field 
  CALL CMISSControlLoop_Initialise(ControlLoopM,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,[ControlLoopMonodomainNumber,CMISS_CONTROL_LOOP_NODE],ControlLoopM,Err)
  CALL CMISSBioelectricsFiniteElasticity_UpdateGeometricField(ControlLoopM,.TRUE.,Err)

  !reset the relative contraction velocity to 0
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3, &
   & 0.0_CMISSDP,Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(RegionM,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"simple_geometryExample_M","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"simple_geometryExample_M","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)

    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(RegionFE,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"simple_geometryExample_FE","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"simple_geometryExample_FE","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF



  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem


  !Solve the problem -- bring to new length before applying the stimulus
  WRITE(*,'(A)') "Start solve before stimulation"
  CALL CMISSProblem_Solve(Problem,Err)

  !reset the relative contraction velocity to 0
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldM,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3, &
   & 0.0_CMISSDP,Err)

  CALL CMISSControlLoop_Initialise(ControlLoopFE,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,[ControlLoopElasticityNumber,CMISS_CONTROL_LOOP_NODE],ControlLoopFE,Err)
  CALL CMISSControlLoop_MaximumIterationsSet(ControlLoopFE,1,Err)

  CALL CMISSField_ParameterSetUpdateConstant(MaterialFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,5, &
   & P_max,Err)


! no change for BCs -- fix at this length!!!




  !Set the Stimulus for monodomain at the middle of the fibres
  CALL CMISSCellML_FieldComponentGet(CellML,shortenModelIndex,CMISS_CELLML_PARAMETERS_FIELD, &
   & "wal_environment/I_HH",stimcomponent,Err)

  !update the sarcomere stretch at activation
  CALL CMISSField_ParametersToFieldParametersComponentCopy(IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE, &
   & CMISS_FIELD_VALUES_SET_TYPE,1,IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,Err)



  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  k = -1
  time = 0.0_CMISSDP 
  !first activate without stretching
  do while(time <= TIME_STOP)

  k = k+1


  NodeNumber=(NumberOfNodesPerFibre+1)/2
  DO WHILE(NodeNumber<NumberOfNodesM)
    CALL CMISSDecomposition_NodeDomainGet(DecompositionM,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField, &
     & CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,STIM_VALUE,Err)
    NodeNumber=NodeNumber+NumberOfNodesPerFibre
  ENDDO
  
  CALL CMISSControlLoop_TimesSet(ControlLoopMain,time,time+STIM_STOP,ELASTICITY_TIME_STEP,Err)

  !Solve the problem for the stimulation time
  WRITE(*,'(A)') "Start solve with stimulation"
  CALL CMISSProblem_Solve(Problem,Err)

  time = time+STIM_STOP

  !--------------------------------------------------------------------------------------------------------------------------------
  !Now turn the stimulus off
  NodeNumber=(NumberOfNodesPerFibre+1)/2
  DO WHILE(NodeNumber<NumberOfNodesM)
    CALL CMISSDecomposition_NodeDomainGet(DecompositionM,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField, &
     & CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,0.0_CMISSDP,Err)
    NodeNumber=NodeNumber+NumberOfNodesPerFibre
  ENDDO

  WRITE(*,'(A)') "Start solve without stimulation"

  do while(time <= (k+1)*PERIODD)

    CALL CMISSControlLoop_TimesSet(ControlLoopMain,time,time+STIM_STOP,ELASTICITY_TIME_STEP,Err)

    !Solve the problem for the rest of the period
    CALL CMISSProblem_Solve(Problem,Err)
  
    time = time+STIM_STOP
    
  end do !time <= PERIODD

!  time = time+PERIODD
  end do !time <= TIME_STOP
  !--------------------------------------------------------------------------------------------------------------------------------


  !update the sarcomere stretch at activation
!  CALL CMISSField_ParametersToFieldParametersComponentCopy(IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE, &
!   & CMISS_FIELD_VALUES_SET_TYPE,1,IndependentFieldM,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,Err)




  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! bis hierher kein aktiver stretch!!!
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------





  !--------------------------------------------------------------------------------------------------------------------------------
  VALUE = 0.0_CMISSDP 

  !second, activate and stretch!!!
  do while(time <= TIME_STOP_2)

  k = k+1


!  VALUE = VALUE-ABS(Vmax)/10.0_CMISSDP*STIM_STOP
  VALUE = VALUE-ABS(Vmax)/10.0_CMISSDP*STIM_STOP*2.0_CMISSDP

  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  ENDDO


  NodeNumber=(NumberOfNodesPerFibre+1)/2
  DO WHILE(NodeNumber<NumberOfNodesM)
    CALL CMISSDecomposition_NodeDomainGet(DecompositionM,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField, &
     & CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,STIM_VALUE,Err)
    NodeNumber=NodeNumber+NumberOfNodesPerFibre
  ENDDO
  
!  !stimulate the fibers one after another
!  NodeNumber=(NumberOfNodesPerFibre+1)/2+k*NumberOfNodesPerFibre
!  DO WHILE(NodeNumber>NumberOfNodesM)
!    NodeNumber = NodeNumber-NumberOfNodesM
!  ENDDO
!  
!  WRITE(*,*) NodeNumber
!  
!  CALL CMISSDecomposition_NodeDomainGet(DecompositionM,NodeNumber,1,NodeDomain,Err)
!  IF(NodeDomain==ComputationalNodeNumber) CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField, &
!   & CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,STIM_VALUE,Err)


  CALL CMISSControlLoop_TimesSet(ControlLoopMain,time,time+STIM_STOP,ELASTICITY_TIME_STEP,Err)

  !Solve the problem for the stimulation time
  WRITE(*,'(A)') "Start solve with stimulation"
  CALL CMISSProblem_Solve(Problem,Err)

  time = time+STIM_STOP


  !--------------------------------------------------------------------------------------------------------------------------------
  !Now turn the stimulus off
  NodeNumber=(NumberOfNodesPerFibre+1)/2
  DO WHILE(NodeNumber<NumberOfNodesM)
    CALL CMISSDecomposition_NodeDomainGet(DecompositionM,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) CALL CMISSField_ParameterSetUpdateNode(CellMLParametersField, &
     & CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,stimcomponent,0.0_CMISSDP,Err)
    NodeNumber=NodeNumber+NumberOfNodesPerFibre
  ENDDO



  WRITE(*,'(A)') "Start solve without stimulation"

  do while(time <= (k+1)*PERIODD)

!    VALUE = VALUE-ABS(Vmax)/10.0_CMISSDP*STIM_STOP
    VALUE = VALUE-ABS(Vmax)/10.0_CMISSDP*STIM_STOP*2.0_CMISSDP
    DO node_idx=1,SIZE(LeftSurfaceNodes,1)
      NodeNumber=LeftSurfaceNodes(node_idx)
      CALL CMISSDecomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
          & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      ENDIF
    ENDDO
    
    CALL CMISSControlLoop_TimesSet(ControlLoopMain,time,time+STIM_STOP,ELASTICITY_TIME_STEP,Err)

    !Solve the problem for the rest of the period
    CALL CMISSProblem_Solve(Problem,Err)
  
    time = time+STIM_STOP
    
  end do !time <= PERIODD

!  time = time+PERIODD
  end do !time <= TIME_STOP
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------


  
  !--------------------------------------------------------------------------------------------------------------------------------
!  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(RegionM,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"simple_geometryExample_M","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"simple_geometryExample_M","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)

    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(RegionFE,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"simple_geometryExample_FE","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"simple_geometryExample_FE","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  if(NumberOfComputationalNodes==1) then
    open(unit=888,file="times.out",iostat=stat)
    write(888,*) 'single processor'
!!! time elapsed
    total = etime(elapsed)
    write(888,*) 'End: total=', total, ' user=', elapsed(1), ' system=', elapsed(2)
!!!
    close(unit=888)
  else
    open(unit=888,file="times.out",iostat=stat,access='append')
    write(888,*) 'two processors'
!!! time elapsed
    total = etime(elapsed)
    write(888,*) 'End: total=', total, ' user=', elapsed(1), ' system=', elapsed(2)
!!!
    close(unit=888)
  endif
  
  STOP
  
END PROGRAM simple_geometryEXAMPLE

