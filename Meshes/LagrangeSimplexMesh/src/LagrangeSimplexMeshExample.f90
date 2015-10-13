!> \file
!> \author Chris Bradley
!> \brief This is an example program which sets up a field which uses a mixed Lagrange and Simplex mesh using OpenCMISS calls.
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

!> \example LagrangeSimplexMesh/src/LagrangeSimplexMeshExample.f90
!! Example program which sets up a field which uses a mixed Lagrange and Simplex mesh using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/LagrangeSimplexMesh/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/LagrangeSimplexMesh/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM LagrangeSimplexMeshExample

  USE OpenCMISS
  USE OpenCMISS_Iron

  IMPLICIT NONE

  !Test program parameters

  !Program types
  
  !Program variables
  TYPE(cmfe_BasisType) :: basis1,basis2
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem,worldCoordinates
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_MeshElementsType) :: meshElements
  TYPE(cmfe_NodesType) :: nodes
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_FieldType) :: geometricField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_RegionType) :: region,worldRegion
     
  !Generic CMISS variables
  INTEGER(CMISSIntg) :: err
 
  !Intialise cmiss
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_CoordinateSystem_Initialise(worldCoordinates,err)
  CALL cmfe_Initialise(worldCoordinates,WorldRegion,err)
  
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(1,coordinateSystem,Err)
  !Set the coordinate system to be 2D
  CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,2,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,Err)
 
  !Start the creation of a region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(1,worldRegion,region,err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(region,err)

  !Start the creation of a linear-quadratic Lagrange basis
  CALL cmfe_Basis_Initialise(basis1,err)
  CALL cmfe_Basis_CreateStart(1,basis1,err)
  !Set the basis to be a 2D basis
  CALL cmfe_Basis_NumberOfXiSet(basis1,2,err)
  !Set the interpolation to be linear-quadratic
  CALL cmfe_Basis_InterpolationXiSet(basis1,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(basis1,err)

  !Start the creation of a quadratic Simplex triangle basis
  CALL cmfe_Basis_Initialise(basis2,err)
  CALL cmfe_Basis_CreateStart(2,basis2,err)
  !Set the basis to be of Simplex type
  CALL cmfe_Basis_TypeSet(basis2,CMFE_BASIS_SIMPLEX_TYPE,err)
  !Set the basis to be a triangluar basis
  CALL cmfe_Basis_NumberOfXiSet(basis2,2,err)
  !Set the interpolation to be quadratic
  CALL cmfe_Basis_InterpolationXiSet(basis2,[CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION],err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(basis2,err)

  !Create a mesh. The mesh will consist of a linear-quadratic Lagrange element and a quadratic Simplex element i.e.,
  !
  !  8-------9
  !  |       |\
  !  |       | \
  !  |       |  \
  !  5   e1  6   7
  !  |       |    \
  !  |       | e2  \
  !  |       |      \
  !  1-------2---3---4
  !
  CALL cmfe_Nodes_Initialise(nodes,err)
  CALL cmfe_Nodes_CreateStart(region,9,nodes,err)
  CALL cmfe_Nodes_CreateFinish(nodes,err)
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_MeshElements_Initialise(meshElements,err)
  CALL cmfe_Mesh_CreateStart(1,region,2,mesh,err)
  CALL cmfe_Mesh_NumberOfElementsSet(mesh,2,err)
  CALL cmfe_MeshElements_CreateStart(mesh,1,basis1,meshElements,err)
  CALL cmfe_MeshElements_NodesSet(meshElements,1,[1,2,5,6,8,9],err)
  CALL cmfe_MeshElements_BasisSet(meshElements,2,basis2,err)
  CALL cmfe_MeshElements_NodesSet(meshElements,2,[9,2,4,6,3,7],err)
  CALL cmfe_MeshElements_CreateFinish(meshElements,err)
  CALL cmfe_Mesh_CreateFinish(mesh,err)

  !Create a decomposition for mesh
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_CreateStart(1,mesh,decomposition,err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(decomposition,1,err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)
 
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(geometricField,Err)
  CALL cmfe_Field_CreateStart(1,region,geometricField,err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(geometricField,decomposition,err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(geometricField,Err)

  !Set the geometric field values
  !X values
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,1,1,0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,2,1,1.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,3,1,1.5_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,4,1,2.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,5,1,0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,6,1,1.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,7,1,1.5_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,8,1,0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,9,1,1.0_CMISSRP,err)
  !Y values
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,1,2,0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,2,2,0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,3,2,0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,4,2,0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,5,2,0.5_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,6,2,0.5_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,7,2,0.5_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,8,2,1.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1,9,2,1.0_CMISSRP,err)

  CALL cmfe_Field_ParameterSetUpdateStart(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_Field_ParameterSetUpdateFinish(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,err)
 
  !Export results
  CALL cmfe_Fields_Initialise(fields,err)
  CALL cmfe_Fields_Create(region,fields,err)
  CALL cmfe_Fields_NodesExport(fields,"LagrangeSimplexMeshExample","FORTRAN",err)
  CALL cmfe_Fields_ElementsExport(fields,"LagrangeSimplexMeshExample","FORTRAN",err)
  CALL cmfe_Fields_Finalise(fields,err)

  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM LagrangeSimplexMeshExample
