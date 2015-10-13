# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# This python file is a library of python functions which provides modular
# common set-up commands for solving a problem in OpenCMISS. 

# Each function has a range of input options and calls the appropriate
# OpenCMISS linked commands to set up the problem. This is a high 
# level library that will allow shorter scripting for solving cardiac mechanics
# simulations and also making it easier to debug.  

# Author: Zhinuo Jenny Wang
# Start Date: 20th October 2014

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

from opencmiss import iron
import numpy
import math

# =================================================================================#
def BasicSetUp(regionUserNumber, coordinateSystemUserNumber):
    # This function sets up the world region, 3D CS, parallel computing nodes, and
    # diagnostics. 

    # Set up diagnostics/debug
    #iron.DiagnosticsSetOn(iron.DiagnosticTypes.IN,[1,2,3,4,5],
                           #"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

    # Get computational node information for parallel computing
    numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
    computationalNodeNumber = iron.ComputationalNodeNumberGet()

    # Set up 3D RC coordinate system
    coordinateSystem = iron.CoordinateSystem()
    coordinateSystem.CreateStart(coordinateSystemUserNumber)
    coordinateSystem.dimension = 3
    coordinateSystem.CreateFinish()

    # Create world region
    region = iron.Region()
    region.CreateStart(regionUserNumber, iron.WorldRegion)
    region.label = "Region"
    region.coordinateSystem = coordinateSystem
    region.CreateFinish()

    # Output for diagnostics
    print "----> Set up coordinate system and world region <----\n"

    return numberOfComputationalNodes, computationalNodeNumber, coordinateSystem, region


# =================================================================================#

#=================================================================================#
def BasisFunction(basisUserNumber, numOfXi, option, collapsed):
    # This function sets up the basis function depending on the option given.
    if option[0] == 1:
        # Trilinear basis function for interpolation of geometry.
        basis = iron.Basis()
        basis.CreateStart(basisUserNumber)
        basis.numberOfXi = numOfXi
        basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
        basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE] * numOfXi
        basis.QuadratureLocalFaceGaussEvaluateSet(True)
        basis.quadratureNumberOfGaussXi = [2,2,2]
        basis.CreateFinish()
        # Output for diagnostics
        print "----> Set up trilinear basis functions for geometry, use element based interpolation for pressure <----\n"
        if collapsed:
            basisCol = iron.Basis()
            basisCol.CreateStart(basisUserNumber+1)
            basisCol.numberOfXi = numOfXi
            basisCol.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
            basisCol.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE] * numOfXi
            basisCol.QuadratureLocalFaceGaussEvaluateSet(True)
            basisCol.quadratureNumberOfGaussXi = [2,2,2]
            basisCol.CollapsedXiSet([iron.BasisXiCollapse.XI_COLLAPSED, iron.BasisXiCollapse.COLLAPSED_AT_XI0, iron.BasisXiCollapse.NOT_COLLAPSED])
            print "---> Set up collapsed basis functions for apical elements"
            basisCol.CreateFinish()
            return basis, basisCol
        return basis
    elif option[0] == 2:
        quadBasis = iron.Basis()
        quadBasis.CreateStart(basisUserNumber[0])
        quadBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*numOfXi)
        quadBasis.QuadratureNumberOfGaussXiSet([4]*numOfXi)
        quadBasis.QuadratureLocalFaceGaussEvaluateSet(True)
        quadBasis.CreateFinish()

        # Tricubic Hermite basis function for interpolation of geometry.
        cubicBasis = iron.Basis()  # For geometry.
        cubicBasis.CreateStart(basisUserNumber[1])
        cubicBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.CUBIC_HERMITE] * numOfXi)
        cubicBasis.QuadratureNumberOfGaussXiSet([4] * numOfXi)
        cubicBasis.QuadratureLocalFaceGaussEvaluateSet(True)
        cubicBasis.CreateFinish()
        # Output for diagnostics
        print "----> Set up tricubic hermite basis function for geometry and trilinear for hydrostatic pressure <----\n"
        return quadBasis, cubicBasis


#=================================================================================#

#=================================================================================#
def GeneratedMesh(generatedMeshUserNumber, meshUserNumber, region, bases, dimensions, elements):
    # This function sets up a generated mesh using user specified dimensions. 
    generatedMesh = iron.GeneratedMesh()
    generatedMesh.CreateStart(generatedMeshUserNumber, region)
    generatedMesh.TypeSet(iron.GeneratedMeshTypes.REGULAR)
    generatedMesh.BasisSet(bases)
    generatedMesh.ExtentSet(dimensions)
    generatedMesh.NumberOfElementsSet(elements)
    mesh = iron.Mesh()
    generatedMesh.CreateFinish(meshUserNumber, mesh)

    return generatedMesh, mesh


#=================================================================================#

#=================================================================================#
def DecompositionSetUp(decompositionUserNumber, mesh, numberOfComputationalNodes):
    # This function sets up the decomposition of the mesh. 
    decomposition = iron.Decomposition()
    decomposition.CreateStart(decompositionUserNumber, mesh)
    decomposition.type = iron.DecompositionTypes.CALCULATED
    decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
    decomposition.CalculateFacesSet(True)
    decomposition.CreateFinish()

    # Output for diagnostics
    print "----> Set up decomposition <----\n"
    return decomposition


#=================================================================================#

#=================================================================================#
def GeometricFieldSetUp(geometricFieldUserNumber, region, decomposition, option):
    # Set up geometry field 
    geometricField = iron.Field()  # Initialise
    geometricField.CreateStart(geometricFieldUserNumber, region)
    geometricField.MeshDecompositionSet(decomposition)
    geometricField.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")

    if option[0] == 2:
        # Tricubic Hermite
        if option[1] == 1:
            geometricField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
            # Output for diagnostics
            print "----> Set up tricubic Hermite geometric field with unit scaling <----\n"
        elif option[1] == 2:
            geometricField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
            # Output for diagnostics
            print "----> Set up tricubic Hermite geometric field with arithmetic mean scaling <----\n"

    geometricField.CreateFinish()

    return geometricField


#=================================================================================#

#=================================================================================#
def GeometricFieldInitialise(xNodes, yNodes, zNodes, geometricField, numNodes, option):
    # This function initialises the geometric field with user specified coordinates. 
    # Initialise nodal values.
    for node, value in enumerate(xNodes, 1):
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 1, value)
    for node, value in enumerate(yNodes, 1):
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 2, value)
    for node, value in enumerate(zNodes, 1):
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 3, value)

    # Initialise first derivatives.
    if option[0] == 2:
        # Tricubic Hermite basis. 
        for node in range(numNodes):
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, node + 1, 1, max(xNodes))
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, node + 1, 2, max(yNodes))
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3, node + 1, 3, max(zNodes))

    # Output
    print "----> Initialised geometric nodal values <----\n"

    return geometricField


#=================================================================================#

#=================================================================================#
def GeometricFieldExport(region, filename):
    # This function exports the undeformed geometric field. 
    exportField = iron.Fields()
    exportField.CreateRegion(region)
    exportField.NodesExport(filename, "FORTRAN")
    exportField.ElementsExport(filename, "FORTRAN")
    exportField.Finalise()

    # Output
    print "----> Export undeformed geometry <----\n"

#=================================================================================#

#=================================================================================#
def ExtractNodesElements(filename):
    # This function extracts nodes and element connectivity information from 
    # exnode and exelem files. 
    try:
        fid_node = open(filename+'.exnode', 'r')
    except IOError:
        print 'ERROR: Unable to open '+filename+'.exnode'
        return

    try:
        fid_elem = open(filename+'.exelem', 'r')
    except IOError:
        print 'ERROR: Unable to open '+filename+'.exelem'
        return

    for i in range(1,86):
        junk = fid_elem.readline()

    nodesX = []
    nodesY = []
    nodesZ = []
    elements = []
    for i in [1,2,3,4,5,6]:
        junk = fid_node.readline()

    # Read nodal information.
    i = 0
    temp = fid_node.readline()
    while temp != '':
        currentNode = temp.split()[1]
        temp = fid_node.readline()
        nodesX.append(temp.split())
        temp = fid_node.readline()
        nodesY.append(temp.split())
        temp = fid_node.readline()
        nodesZ.append(temp.split())
        i = i+1
        temp = fid_node.readline()
    nodesX = numpy.array(nodesX)
    nodesY = numpy.array(nodesY)
    nodesZ = numpy.array(nodesZ)
    nodes = [nodesX, nodesY, nodesZ]
    nodes = numpy.array(nodes)

    # Read element connectivity
    temp = fid_elem.readline()
    #print temp.split()[0]
    while temp.split() != []:
        currentElem = temp.split()[1]
        junk = fid_elem.readline()
        temp = fid_elem.readline()
        elements.append(temp.split())
        junk = fid_elem.readline()
        junk = fid_elem.readline()
        temp = fid_elem.readline()
    elements = numpy.array(elements)

    return nodes, elements


#=================================================================================#

#=================================================================================#
def FibreFieldSetUp(fibreFieldUserNumber, region, decomposition, geometricField, option, microstructure, inputNodes):
    # This function sets up the fibre field and initialises the values. 
    # Sets up the fibre field.
    fibreField = iron.Field()
    fibreField.CreateStart(fibreFieldUserNumber, region)
    fibreField.TypeSet(iron.FieldTypes.FIBRE)
    fibreField.MeshDecompositionSet(decomposition)
    fibreField.GeometricFieldSet(geometricField)
    fibreField.VariableLabelSet(iron.FieldVariableTypes.U, "Fibre")

    if option[0] == 1:
        fibreField.NumberOfVariablesSet(1)
        fibreField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 3)
        if microstructure == 1:
            for component in [1, 2, 3]:
                fibreField.ComponentInterpolationSet(iron.FieldVariableTypes.U, component,
                                                     iron.FieldInterpolationTypes.CONSTANT)
        elif microstructure == 2:
            for component in [1, 2, 3]:
                fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, component, 1)


    elif option[0] == 2:
        # Tricubic Hermite interpolation
        if option[1] == 1:
            fibreField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
            # Output
            print "----> Set up tricubic hermite fibre field with unit scaling <----\n"
        elif option[1] == 2:
            fibreField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
            # Output
            print "----> Set up tricubic hermite fibre field with arithmetic mean scaling <----\n"

        if microstructure == 1:
            # Homogeneous fibre field.
            for component in [1, 2, 3]:
                fibreField.ComponentInterpolationSet(iron.FieldVariableTypes.U, component,
                                                     iron.FieldInterpolationTypes.CONSTANT)
        elif microstructure == 2:
            # Heterogeneous fibre field using linear interpolation.
            for component in [1, 2, 3]:
                fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, component, 1)

    fibreField.CreateFinish()

    #####################################################
    if microstructure == 2:
        # Inhomogeneous fibre field using linear interpolation.
        for n in range(1, inputNodes.num_nodes+1):
            for component in [1,2,3]:
                component_name = ["x","y","z"][component-1]
                angle = inputNodes.node_values("fibers", component_name, n)
                angle = float(angle[0])
                angle = angle*math.pi/180
                fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                                    1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, n,
                                                    component, angle)
        print "----> Initialised heterogeneous fibre angles <----\n"
    return fibreField


#=================================================================================#

#=================================================================================#
def MaterialFieldSetUpAuto(materialFieldUserNumber, equationsSet, params, cellMLOption):
    # This function is used for setting up material field when using CellML
    # description of constitutive model. 
    # Sets up material field, and apply field to mesh component.
    materialField = iron.Field()
    equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
    materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
    if cellMLOption[0]:
        print "----> CellML Material Field using gauss point interpolation <----\n"
        for component, param in enumerate(params, 1):
            materialField.ComponentInterpolationSet(iron.FieldVariableTypes.U, component,
                                                    iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    materialField.CreateFinish()

    #########################################################################
    # Initialise parameter values.
    for component, param in enumerate(params, 1):
        materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                                  component, param)

    # Output
    print "----> Initialised " + str(len(params)) + " material parameters <----\n"
    return materialField, equationsSet


#=================================================================================#

#=================================================================================#
def MaterialFieldSetUp(materialFieldUserNumber, region, decomposition, geometricField, params, option, cellMLOption):
    # Sets up material field, and apply field to mesh component.
    materialField = iron.Field()
    materialField.CreateStart(materialFieldUserNumber, region)
    materialField.TypeSet(iron.FieldTypes.MATERIAL)
    materialField.MeshDecompositionSet(decomposition)
    materialField.GeometricFieldSet(geometricField)
    materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
    materialField.NumberOfVariablesSet(1)
    materialField.NumberOfComponentsSet(iron.FieldVariableTypes.U,len(params))
    materialField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)

    if cellMLOption[0]:
        print "----> CellML Material Field using gauss point interpolation <----\n"
        for component, param in enumerate(params, 1):
            materialField.ComponentInterpolationSet(iron.FieldVariableTypes.U, component,
                                                    iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    else:
        print "----> Material Field using constant interpolation <----\n"
        for component, param in enumerate(params, 1):
            materialField.ComponentInterpolationSet(iron.FieldVariableTypes.U, component,
                                                    iron.FieldInterpolationTypes.CONSTANT)

    if option[0] == 2:
        # Tricubic Hermite
        if option[1] == 1:
            materialField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
        elif option[1] == 2:
            materialField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    materialField.CreateFinish()

    for component, param in enumerate(params, 1):
        materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                                  component, param)


    materialField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    materialField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    # Output
    print "----> Initialised " + str(len(params)) + " material parameters <----\n"
    return materialField


#=================================================================================#

#=================================================================================#
def DependentFieldSetUp(dependentFieldUserNumber, equationsSet, option, cellMLOption):
    # Set up dependent field
    dependentField = iron.Field()
    equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
    dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")

    if cellMLOption[0]:
        print '----> Labelling dependent field strain and stress <----\n'
        dependentField.VariableLabelSet(iron.FieldVariableTypes.U1, "Strain")
        dependentField.VariableLabelSet(iron.FieldVariableTypes.U2, "Stress")


    if option[0] == 1:
        # Trilinear
        for i in [1, 2, 3]:
            dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, i, 1)
            dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, i, 1)

        dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 4,
                                                 iron.FieldInterpolationTypes.ELEMENT_BASED)
        dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN, 4,
                                                 iron.FieldInterpolationTypes.ELEMENT_BASED)
        # Output
        print "----> Use element based interpolation for hydrostatic pressure <----\n"

    elif option[0] == 2:
        # Tricubic Hermite
        for i in [1, 2, 3]:
            dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, i, 1)
            dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, i, 1)
        dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 4, 2)
        dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 4, 2)

        dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 4,
                                                 iron.FieldInterpolationTypes.NODE_BASED)
        dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN, 4,
                                                 iron.FieldInterpolationTypes.NODE_BASED)

        # Output
        print "----> Interpolate hydrostatic pressure linearly <----\n"
        if option[1] == 1:
            dependentField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
            # Output
            print "----> Set up dependent field with unit scaling <----\n"
        elif option[1] == 2:
            dependentField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
            # Output
            print "----> Set up dependent field with arithmetic mean scaling <----\n"

    if cellMLOption[0]:
        for i in [1,2,3,4,5,6]:
            dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1, i, 1)
            dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2, i, 1)

    equationsSet.DependentCreateFinish()

    return dependentField, equationsSet


#=================================================================================#

#=================================================================================#
def DependentFieldInitialise(dependentField, geometricField, hydroInit):
    # This function initialises the dependent field with reference geometry and 
    # initial guess for hydrostatic pressure. 
    # Copy over undeformed geometry to initialise dependent field.
    for i in [1, 2, 3]:
        iron.Field.ParametersToFieldParametersComponentCopy(geometricField, iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES, i, dependentField,
                                                             iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES, i)

    # Output
    print "----> Initialised dependent field with undeformed geometry <----\n"

    # Set hydrostatic pressure initial guess.
    dependentField.ComponentValuesInitialise(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 4,
                                             hydroInit)

    dependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    dependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    # Output
    print "----> Initialised hydrostatic pressure guess of " + str(hydroInit) + " <----\n"


#=================================================================================#

#=================================================================================#
def DependentFieldWarmStart(dependentField, deformedGeomDOFs, deformedHydro, option):
    # This function reads in warm-start solution to the dependent field. 
    if option[0] == 1:
        # Trilinear elements
        # deformedGeomDOFs indices: component, node
        numNodes = len(deformedGeomDOFs[0,:])
        for component in [1,2,3]:
            for node in range(1, numNodes+1):
                value = deformedGeomDOFs[component-1, node-1]
                dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                                        1, 1, node, component, value)
        for e in range(1, len(deformedHydro)+1):
            dependentField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, e,
                                                    4, deformedHydro[e-1])
    elif option[0] == 2:
        # Initialise dependent field to deformed warmstart solution.
        numNodes = len(deformedGeomDOFs[0,:,0])
        for component in [1,2,3]:
            print 'Component: ', component
            for node in range(1, numNodes+1):
                print '  Node number: ', node
                for deriv in [1,2,3,4,5,6,7,8]:
                    value = deformedGeomDOFs[component-1,node-1,deriv-1]
                    print '    value: ', value
                    dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                                            deriv, node, component, value)


        # Output
        print "----> Initialised dependent field with warm-start geometry <----\n"


        # Set hydrostatic pressure initial guess.
        for node in range(1,numNodes+1):
            dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1,
                                                     1, node, 4, deformedHydro[node-1])

    dependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    dependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    # Output
    print "----> Initialised warm-start hydrostatic pressure of " + str(deformedHydro) + " <----\n"


#=================================================================================#

#=================================================================================#
def ParseWarmStart(filename, option):
    # Read in warmstart solutions from ipinit or exnode files.
    temp = filename.split('.')[1]
    nodesX = []
    nodesY = []
    nodesZ  = []
    if option[0] == 2:
        if temp == 'ipinit':
            try:
                fid = open(filename, 'r')

            except IOError:
                print 'ERROR: Unable to open ', filename
                return

            try:
                junk = fid.readline()
                toggle = True
                while toggle:
                    if junk == " Dependent variable initial conditions:\n":
                        toggle = False
                    else:
                        junk = fid.readline()
                junk = fid.readline()
                junk = fid.readline()

                # Read DOF's in X direction
                temp = fid.readline()
                temp = temp.split()
                num = temp[len(temp)-1]
                while num != '0':
                    derivs = []
                    for i in [0,1,2,3,4,5,6,7]:
                        temp = fid.readline()
                        temp = temp.split()
                        derivs.append(float(temp[len(temp)-1]))
                    nodesX.append(derivs)
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
                nodesX = numpy.array(nodesX)

                # Read DOF's in Y direction
                junk = fid.readline()
                junk = fid.readline()
                temp = fid.readline()
                temp = temp.split()
                num = temp[len(temp)-1]
                while num != '0':
                    derivs = []
                    for i in [0,1,2,3,4,5,6,7]:
                        temp = fid.readline()
                        temp = temp.split()
                        derivs.append(float(temp[len(temp)-1]))
                    nodesY.append(derivs)
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
                nodesY = numpy.array(nodesY)

                # Read DOF's in Y direction
                junk = fid.readline()
                junk = fid.readline()
                temp = fid.readline()
                temp = temp.split()
                num = temp[len(temp)-1]
                while num != '0':
                    derivs = []
                    for i in [0,1,2,3,4,5,6,7]:
                        temp = fid.readline()
                        temp = temp.split()
                        derivs.append(float(temp[len(temp)-1]))
                    nodesZ.append(derivs)
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
                nodesZ = numpy.array(nodesZ)

                # The indices for nodes goes: component (x,y,z), node number, derivative number.
                nodes = [nodesX, nodesY, nodesZ]
                nodes = numpy.array(nodes)

                # Read hydrostatic pressure at nodes
                junk = fid.readline()
                junk = fid.readline()

                node_idx = 0
                temp = fid.readline()
                temp = temp.split()
                num = temp[len(temp)-1]
                hydro = []
                while num!='0':
                    temp = fid.readline()
                    temp = temp.split()
                    hydro.append(float(temp[len(temp)-1]))
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
                hydro = numpy.array(hydro)
                print hydro
                return nodes, hydro
            finally:
                fid.close()
        elif temp == 'exnode':
            nodes = ExtractNodesElements(filename.split('.')[1])
            return nodes
    elif option[0] == 1:
        # Trilinear solution
        if temp == 'ipinit':
            try:
                fid = open(filename, 'r')

            except IOError:
                print 'ERROR: Unable to open ', filename
                return

            try:
                junk = fid.readline()
                toggle = True
                while toggle:
                    if junk == " Dependent variable initial conditions:\n":
                        toggle = False
                    else:
                        junk = fid.readline()
                junk = fid.readline()
                junk = fid.readline()

                # Read DOF's in X direction
                temp = fid.readline()
                temp = temp.split()
                num = temp[len(temp)-1]
                while num != '0':
                    temp = fid.readline()
                    temp = temp.split()
                    nodesX.append(float(temp[len(temp)-1]))
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
                nodesX = numpy.array(nodesX)

                # Read DOF's in Y direction
                junk = fid.readline()
                junk = fid.readline()
                temp = fid.readline()
                temp = temp.split()
                num = temp[len(temp)-1]
                while num != '0':
                    temp = fid.readline()
                    temp = temp.split()
                    nodesY.append(float(temp[len(temp)-1]))
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
                nodesY = numpy.array(nodesY)

                # Read DOF's in Y direction
                junk = fid.readline()
                junk = fid.readline()
                temp = fid.readline()
                temp = temp.split()
                num = temp[len(temp)-1]
                while num != '0':
                    temp = fid.readline()
                    temp = temp.split()
                    nodesZ.append(float(temp[len(temp)-1]))
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
                nodesZ = numpy.array(nodesZ)

                # The indices for nodes goes: component (x,y,z), node number, derivative number.
                nodes = [nodesX, nodesY, nodesZ]
                nodes = numpy.array(nodes)

                # Read hydrostatic pressure at nodes
                junk = fid.readline()
                junk = fid.readline()

                node_idx = 0
                temp = fid.readline()
                temp = temp.split()
                num = temp[len(temp)-1]
                hydro = []
                while num!='0':
                    temp = fid.readline()
                    temp = temp.split()
                    hydro.append(float(temp[len(temp)-1]))
                    temp = fid.readline()
                    temp = temp.split()
                    while (temp[0] != 'Enter'):
                        temp = fid.readline()
                        temp = temp.split()
                    num = temp[len(temp)-1]
                hydro = numpy.array(hydro)
                print hydro
                return nodes, hydro
            finally:
                fid.close()
        elif temp == 'exnode':
            nodes = ExtractNodesElements(filename.split('.')[1])
            return nodes

#=================================================================================#

#=================================================================================#
def CellMLSetUp(cellMLUserNumber, cellMLModelsFieldUserNumber, cellMLParametersFieldUserNumber,
                cellMLIntermediateFieldUserNumber, region, materialField, dependentField, parameters, filename, option):

    # This function sets up the CellML environment for defining constitutive models. 
    cellMLModelIndex = 1
    cellML = iron.CellML()
    cellML.CreateStart(cellMLUserNumber, region)
    cellML.ModelImport(filename)
    strain = ["E11", "E12", "E13", "E22", "E23", "E33"]
    stress2PK = ["Tdev11", "Tdev12", "Tdev13", "Tdev22", "Tdev23", "Tdev33"]

    # Set strains as known in CellML. These will be fed into the model from iron.
    for i in range(0, 6):
        cellML.VariableSetAsKnown(cellMLModelIndex, "equations/" + strain[i])
    for component, parameter in enumerate(parameters):
        cellML.VariableSetAsKnown(cellMLModelIndex, "equations/" + parameter)

    # Set stresses as unknown/wanted in CellML. These will be calculated using the transversely isotropic constitutive model
    for i in range(0, 6):
        cellML.VariableSetAsWanted(cellMLModelIndex, "equations/" + stress2PK[i])
    cellML.CreateFinish()

    # ## Step 13: Map the variables to CellML model ###################################
    cellML.FieldMapsCreateStart()

    # Map the strain from dependentField U1 variable to CellML.
    for component, variable in enumerate(strain, 1):
        #print "----> Mapping strain ", str(variable)+ " to CellML <----\n"
        cellML.CreateFieldToCellMLMap(dependentField, iron.FieldVariableTypes.U1, component,
                                      iron.FieldParameterSetTypes.VALUES, cellMLModelIndex, "equations/" + variable,
                                      iron.FieldParameterSetTypes.VALUES)

    # Map the material parameters from material field to CellML.
    for component, parameter in enumerate(parameters, 1):
        #print "----> Mapping parameter ", str(parameter)+ " to CellML <----\n"
        cellML.CreateFieldToCellMLMap(materialField, iron.FieldVariableTypes.U, component,
                                      iron.FieldParameterSetTypes.VALUES, cellMLModelIndex, "equations/" + parameter,
                                      iron.FieldParameterSetTypes.VALUES)

    # Map the stress from CellML to dependentFieldU2 variable
    for component, variable in enumerate(stress2PK, 1):
        #print "----> Mapping stress ", str(variable)+ " to CellML <----\n"
        cellML.CreateCellMLToFieldMap(cellMLModelIndex, "equations/" + variable, iron.FieldParameterSetTypes.VALUES,
                                      dependentField, iron.FieldVariableTypes.U2, component,
                                      iron.FieldParameterSetTypes.VALUES)
    cellML.FieldMapsCreateFinish()
    print "----> Finished mapping variables to CellML <----\n"

    # Create models field for CellML
    CellMLModelsField = iron.Field()
    cellML.ModelsFieldCreateStart(cellMLModelsFieldUserNumber, CellMLModelsField)
    if option[0] == 2:
        # Tricubic Hermite
        if option[1] == 1:
            CellMLModelsField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
        elif option[1] == 2:
            CellMLModelsField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    cellML.ModelsFieldCreateFinish()
    print "----> Finished creating models field for CellML <----\n"
    # No need to create a state field since we aren't integrating.

    # Create parameters field for CellML, this is used as the strain field.
    CellMLParametersField = iron.Field()
    cellML.ParametersFieldCreateStart(cellMLParametersFieldUserNumber, CellMLParametersField)
    if option[0] == 2:
        # Tricubic Hermite
        if option[1] == 1:
            CellMLParametersField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
        elif option[1] == 2:
            CellMLParametersField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    cellML.ParametersFieldCreateFinish()
    print "----> Finished creating parameters field for CellML <----\n"
    # Create intermediate field for CellML, this is used as the stress field.
    CellMLIntermediateField = iron.Field()
    cellML.IntermediateFieldCreateStart(cellMLIntermediateFieldUserNumber, CellMLIntermediateField)
    if option[0] == 2:
        # Tricubic Hermite
        if option[1] == 1:
            CellMLIntermediateField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
        elif option[1] == 2:
            CellMLIntermediateField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    cellML.IntermediateFieldCreateFinish()
    print "----> Finished creating intermediate field for CellML <----\n"

    return cellML, CellMLModelsField, CellMLParametersField, CellMLIntermediateField
#=================================================================================#

#=================================================================================#
def StrainFieldSetUp(strainFieldUserNumber, region, decomposition, geometricField, equationsSet, option):
    # Set up strain field for output
    strainField = iron.Field()
    strainField.CreateStart(strainFieldUserNumber, region)
    strainField.MeshDecompositionSet(decomposition)
    strainField.TypeSet(iron.FieldTypes.GENERAL)
    strainField.GeometricFieldSet(geometricField)
    strainField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
    strainField.VariableTypesSet([iron.FieldVariableTypes.U])
    strainField.VariableLabelSet(iron.FieldVariableTypes.U, "Strain")
    strainField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 6)

    for component in [1,2,3,4,5,6]:
        strainField.ComponentInterpolationSet(iron.FieldVariableTypes.U, component,
                                              iron.FieldInterpolationTypes.GAUSS_POINT_BASED)

    if option[0]==2:
        if option[1]==1:
            strainField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
        elif option[1]==2:
            strainField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)

    strainField.CreateFinish()

    equationsSet.DerivedCreateStart(strainFieldUserNumber, strainField)
    equationsSet.DerivedVariableSet(iron.EquationsSetDerivedTypes.STRAIN, iron.FieldVariableTypes.U)
    equationsSet.DerivedCreateFinish()

    return strainField
#=================================================================================#

#=================================================================================#
def matrixFromSymmetricComponents(components):
    return numpy.array([
        [components[0], components[1], components[2]],
        [components[1], components[3], components[4]],
        [components[2], components[4], components[5]],
        ])
#=================================================================================#

#=================================================================================#
def EquationsSetSetUp(equationsSet):
    # Set up standard options for problem and solvers.
    # Create equations
    equations = iron.Equations()
    equationsSet.EquationsCreateStart(equations)
    equations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
    equations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
    #equations.OutputTypeSet(iron.EquationsOutputTypes.MATRIX)
    equationsSet.EquationsCreateFinish()

#=================================================================================#

#=================================================================================#
def ProblemSolverSetup(equationsSet,problemUserNumber,maxIter, TOL, cellMLOption):
    # This function sets up the problem as well as the solver options. 
    print "----> Set up equations <----\n"
    # Define problem
    problem = iron.Problem()

    if cellMLOption[0]:
        problemSpecification = [iron.ProblemClasses.ELASTICITY, 
                                iron.ProblemTypes.FINITE_ELASTICITY,
                                iron.ProblemSubtypes.FINITE_ELASTICITY_CELLML]
    else:
        problemSpecification = [iron.ProblemClasses.ELASTICITY, 
                                iron.ProblemTypes.FINITE_ELASTICITY,
                                iron.ProblemSubtypes.NONE]

    problem.CreateStart(problemUserNumber, problemSpecification)
    problem.CreateFinish()
    # Output
    print "----> Set up problem <----\n"
    # Create control loops
    problem.ControlLoopCreateStart()
    controlLoop = iron.ControlLoop()
    problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], controlLoop)
    #controlLoop.TypeSet(iron.ProblemControlLoopTypes.WHILE_LOOP)
    #controlLoop.IterationsSet(1,1,1)
    controlLoop.MaximumIterationsSet(maxIter)
    #controlLoop.MaximumIterationsSet(3)
    problem.ControlLoopCreateFinish()
    # Output
    print "----> Set up control loop <----\n"
    # Create nonlinear numerical solver
    linearSolver = iron.Solver()
    nonLinearSolver = iron.Solver()
    problem.SolversCreateStart()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, nonLinearSolver)
    nonLinearSolver.OutputTypeSet(iron.SolverOutputTypes.PROGRESS)
    nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS)
    nonLinearSolver.NewtonAbsoluteToleranceSet(1e-3)
    nonLinearSolver.NewtonSolutionToleranceSet(1e-2)
    #nonLinearSolver.NewtonRelativeToleranceSet(1e-6)
    nonLinearSolver.NewtonConvergenceTestTypeSet(iron.NewtonConvergenceTypes.PETSC_DEFAULT)
    nonLinearSolver.NewtonLinearSolverGet(linearSolver)
    nonLinearSolver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.LINEAR)
    #nonLinearSolver.NewtonLineSearchAlphaSet(1e-6)
    #nonLinearSolver.NewtonLineSearchMaxStepSet(1e5)
    #nonLinearSolver.NewtonLineSearchMonitorOutputSet()
    #nonLinearSolver.NewtonLineSearchStepTolSet(1e-5)
    linearSolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
    #linearSolver.LinearDirectTypeSet(iron.DirectLinearSolverTypes.LU)
    linearSolver.LibraryTypeSet(iron.SolverLibraries.MUMPS)
    problem.SolversCreateFinish()

    if cellMLOption[0]:
        cellMLSolver = iron.Solver()
        cellMLEquations = iron.CellMLEquations()
        problem.CellMLEquationsCreateStart()
        nonLinearSolver.NewtonCellMLSolverGet(cellMLSolver)
        cellMLSolver.CellMLEquationsGet(cellMLEquations)
        cellMLEquations.CellMLAdd(cellMLOption[1])
        problem.CellMLEquationsCreateFinish()

    # Output
    print "----> Set up linear and nonlinear solvers <----\n"
    # Add solver equations sets which encompass the physics
    solverEquations = iron.SolverEquations()
    solver = iron.Solver()
    problem.SolverEquationsCreateStart()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
    solver.SolverEquationsGet(solverEquations)
    solverEquations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
    equationSetIndex = solverEquations.EquationsSetAdd(equationsSet)
    problem.SolverEquationsCreateFinish()
    # Output
    print "----> Set up solver with equations <----\n"

    return problem, solverEquations
#=================================================================================#

#=================================================================================#
def BCCubeSingleFace(solverEquations, dependentField, appliedFace, faceNormal, appliedDirection, increm, optionBC,
                     fixXFace, fixYFace, fixZFace, numNodes, option):
    # This function sets up the boundary conditions for dealing with BC's on a 
    # single face of a cube. 
    # Set up
    boundaryConditions = iron.BoundaryConditions()
    solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

    # Initialise fixed faces node values.
    for node in fixXFace:
        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1,
                                   iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 1,
                                   iron.BoundaryConditionsTypes.FIXED, 0.0)
    for node in fixYFace:
        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1,
                                   iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 2,
                                   iron.BoundaryConditionsTypes.FIXED, 0.0)
    for node in fixZFace:
        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1,
                                   iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 3,
                                   iron.BoundaryConditionsTypes.FIXED, 0.0)

    if option[0] == 2:

        # Fix derivatives
        if faceNormal == 1:
            derivFix = [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,
                        iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3]
        elif faceNormal == 2:
            derivFix = [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                        iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3]
        else:
            derivFix = [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                        iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2]

        for node in range(1,numNodes+1):
            for j in derivFix:
                for component in [1,2,3]:
                    boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, j, node,
                                               component, iron.BoundaryConditionsTypes.FIXED, 0.0)



        # Fix all second and third derivatives.
        for i in range(1, numNodes + 1):
            for j in [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,
                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,
                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,
                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3]:
                for k in [1, 2, 3]:
                    boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, j, i, k,
                                               iron.BoundaryConditionsTypes.FIXED, 0.0)



    # Output
    print "----> Implemented fixed boundary conditions <----\n"

    # Initialise applied faces.
    if optionBC == 1:
        # Option 1: Compression/extension
        for node in appliedFace:
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1,
                                       iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, appliedDirection,
                                       iron.BoundaryConditionsTypes.FIXED, increm)
        # Output
        print "----> Implemented compression/extension boundary condition of " + str(increm) + " <----\n"

    elif optionBC == 2:
        # Option 2: Force
        for node in appliedFace:
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.DELUDELN, 1,
                                       iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, appliedDirection,
                                       iron.BoundaryConditionsTypes.FIXED_INCREMENTED, increm)
        # Output
        print "----> Implemented force boundary condition of " + str(increm) + "N <----\n"
    elif optionBC == 3:
        # Option 3: Pressure
        for node in appliedFace:
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.DELUDELN, 1,
                                       iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, faceNormal,
                                       iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED, increm)
        # Output
        print "----> Implemented pressure boundary condition of " + str(increm) + " kPa <----\n"

    solverEquations.BoundaryConditionsCreateFinish()


#=================================================================================#

#=================================================================================#
def BCCantilever(solverEquations, dependentField, appliedFace, faceNormal, appliedDirection, increm, optionBC,
                     fixBackFace, fixedFaceNormal, option):
    # This function sets up the BC for a cantilever problem. 
    # Set up
    boundaryConditions = iron.BoundaryConditions()
    solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

    # Initialise fixed faces node values.
    for component in [1, 2, 3]:
        for node in fixBackFace:
            for component in [1,2,3]:
                boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1,
                                           iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, component,
                                           iron.BoundaryConditionsTypes.FIXED, 0.0)
                if option[0] == 2:
                    # print 'Node number ', node

                    # Fix derivatives
                    if fixedFaceNormal == 1:
                        #print "Fixed back normal is 1. "
                        derivFix = [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3]
                    elif fixedFaceNormal == 2:
                        derivFix = [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3]
                    else:
                        derivFix = [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,
                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3]
                    for deriv in derivFix:
                        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, deriv, node,
                                                   component, iron.BoundaryConditionsTypes.FIXED, 0.0)


    # Output
    print "----> Implemented fixed boundary conditions <----\n"

    # Initialise applied faces.
    if optionBC == 1:
        # Option 1: Compression/extension
        for node in appliedFace:
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1,
                                       iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, appliedDirection,
                                       iron.BoundaryConditionsTypes.FIXED, increm)
        # Output
        print "----> Implemented compression/extension boundary condition of " + str(increm) + " <----\n"

    elif optionBC == 2:
        # Option 2: Force
        for node in appliedFace:
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.DELUDELN, 1,
                                       iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, appliedDirection,
                                       iron.BoundaryConditionsTypes.FIXED_INCREMENTED, increm)
        # Output
        print "----> Implemented force boundary condition of " + str(increm) + "N <----\n"
    elif optionBC == 3:
        # Option 3: Pressure
        print 'Pressure applied on: '
        for node in appliedFace:
            print 'Node ', node
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.DELUDELN, 1,
                                       iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, faceNormal,
                                       iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED, increm)
        print 'Face normal ', faceNormal
        # Output
        print "----> Implemented pressure boundary condition of " + str(increm) + " kPa <----\n"

    solverEquations.BoundaryConditionsCreateFinish()


#=================================================================================#

#=================================================================================#
def BCEndoPressure(solverEquations, dependentField, endoFace, pressure, basalFace, option):
    # This function sets up the BC for a LV inflation problem where pressure is applied
    # on the endocardial surface. 
    # Set up
    boundaryConditions = iron.BoundaryConditions()
    solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

    if option[0] == 1:
        derivFix = [iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV]
    else:
        derivFix = [iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,
                iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3]

    # Fix basal nodes and derivatives.
    for component in [1, 2, 3]:
        for node in basalFace:
            for deriv in derivFix:
                boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1,
                                           deriv, node, component,iron.BoundaryConditionsTypes.FIXED, 0.0)


    # Apply pressure BC on endocardial nodes.
    for node in endoFace:
        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.DELUDELN, 1,
                                   iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 3,
                                   iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED, pressure)

    """
    for component in [1,2,3]:
        boundaryConditions.ConstrainNodeDofsEqual(dependentField, iron.FieldVariableTypes.U, 1,
                                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, component, apexEndoNodes)
        boundaryConditions.ConstrainNodeDofsEqual(dependentField, iron.FieldVariableTypes.U, 1,
                                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, component, apexEpiNodes)
    """
    solverEquations.BoundaryConditionsCreateFinish()


#=================================================================================#

#=================================================================================#
def ExportResults(dependentField, deformedFieldUserNumber, decomposition, region, filename, option):
    # This function exports the results of simulation to exnode and exelem files. 
    # Copy over deformed field.
    deformedField = iron.Field()
    deformedField.CreateStart(deformedFieldUserNumber, region)
    deformedField.MeshDecompositionSet(decomposition)
    deformedField.TypeSet(iron.FieldTypes.GEOMETRIC)
    deformedField.VariableLabelSet(iron.FieldVariableTypes.U, "DeformedGeometry")

    if option[0] == 1:
        # Trilinear.
        for component in [1, 2, 3]:
            deformedField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, component, 1)
        deformedField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    elif option[0] == 2:
        # Tricubic hermite. Geometry interpolated using cubic hermite basis (2nd mesh component).
        for component in [1, 2, 3]:
            deformedField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, component, 1)
        if option[1] == 1:
            deformedField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
        elif option[1] == 2:
            deformedField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)

    deformedField.CreateFinish()
    for component in [1, 2, 3]:
        dependentField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES, component,
                                                                deformedField, iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES, component)
    dependentField.Destroy()
    #deformedField.Destroy()
    # Export deformation.
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport(filename, "FORTRAN")
    fields.ElementsExport(filename, "FORTRAN")
    fields.Finalise()

    # Output
    print "----> Export deformed geometric solutions <----\n"


#=================================================================================#

#=================================================================================#
def ExportStressStrain(elements, xiPositions, cellML, equationsSet, filename_disp, filename_strain, filename_stress2PK,
                       filename_stressCauchy, groupname_disp, groupname_strain, groupname_stress):
    # Evaluates the Cauchy strain and 2PK stress at user-specified xi positions
    # for each element.
    # Writes evaluated strain and stress out to exdata file format for visualisation
    # in CMGUI.

    try:
        file_disp = open(filename_disp, 'w')

    except IOError:
        print 'ERROR: Unable to open ', filename_disp
        return

    try:
        file_strain = open(filename_strain, 'w')

    except IOError:
        print 'ERROR: Unable to open ', filename_strain
        return

    try:
        file_stress2PK = open(filename_stress2PK, 'w')

    except IOError:
        print 'ERROR: Unable to open ', filename_stress2PK
        return

    try:
        file_stressCauchy = open(filename_stressCauchy, 'w')

    except IOError:
        print 'ERROR: Unable to open ', filename_stressCauchy
        return

    # Write file headers for displacement
    file_disp.write(' Group name: ' + str(groupname_disp) + '\n')
    file_disp.write(' #Fields= 4\n')
    file_disp.write(' 1) element_xi, field, element_xi, #Components=1\n')
    file_disp.write('   1.  Value index= 1, #Derivatives=0\n')
    file_disp.write(' 2) yg1, field, real, #Components=1\n')
    file_disp.write('   1.  Value index= 1, #Derivatives=0\n')
    file_disp.write(' 3) yg2, field, real, #Components=1\n')
    file_disp.write('   1.  Value index= 1, #Derivatives=0\n')
    file_disp.write(' 4) yg3, field, real, #Components=1\n')
    file_disp.write('   1.  Value index= 1, #Derivatives=0\n')

    # Write file headers for strain
    file_strain.write(' Group name: ' + str(groupname_strain) + '\n')
    file_strain.write(' #Fields= 7\n')
    file_strain.write(' 1) element_xi, field, element_xi, #Components=1\n')
    file_strain.write('   1.  Value index= 1, #Derivatives=0\n')
    file_strain.write(' 2) yg1, field, real, #Components=1\n')
    file_strain.write('   1.  Value index= 1, #Derivatives=0\n')
    file_strain.write(' 3) yg2, field, real, #Components=1\n')
    file_strain.write('   1.  Value index= 1, #Derivatives=0\n')
    file_strain.write(' 4) yg3, field, real, #Components=1\n')
    file_strain.write('   1.  Value index= 1, #Derivatives=0\n')
    file_strain.write(' 5) yg4, field, real, #Components=1\n')
    file_strain.write('   1.  Value index= 1, #Derivatives=0\n')
    file_strain.write(' 6) yg5, field, real, #Components=1\n')
    file_strain.write('   1.  Value index= 1, #Derivatives=0\n')
    file_strain.write(' 7) yg6, field, real, #Components=1\n')
    file_strain.write('   1.  Value index= 1, #Derivatives=0\n')

    # Write file headers for stress
    file_stress2PK.write(' Group name: ' + str(groupname_stress) + '\n')
    file_stress2PK.write(' #Fields= 7\n')
    file_stress2PK.write(' 1) element_xi, field, element_xi, #Components=1\n')
    file_stress2PK.write('   1.  Value index= 1, #Derivatives=0\n')
    file_stress2PK.write(' 2) yg1, field, real, #Components=1\n')
    file_stress2PK.write('   1.  Value index= 1, #Derivatives=0\n')
    file_stress2PK.write(' 3) yg2, field, real, #Components=1\n')
    file_stress2PK.write('   1.  Value index= 1, #Derivatives=0\n')
    file_stress2PK.write(' 4) yg3, field, real, #Components=1\n')
    file_stress2PK.write('   1.  Value index= 1, #Derivatives=0\n')
    file_stress2PK.write(' 5) yg4, field, real, #Components=1\n')
    file_stress2PK.write('   1.  Value index= 1, #Derivatives=0\n')
    file_stress2PK.write(' 6) yg5, field, real, #Components=1\n')
    file_stress2PK.write('   1.  Value index= 1, #Derivatives=0\n')
    file_stress2PK.write(' 7) yg6, field, real, #Components=1\n')
    file_stress2PK.write('   1.  Value index= 1, #Derivatives=0\n')

    # Write file headers for stress
    file_stressCauchy.write(' Group name: ' + str(groupname_stress) + '\n')
    file_stressCauchy.write(' #Fields= 7\n')
    file_stressCauchy.write(' 1) element_xi, field, element_xi, #Components=1\n')
    file_stressCauchy.write('   1.  Value index= 1, #Derivatives=0\n')
    file_stressCauchy.write(' 2) yg1, field, real, #Components=1\n')
    file_stressCauchy.write('   1.  Value index= 1, #Derivatives=0\n')
    file_stressCauchy.write(' 3) yg2, field, real, #Components=1\n')
    file_stressCauchy.write('   1.  Value index= 1, #Derivatives=0\n')
    file_stressCauchy.write(' 4) yg3, field, real, #Components=1\n')
    file_stressCauchy.write('   1.  Value index= 1, #Derivatives=0\n')
    file_stressCauchy.write(' 5) yg4, field, real, #Components=1\n')
    file_stressCauchy.write('   1.  Value index= 1, #Derivatives=0\n')
    file_stressCauchy.write(' 6) yg5, field, real, #Components=1\n')
    file_stressCauchy.write('   1.  Value index= 1, #Derivatives=0\n')
    file_stressCauchy.write(' 7) yg6, field, real, #Components=1\n')
    file_stressCauchy.write('   1.  Value index= 1, #Derivatives=0\n')

    data_pt = 1
    for i in range(0, len(xiPositions)):
        file_disp.write(' Node:         '+str(data_pt)+'\n')
        file_strain.write(' Node:         '+str(data_pt)+'\n')
        file_stress2PK.write(' Node:         '+str(data_pt)+'\n')
        file_stressCauchy.write(' Node:         '+str(data_pt)+'\n')
        file_disp.write(' E         '+str(elements[i])+' 3  '+str(xiPositions[i,0])+'  '+str(xiPositions[i,1])+'  '+str(xiPositions[i,2])+'\n')
        file_strain.write(' E         '+str(elements[i])+' 3  '+str(xiPositions[i,0])+'  '+str(xiPositions[i,1])+'  '+str(xiPositions[i,2])+'\n')
        file_stress2PK.write(' E         '+str(elements[i])+' 3  '+str(xiPositions[i,0])+'  '+str(xiPositions[i,1])+'  '+str(xiPositions[i,2])+'\n')
        file_stressCauchy.write(' E         '+str(elements[i])+' 3  '+str(xiPositions[i,0])+'  '+str(xiPositions[i,1])+'  '+str(xiPositions[i,2])+'\n')
        [disp_temp, strain_temp, stress2PK_temp, stressCauchy_temp] = equationsSet.StrainInterpolateXi(elements[i], xiPositions[i,:], cellML)

        for k in range(0,6):
            file_strain.write('   '+str(strain_temp[k])+'\n')
            file_stress2PK.write('   '+str(stress2PK_temp[k])+'\n')
            file_stressCauchy.write('   '+str(stressCauchy_temp[k])+'\n')
        for m in range(0,3):
            file_disp.write('   '+str(disp_temp[m])+'\n')
        data_pt = data_pt + 1

    file_disp.close()
    file_strain.close()
    file_stress2PK.close()
    file_stressCauchy.close()

    # Output
    print "----> Export stresses and strains of deformed solution <----\n"
