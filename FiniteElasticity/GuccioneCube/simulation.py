#!/usr/bin/env python

# Based on the OpenCMISS-Iron uniaxial extension example, trying to replicate example 524 from classic-cm, using 
# Guccione in CellML.

#> Main script
# Add Python bindings directory to PATH
import sys, os
from numpy import pi

# Intialise OpenCMISS
from opencmiss.iron import iron

def usage(progName):
    print("Usage: " + progName + " <fibre angle>")
    print("\tAngle should be given in degrees.")
    print("\tReplicating the examples 524-527 (0, 30, 60, 90 degree fibre angles) from cm-classic.")

def simulate(fibreAngleIn, materialParameters):
    '''
    Main function to run a finite elasticity simulation using OpenCMISS with a unit cube and the Guccione material law (in CellML).
    fibreAngle is a scalar angle in degrees
    materialParameters is a list of the four Guccione material parameters. 
    '''
    # Set problem parameters - Unit cube
    height = 1.0
    width = 1.0
    length = 1.0
    
#     if len(sys.argv) != 2:
#         usage(sys.argv[0])
#         exit(-1)
    
    fa = fibreAngleIn
    print("fa = " + str(fa))
    
    # Fibre angles in radians:
    fibreAngle = fa * pi / 180.0
    # (transversly isotropic, so assume 0 cross-fibre angle?) 
    sheetAngle = 90.0 * pi / 180.0
    fibreAngles = [fibreAngle, 0.0, sheetAngle]
    
    # Guccione constitutive relation:
    constitutiveRelation = iron.EquationsSetSubtypes.CONSTITUTIVE_LAW_IN_CELLML_EVALUATE
    constitutiveParameters = [0.88, 18.5, 3.58, 3.26]
    constitutiveParameters = materialParameters
    initialHydrostaticPressure = 0.0
    
    UsePressureBasis = False
    NumberOfGaussXi = 2
    
    coordinateSystemUserNumber = 1
    regionUserNumber = 1
    basisUserNumber = 1
    pressureBasisUserNumber = 2
    generatedMeshUserNumber = 1
    meshUserNumber = 1
    decompositionUserNumber = 1
    geometricFieldUserNumber = 1
    fibreFieldUserNumber = 2
    materialFieldUserNumber = 3
    dependentFieldUserNumber = 4
    equationsSetFieldUserNumber = 5
    deformedFieldUserNumber = 6
    equationsSetUserNumber = 1
    problemUserNumber = 1
    
    CellMLUserNumber = 1
    CellMLModelsFieldUserNumber = 7
    CellMLParametersFieldUserNumber = 8
    CellMLIntermediateFieldUserNumber = 9
    
    # Set all diganostic levels on for testing
    #iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])
    
    numberGlobalXElements = 1
    numberGlobalYElements = 1
    numberGlobalZElements = 1
    totalNumberOfNodes=8
    totalNumberOfElements=1
    InterpolationType = 1
    if(UsePressureBasis):
      numberOfMeshComponents = 2
    else:
      numberOfMeshComponents = 1
    if(numberGlobalZElements==0):
        numberOfXi = 2
    else:
        numberOfXi = 3
    
    # Get the number of computational nodes and this computational node number
    numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
    computationalNodeNumber = iron.ComputationalNodeNumberGet()
    
    # Create a 3D rectangular cartesian coordinate system
    coordinateSystem = iron.CoordinateSystem()
    coordinateSystem.CreateStart(coordinateSystemUserNumber)
    coordinateSystem.DimensionSet(3)
    coordinateSystem.CreateFinish()
    
    # Create a region and assign the coordinate system to the region
    region = iron.Region()
    region.CreateStart(regionUserNumber,iron.WorldRegion)
    region.LabelSet("Region")
    region.coordinateSystem = coordinateSystem
    region.CreateFinish()
    
    # Define basis
    basis = iron.Basis()
    basis.CreateStart(basisUserNumber)
    if InterpolationType in (1,2,3,4):
        basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
    elif InterpolationType in (7,8,9):
        basis.type = iron.BasisTypes.SIMPLEX
    basis.numberOfXi = numberOfXi
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
    if(NumberOfGaussXi>0):
        basis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
    basis.CreateFinish()
    
    if(UsePressureBasis):
        # Define pressure basis
        pressureBasis = iron.Basis()
        pressureBasis.CreateStart(pressureBasisUserNumber)
        if InterpolationType in (1,2,3,4):
            pressureBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
        elif InterpolationType in (7,8,9):
            pressureBasis.type = iron.BasisTypes.SIMPLEX
        pressureBasis.numberOfXi = numberOfXi
        pressureBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
        if(NumberOfGaussXi>0):
            pressureBasis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
        pressureBasis.CreateFinish()
    
    # Start the creation of a manually generated mesh in the region
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber,region,numberOfXi)
    mesh.NumberOfComponentsSet(numberOfMeshComponents)
    mesh.NumberOfElementsSet(totalNumberOfElements)
    
    #Define nodes for the mesh
    nodes = iron.Nodes()
    nodes.CreateStart(region,totalNumberOfNodes)
    nodes.CreateFinish()
    
    elements = iron.MeshElements()
    meshComponentNumber=1
    elements.CreateStart(mesh,meshComponentNumber,basis)
    elements.NodesSet(1,[1,2,3,4,5,6,7,8])
    elements.CreateFinish()
    
    mesh.CreateFinish() 
    
    # Create a decomposition for the mesh
    decomposition = iron.Decomposition()
    decomposition.CreateStart(decompositionUserNumber,mesh)
    decomposition.type = iron.DecompositionTypes.CALCULATED
    decomposition.numberOfDomains = numberOfComputationalNodes
    decomposition.CreateFinish()
    
    # Create a field for the geometry
    geometricField = iron.Field()
    geometricField.CreateStart(geometricFieldUserNumber,region)
    geometricField.MeshDecompositionSet(decomposition)
    geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
    geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
    if InterpolationType == 4:
        geometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
    geometricField.CreateFinish()
    
    # Update the geometric field parameters manually
    geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    # node 1
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,1,1,0.0)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,1,2,0.0)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,1,3,0.0)
    # node 2
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,2,1,height)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,2,2,0.0)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,2,3,0.0)
    # node 3
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,3,1,0.0)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,3,2,width)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,3,3,0.0)
    # node 4
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,4,1,height)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,4,2,width)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,4,3,0.0)
    # node 5
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,5,1,0.0)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,5,2,0.0)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,5,3,length)
    # node 6
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,6,1,height)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,6,2,0.0)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,6,3,length)
    # node 7
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,7,1,0.0)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,7,2,width)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,7,3,length)
    # node 8
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,8,1,height)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,8,2,width)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,8,3,length)
    geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    
    # Create a fibre field and attach it to the geometric field
    fibreField = iron.Field()
    fibreField.CreateStart(fibreFieldUserNumber,region)
    fibreField.TypeSet(iron.FieldTypes.FIBRE)
    fibreField.MeshDecompositionSet(decomposition)
    fibreField.GeometricFieldSet(geometricField)
    fibreField.VariableLabelSet(iron.FieldVariableTypes.U,"Fibre")
    if InterpolationType == 4:
        fibreField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
    fibreField.CreateFinish()
    
    iron.Field.ComponentValuesInitialiseDP(fibreField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                           1, fibreAngle)
    iron.Field.ComponentValuesInitialiseDP(fibreField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                           2, 0.0)
    iron.Field.ComponentValuesInitialiseDP(fibreField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                           3, sheetAngle)
    
    # # Create the material field
    # materialField = iron.Field()
    # materialField.CreateStart(materialFieldUserNumber,region)
    # materialField.TypeSet(iron.FieldTypes.MATERIAL)
    # materialField.MeshDecompositionSet(decomposition)
    # materialField.GeometricFieldSet(geometricField)
    # materialField.NumberOfVariablesSet(1)
    # materialField.NumberOfComponentsSet(iron.FieldVariableTypes.U,len(constitutiveParameters))
    # materialField.VariableLabelSet(iron.FieldVariableTypes.U,"Material")
    # materialField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
    # materialField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
    # if InterpolationType == 4:
    #     materialField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
    # materialField.CreateFinish()
    
    # Set constant material parameters:
    # for (component, value) in enumerate(constitutiveParameters, 1):
    #     materialField.ComponentValuesInitialise(
    #             iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
    #             component, value)
        
    # Create the dependent field
    dependentField = iron.Field()
    dependentField.CreateStart(dependentFieldUserNumber,region)
    dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
    dependentField.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)  
    dependentField.MeshDecompositionSet(decomposition)
    dependentField.GeometricFieldSet(geometricField) 
    dependentField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT) 
    dependentField.NumberOfVariablesSet(4)
    dependentField.VariableTypesSet([iron.FieldVariableTypes.U,iron.FieldVariableTypes.DELUDELN,iron.FieldVariableTypes.U1,iron.FieldVariableTypes.U2])
    dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,4)
    dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,4)
    dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U1,6)
    dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U2,6)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)  
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,1,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,2,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,3,1)  
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,1,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,2,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,3,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,4,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,5,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,6,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,1,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,2,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,3,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,4,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,5,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,6,1)
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
    if(UsePressureBasis):
        # Set the pressure to be nodally based and use the second mesh component
        if InterpolationType == 4:
            dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.NODE_BASED)
            dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.NODE_BASED)
        dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
        dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)
    if InterpolationType == 4:
        dependentField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
        
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,1,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,2,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,3,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,4,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,5,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,6,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,1,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,2,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,3,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,4,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,5,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,6,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    
    dependentField.VariableLabelSet(iron.FieldVariableTypes.U1,"strain")    
    dependentField.VariableLabelSet(iron.FieldVariableTypes.U2,"stress")
    
    dependentField.CreateFinish()
    
    # Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
    iron.Field.ParametersToFieldParametersComponentCopy(
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
    iron.Field.ParametersToFieldParametersComponentCopy(
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
    iron.Field.ParametersToFieldParametersComponentCopy(
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
    iron.Field.ComponentValuesInitialiseDP(
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,-8.0)
    
    # Create a deformed geometry field, as cmgui doesn't like displaying
    # deformed fibres from the dependent field because it isn't a geometric field.
    deformedField = iron.Field()
    deformedField.CreateStart(deformedFieldUserNumber, region)
    deformedField.MeshDecompositionSet(decomposition)
    deformedField.TypeSet(iron.FieldTypes.GEOMETRIC)
    deformedField.VariableLabelSet(iron.FieldVariableTypes.U, "DeformedGeometry")
    for component in [1, 2, 3]:
        deformedField.ComponentMeshComponentSet(
                iron.FieldVariableTypes.U, component, 1)
    if InterpolationType == 4:
        deformedField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    deformedField.CreateFinish()
    
    # Create the equations_set
    equationsSetField = iron.Field()
    equationsSet = iron.EquationsSet()
    equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
        iron.EquationsSetTypes.FINITE_ELASTICITY,
        constitutiveRelation]
    equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,
        equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
    equationsSet.CreateFinish()
    
    # Create the CellML environment
    CellML = iron.CellML()
    CellML.CreateStart(CellMLUserNumber, region)
    
    # here we get the current path of this script so that we can specify the CellML model document relative to
    # this script, rather than the execution folder.
    from os.path import dirname, join
    script_path = dirname(__file__)
    cellmlFile = join(script_path, "guccione.cellml")
    # guccioneCellMLParameters = [0.88, 0.0, 18.5, 3.58, 3.26] # default values in CellML model
    guccioneCellMLParameters = [1.0, 0.0, 5.0, 10.0, 5.0]
    guccioneCellMLParameters = [materialParameters[0], 0.0, materialParameters[1], materialParameters[2], materialParameters[3]]
    # the names of the variables in the CellML model for the parameters in the same order as the values above
    guccioneCellMLParameterIds = ["interface/c1", "interface/c2", "interface/c3", "interface/c4", "interface/c5"]
    # Import a Guccione material law from a file
    GuccioneModel = CellML.ModelImport(cellmlFile)
    # Now we have imported the model we are able to specify which variables from the model we want to set from openCMISS
    CellML.VariableSetAsKnown(GuccioneModel, "equations/E11")
    CellML.VariableSetAsKnown(GuccioneModel, "equations/E12")
    CellML.VariableSetAsKnown(GuccioneModel, "equations/E13")
    CellML.VariableSetAsKnown(GuccioneModel, "equations/E22")
    CellML.VariableSetAsKnown(GuccioneModel, "equations/E23")
    CellML.VariableSetAsKnown(GuccioneModel, "equations/E33")
    for component, parameter in enumerate(guccioneCellMLParameterIds):
        CellML.VariableSetAsKnown(GuccioneModel, parameter)
    # and variables to get from the CellML 
    CellML.VariableSetAsWanted(GuccioneModel, "equations/Tdev11")
    CellML.VariableSetAsWanted(GuccioneModel, "equations/Tdev12")
    CellML.VariableSetAsWanted(GuccioneModel, "equations/Tdev13")
    CellML.VariableSetAsWanted(GuccioneModel, "equations/Tdev22")
    CellML.VariableSetAsWanted(GuccioneModel, "equations/Tdev23")
    CellML.VariableSetAsWanted(GuccioneModel, "equations/Tdev33")
    
    CellML.CreateFinish()
    
    # Start the creation of CellML <--> OpenCMISS field maps
    CellML.FieldMapsCreateStart()
    #Now we can set up the field variable component <--> CellML model variable mappings.
    #Map the strain components
    CellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,1,
                                  iron.FieldParameterSetTypes.VALUES,GuccioneModel,"equations/E11",
                                  iron.FieldParameterSetTypes.VALUES)
    CellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,2, iron.FieldParameterSetTypes.VALUES,GuccioneModel,"equations/E12", iron.FieldParameterSetTypes.VALUES)
    CellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,3, iron.FieldParameterSetTypes.VALUES,GuccioneModel,"equations/E13", iron.FieldParameterSetTypes.VALUES)
    CellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,4, iron.FieldParameterSetTypes.VALUES,GuccioneModel,"equations/E22", iron.FieldParameterSetTypes.VALUES)
    CellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,5, iron.FieldParameterSetTypes.VALUES,GuccioneModel,"equations/E23", iron.FieldParameterSetTypes.VALUES)
    CellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,6, iron.FieldParameterSetTypes.VALUES,GuccioneModel,"equations/E33", iron.FieldParameterSetTypes.VALUES)
    #DOC-END map strain components
        
    #DOC-START map stress components
    #Map the stress components
    CellML.CreateCellMLToFieldMap(GuccioneModel,"equations/Tdev11", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,1,iron.FieldParameterSetTypes.VALUES)
    CellML.CreateCellMLToFieldMap(GuccioneModel,"equations/Tdev12", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,2,iron.FieldParameterSetTypes.VALUES)
    CellML.CreateCellMLToFieldMap(GuccioneModel,"equations/Tdev13", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,3,iron.FieldParameterSetTypes.VALUES)
    CellML.CreateCellMLToFieldMap(GuccioneModel,"equations/Tdev22", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,4,iron.FieldParameterSetTypes.VALUES)
    CellML.CreateCellMLToFieldMap(GuccioneModel,"equations/Tdev23", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,5,iron.FieldParameterSetTypes.VALUES)
    CellML.CreateCellMLToFieldMap(GuccioneModel,"equations/Tdev33", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,6,iron.FieldParameterSetTypes.VALUES)
    
    #Finish the creation of CellML <--> OpenCMISS field maps
    CellML.FieldMapsCreateFinish()
    
    #Create the CellML models field
    CellMLModelsField = iron.Field()
    CellML.ModelsFieldCreateStart(CellMLModelsFieldUserNumber,CellMLModelsField)
    CellML.ModelsFieldCreateFinish()
    
    xidiv = 1.0/(NumberOfGaussXi+1)
    for elem in [1]:
    #Gauss point number counter    
        ctr = 0
    #Assign model for each quadraturePoint:
        for xi in range(0,NumberOfGaussXi):
            xi1 = (1.0+xi)*xidiv
            for xj in range(0,NumberOfGaussXi):
                xi2 = (1.0+xj)*xidiv
                for xk in range(0,NumberOfGaussXi):
                    xi3 = (1.0+xk)*xidiv
                    ctr = ctr + 1
                    CellMLModelsField.ParameterSetUpdateGaussPoint(iron.FieldVariableTypes.U,
                                                                   iron.FieldParameterSetTypes.VALUES,
                                                                   ctr,
                                                                   1,
                                                                   1,
                                                                   GuccioneModel)
    #Create the CellML parameters field --- the strain field
    CellMLParametersField = iron.Field()
    CellML.ParametersFieldCreateStart(CellMLParametersFieldUserNumber,CellMLParametersField)
    CellML.ParametersFieldCreateFinish()
    
    #  Create the CellML intermediate field --- the stress field
    CellMLIntermediateField = iron.Field()
    CellML.IntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber,CellMLIntermediateField)
    CellML.IntermediateFieldCreateFinish()
    
    for valueIndex, parameter in enumerate(guccioneCellMLParameterIds, 0):
        component = CellML.FieldComponentGet(GuccioneModel, iron.CellMLFieldTypes.PARAMETERS, parameter)
        print("Setting parameter: " + parameter + "; to value: " + str(guccioneCellMLParameters[valueIndex]) + "; field component: " + str(component))
        iron.Field.ComponentValuesInitialiseDP(CellMLParametersField, iron.FieldVariableTypes.U,
            iron.FieldParameterSetTypes.VALUES, component, guccioneCellMLParameters[valueIndex])
        
                    
    #equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
    #equationsSet.MaterialsCreateFinish()
    
    equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
    equationsSet.DependentCreateFinish()
    
    # Create equations
    equations = iron.Equations()
    equationsSet.EquationsCreateStart(equations)
    equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
    equations.outputType = iron.EquationsOutputTypes.NONE
    equationsSet.EquationsCreateFinish()
    
    def defineProblemSolver():
        # Define the problem
        problem = iron.Problem()
        problemSpecification = [iron.ProblemClasses.ELASTICITY,
                iron.ProblemTypes.FINITE_ELASTICITY,
                iron.ProblemSubtypes.FINITE_ELASTICITY_CELLML]
        problem.CreateStart(problemUserNumber, problemSpecification)
        problem.CreateFinish()
        
        # Create control loops
        problem.ControlLoopCreateStart()
        problem.ControlLoopCreateFinish()
        
        # Create problem solver
        nonLinearSolver = iron.Solver()
        linearSolver = iron.Solver()
        problem.SolversCreateStart()
        problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,nonLinearSolver)
        nonLinearSolver.outputType = iron.SolverOutputTypes.PROGRESS
        nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
        nonLinearSolver.NewtonLinearSolverGet(linearSolver)
        linearSolver.linearType = iron.LinearSolverTypes.DIRECT
        #linearSolver.libraryType = iron.SolverLibraries.LAPACK
        problem.SolversCreateFinish()
        
        #Create the problem solver CellML equations
        CellMLSolver = iron.Solver()
        problem.CellMLEquationsCreateStart()
        nonLinearSolver.NewtonCellMLSolverGet(CellMLSolver)
        CellMLEquations = iron.CellMLEquations()
        CellMLSolver.CellMLEquationsGet(CellMLEquations)
        CellMLEquations.CellMLAdd(CellML)
        problem.CellMLEquationsCreateFinish()
        
        # Create solver equations and add equations set to solver equations
        solver = iron.Solver()
        solverEquations = iron.SolverEquations()
        problem.SolverEquationsCreateStart()
        problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
        solver.SolverEquationsGet(solverEquations)
        solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
        equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
        problem.SolverEquationsCreateFinish()
        
        return [problem, solverEquations]
    
    def defineBoundaryConditions(solverEquations, increment):
        # Prescribe boundary conditions (absolute nodal parameters)
        boundaryConditions = iron.BoundaryConditions()
        solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
        
        #Set x=0 nodes to no x displacment in x. Set x=width nodes to 10% x displacement
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,1,1,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,3,1,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,5,1,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,7,1,iron.BoundaryConditionsTypes.FIXED,0.0)
        
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,2,1,iron.BoundaryConditionsTypes.FIXED,increment)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,4,1,iron.BoundaryConditionsTypes.FIXED,increment)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,6,1,iron.BoundaryConditionsTypes.FIXED,increment)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,8,1,iron.BoundaryConditionsTypes.FIXED,increment)
        
        # Set y=0 nodes to no y displacement
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,1,2,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,2,2,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,5,2,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,6,2,iron.BoundaryConditionsTypes.FIXED,0.0)
        
        # Set z=0 nodes to no y displacement
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,1,3,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,2,3,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,3,3,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,4,3,iron.BoundaryConditionsTypes.FIXED,0.0)
        
        solverEquations.BoundaryConditionsCreateFinish()
    
    # loop over load steps
    
    numberOfLoadSteps = 70
    displacementIncrement = 0.01 # 1%
    displacementIncrementDimension = displacementIncrement*width # length units
    resultRecord = {}
    resultRecord["strain"] = [0.0]
    resultRecord["stress"] = [0.0] 
    for counter in range(1, numberOfLoadSteps+1):
        # define the problem, solver, control loops, etc.
        [problem, solverEquations] = defineProblemSolver()
        # define the boundary conditions
        defineBoundaryConditions(solverEquations, displacementIncrementDimension)
        # execute the experiment
        problem.Solve()
        # clean up
        problem.Finalise()
        solverEquations.Finalise()
                
        # export the results
        filename = "results-{:03d}".format(counter)
    
        # Copy deformed geometry into deformed field
        for component in [1, 2, 3]:
            dependentField.ParametersToFieldParametersComponentCopy(
                                                                    iron.FieldVariableTypes.U,
                                                                    iron.FieldParameterSetTypes.VALUES, component,
                                                                    deformedField, iron.FieldVariableTypes.U,
                                                                    iron.FieldParameterSetTypes.VALUES, component)
    
        # Export results
        fields = iron.Fields()
        fields.CreateRegion(region)
        fields.NodesExport(filename, "FORTRAN")
        fields.ElementsExport(filename, "FORTRAN")
        fields.Finalise()
        
        versionNumber = 1
        derivativeNumber = 1
        nodeNumber = 8
        componentNumber = 1 # x-dirn
        reactionForceX = dependentField.ParameterSetGetNode(iron.FieldVariableTypes.DELUDELN, iron.FieldParameterSetTypes.VALUES,
                                                            versionNumber, derivativeNumber, nodeNumber, componentNumber)
        print("Reaction force (x) at counter: " + str(counter) + ": " + str(reactionForceX))
        resultRecord["strain"].append(counter * displacementIncrement)
        resultRecord["stress"].append(reactionForceX)
    
    coordinateSystem.Destroy()
    region.Destroy()
    basis.Destroy()
    
    return resultRecord