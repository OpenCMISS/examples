'''
Created on 24/03/2014

@author: rjag008
'''

from opencmiss.zinc.context import Context
# from opencmiss.zinc.field import Field
from opencmiss.zinc.element import Element, Elementbasis
from copy import deepcopy
import math

openCMISSOutput = True

context = Context("Quad cylinder")
region = context.getDefaultRegion()

#Number of nodes in each direction (num elements = num nodes -1)
xitref = 11 #Along theta
xihref = 11 #Along height
xirref = 5  #Along radius


#Setup for a unit cube
tstep = 1.0 / (xitref-1)
hstep = 5.0 / (xihref-1)
rstep = 1.0 / (xirref-1) 


nodes = []
nodeAngle = dict() #Keep track of the node angle
coordinates = []
nodectr = 0

for r in range(0, xirref):
    #x = r*rstep
    radius = (r+1)*rstep*0.5 + 0.5
    for t in range(0, xitref):
        #y = t*tstep
        angle = t*tstep*2*math.pi
        x = radius*math.cos(angle)
        y = radius*math.sin(angle)
        for h in range(0, xihref):
            nodes.append(nodectr + 1)
            z = h*hstep
            nodeAngle[nodectr+1] = angle
            coordinates.append([x,y,z])
            nodectr = nodectr + 1

sheet = xihref * xitref

# firstelement = [1, sheet + 1, xihref + 1, sheet + xihref + 1, 2, sheet + 2, xihref + 2, sheet + xihref + 2]
# Quad has 27 nodes
    
firstelement = [1, sheet + 1, 2 * sheet + 1,
                    xihref + 1, sheet + xihref + 1, 2 * sheet + xihref + 1,
                    2 * xihref + 1, sheet + 2 * xihref + 1, 2 * sheet + 2 * xihref + 1,
                    2, sheet + 2, 2 * sheet + 2,
                    xihref + 2, sheet + xihref + 2, 2 * sheet + xihref + 2,
                    2 * xihref + 2, sheet + 2 * xihref + 2, 2 * sheet + 2 * xihref + 2,
                    3, sheet + 3, 2 * sheet + 3,
                    xihref + 3, sheet + xihref + 3, 2 * sheet + xihref + 3,
                    2 * xihref + 3, sheet + 2 * xihref + 3, 2 * sheet + 2 * xihref + 3]

elements = []
selectedNodes = dict()
helem = (xihref - 1) / 2
telem = (xitref - 1) / 2
relem = (xirref - 1) / 2
for z in range(0, helem):
    for y in range(0, telem):
        for x in range(0, relem):
            element = deepcopy(firstelement)
            offset = 2 * (z + y * xihref + x * sheet)
            for i in range(0, 27):
                element[i] = element[i] + offset
            if y == telem -1:
                for i in range(6,9):
                    element[i] = element[i] - (xitref-1)*xihref
                for i in range(15,18):
                    element[i] = element[i] - (xitref-1)*xihref
                for i in range(24,27):
                    element[i] = element[i] - (xitref-1)*xihref
            elements.append(element)
            for nd in element:
                selectedNodes[nd] = True                    
          
#Ensure that elements with more than 50% of the material, all of  its material assigned the same material            
            

field_module = region.getFieldmodule()
        
field_module.beginChange()
finite_element_field = field_module.createFieldFiniteElement(3)
# Set the name of the field, we give it label to help us understand it's purpose
finite_element_field.setName('coordinates')
finite_element_field.setTypeCoordinate(True)

material_field = field_module.createFieldFiniteElement(1)
# Set the name of the field, we give it label to help us understand it's purpose
material_field.setName('material')


# Find a special node set named 'cmiss_nodes'
nodeset = field_module.findNodesetByName('nodes')
node_template = nodeset.createNodetemplate()
# Set the finite element coordinate field for the nodes to use
node_template.defineField(finite_element_field)
node_template.defineField(material_field)

field_cache = field_module.createFieldcache()
# Create nodes
nodeHandle = []
numNodes = len(nodes)
for index in range(numNodes):
# Create a node if it is accesed
    nid = index + 1
    if nid in selectedNodes:
        node = nodeset.createNode(nid, node_template)
        nodeHandle.append(node)
# Set the node coordinates, first set the field cache to use the current node
        field_cache.setNode(node)
# Pass in floats as an array
        finite_element_field.assignReal(field_cache, coordinates[index])
        mtype = 1
        if coordinates[index][2] > 5.0/3 and coordinates[index][2] < 10.0/3:
            mtype = 2 
        elif coordinates[index][2] > 10.0/3:
            mtype = 3    
        material_field.assignReal(field_cache, [mtype])
    else:
# Append dummy to maintain indexing        
        nodeHandle.append("")

mesh = field_module.findMeshByDimension(3)
element_template = mesh.createElementtemplate()
element_template.setElementShapeType(Element.SHAPE_TYPE_CUBE)
element_template.setNumberOfNodes(27)
# Specify the dimension and the interpolation function for the element basis function. 
quadratic_basis = field_module.createElementbasis(3, Elementbasis.FUNCTION_TYPE_QUADRATIC_LAGRANGE)
# The indexes of the nodes in the node template we want to use
node_indexes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]
# Define a nodally interpolated element field or field component in the
# element_template. Only Lagrange, simplex and constant basis function types
# may be used with this function, i.e. where only a simple node value is
# mapped. Shape must be set before calling this function.  The -1 for the component number
# defines all components with identical basis and nodal mappings.
element_template.defineFieldSimpleNodal(finite_element_field, -1, quadratic_basis, node_indexes)
element_template.defineFieldSimpleNodal(material_field, -1, quadratic_basis, node_indexes)

#Track the active nodes
activeNodes = dict()
for index in range(len(elements)):
    elem = elements[index]
    try:
        for i, el in enumerate(elem):
            inode = int(el)
            node = nodeHandle[inode - 1]
            activeNodes[inode] = node
            element_template.setNode(i + 1, node)
    except:
        print elem
    mesh.defineElement(-1, element_template)

offset = xirref*xihref*xitref + 100;
#Reset node ids to be continous
#First offset it to avoid collisons
for node in activeNodes:
    nid = offset + activeNodes[node].getIdentifier()
    activeNodes[node].setIdentifier(offset+nid)
    
for i,node in enumerate(activeNodes):
    activeNodes[node].setIdentifier(i+1)    

     
if openCMISSOutput==False:     
    field_module.defineAllFaces()
     
field_module.endChange()

# Use the following commands in cmgui to reorder numbering
# gfx define field sort_field composite xi_texture.2 xi_texture.3 xi_texture.1
# Sort the nodes and elements based on this field
# gfx change_identifier sort_by sort_field node_offset 1 element_offset 1 line_offset 1 face_offset 1;
  
      
# Output mesh
sir = region.createStreaminformationRegion()
fr = sir.createStreamresourceFile("hetrogenouscylinder.exregion")
region.write(sir)       

