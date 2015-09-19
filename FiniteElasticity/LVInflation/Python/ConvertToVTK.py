__author__ = 'zwan145'

# This script takes user input for mesh resolutions and converts the corresponding
# solution exnode and exelem files to vtk format for display in paraview. 

# Author: ZJW

import os, sys
import exfile


#===============================================================================#
def Convert2VTK(f_exnode, f_exelem):

    inputNodes = exfile.Exnode(f_exnode)
    inputElems = exfile.Exelem(f_exelem)

    def ExtractNodeCoords(nodes, field_name):
        coords = []
        for node_num in range(1, nodes.num_nodes+1):
            temp = []
            for component in [1,2,3]:
                component_name = ["x","y","z"][component-1]
                value = nodes.node_value(field_name, component_name, node_num, 1)
                temp.append(value)
            coords.append(temp)
        coords = tuple(coords)
        return coords

    def ExtractElems(elems):
        coords = []
        for e in range(0, elems.num_elements):
            temp = elems.elements[e].nodes
            coords.append(temp)
        coords = tuple(coords)
        return coords

    nodes = ExtractNodeCoords(inputNodes, "DeformedGeometry")
    elems = ExtractElems(inputElems)

    filename = "output_"+str(dimensions[0])+"-"+str(dimensions[1])+"-"+str(dimensions[2])+"'.vtk"
    fid = open(filename,'w')
    fid.write('# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS '+str(inputNodes.num_nodes)+' float\n')

    for i in range(0, inputNodes.num_nodes):
        fid.write(str(inputNodes.node_value("DeformedGeometry", "x", i+1, 1))+' '+
                  str(inputNodes.node_value("DeformedGeometry", "y", i+1, 1))+' '+
                  str(inputNodes.node_value("DeformedGeometry", "z", i+1, 1))+'\n')

    fid.write('CELLS '+str(inputElems.num_elements)+' '+str(9*inputElems.num_elements)+'\n')

    for e in range(0, inputElems.num_elements):
        fid.write('8\t'+str(elems[e][0]-1)+' '+str(elems[e][1]-1)+' '+str(elems[e][3]-1)+' '+str(elems[e][2]-1)+' '+str(elems[e][4]-1)+' '+str(elems[e][5]-1)+' '+str(elems[e][7]-1)+' '+str(elems[e][6]-1)+'\n')

    fid.write('CELL_TYPES '+str(inputElems.num_elements)+'\n')

    for e in range(0, inputElems.num_elements):
        fid.write('12\n')

    fid.close()


#===============================================================================#
def Convert2VTK_undef(f_exnode, f_exelem):

    inputNodes = exfile.Exnode(f_exnode)
    inputElems = exfile.Exelem(f_exelem)

    def ExtractNodeCoords(nodes, field_name):
        coords = []
        for node_num in range(1, nodes.num_nodes+1):
            temp = []
            for component in [1,2,3]:
                component_name = ["x","y","z"][component-1]
                value = nodes.node_value(field_name, component_name, node_num, 1)
                temp.append(value)
            coords.append(temp)
        coords = tuple(coords)
        return coords

    def ExtractElems(elems):
        coords = []
        for e in range(0, elems.num_elements):
            temp = elems.elements[e].nodes
            coords.append(temp)
        coords = tuple(coords)
        return coords

    nodes = ExtractNodeCoords(inputNodes, "Geometry")
    elems = ExtractElems(inputElems)

    filename = 'undeformed.vtk'
    fid = open(filename,'w')
    fid.write('# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS '+str(inputNodes.num_nodes)+' float\n')

    for i in range(0, inputNodes.num_nodes):
        fid.write(str(inputNodes.node_value("Geometry", "x", i+1, 1))+' '+
                  str(inputNodes.node_value("Geometry", "y", i+1, 1))+' '+
                  str(inputNodes.node_value("Geometry", "z", i+1, 1))+'\n')

    fid.write('CELLS '+str(inputElems.num_elements)+' '+str(9*inputElems.num_elements)+'\n')

    for e in range(0, inputElems.num_elements):
        fid.write('8\t'+str(elems[e][0]-1)+' '+str(elems[e][1]-1)+' '+str(elems[e][3]-1)+' '+str(elems[e][2]-1)+' '+str(elems[e][4]-1)+' '+str(elems[e][5]-1)+' '+str(elems[e][7]-1)+' '+str(elems[e][6]-1)+'\n')

    fid.write('CELL_TYPES '+str(inputElems.num_elements)+'\n')

    for e in range(0, inputElems.num_elements):
        fid.write('12\n')

    fid.close()

#===============================================================================#
# User input to select mesh resolution to convert. 
elems = [0,0,0]
elems[0] = int(sys.argv[1])
elems[1] = int(sys.argv[2])
elems[2] = int(sys.argv[3])

f_exnode = "LVInflation_trilinear_"+str(elems[0])+"-"+str(elems[1])+"-"+str(elems[2])+".part0.exnode"
f_exelem = "ellipse_benchmark_lin_"+str(elems[0])+"-"+str(elems[1])+"-"+str(elems[2])+".exelem"
print "Converting solution with dimensions "+str(elems)+" to VTK format\n"
Convert2VTK(f_exnode, f_exelem)
Convert2VTK_undef(f_exnode, f_exelem)

