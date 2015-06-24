__author__ = 'zwan145'

# This script generates linear lagrange hexahedral meshes for a truncated ellipsoid model of the left ventricle.
# It has been adapted from a similar MATLAB script written by Sander Land.

# ZJW 07 May 2015

import os, sys
import numpy


def truncated_ellipsoid_LV(h):
    PI = numpy.pi

    assert(len(h)==3)
    nab = h(1)
    ncirc = h(2)
    nr = h(3)

    # Ellipse geometry
    top = 5
    outer_long_axis = 20
    outer_short_axis = 10
    wall_thickness = 3

    # Fibre angles
    fibre_angle_epi = -90
    fibre_angle_endo = 90

    # Generate nodal points
    nn_circ = ncirc
    nn_ab = nab + 1
    nn_r = nr + 1

    iz = 1 # Count transmural layer.
    coords = numpy.zeros([5, nab+1, nn_circ+1, iz])

    r = numpy.r_[(outer_long_axis-wall_thickness):outer_long_axis:(wall_thickness/nr)]

    for long_r in r:
        short_r = outer_short_axis - (outer_long_axis - long_r)

        uu = -numpy.arccos(top/long_r)  # Angle at base.
        u = numpy.r_[-PI:uu:((PI+uu)/(nab))]
        v = numpy.r_[-PI:PI:(2*PI/ncirc)]

        coords[1,:,:,iz] = short_r*numpy.outer(numpy.sin(u),numpy.cos(v))
        coords[2,:,:,iz] = short_r*numpy.outer(numpy.sin(u),numpy.sin(v))
        coords[3,:,:,iz] = long_r*numpy.outer(numpy.cos(u),numpy.ones([len(v)]))
        coords[4,:,:,iz] = numpy.tile(u, (1, len(v)))
        coords[5,:,:,iz] = numpy.tile(v, (len(u), 1))

        iz = iz + 1 # Next layer out from endocardium.


    # Generate mesh - number nodes naively with duplicates first.
    dof_per_elem = 8
    el_nodes = numpy.zeros([dof_per_elem, 5])
    naive_elem = numpy.zeros([nab*ncirc*nr, dof_per_elem])
    naive_nodes = numpy.zeros([dof_per_elem*nab*ncirc*nr, 5])
    naive_epi1endo0 = numpy.zeros([dof_per_elem*nab*ncirc*nr,1])

    ei = 0
    for el_r in range(0, nr):
        for el_ab in range(0, nab):
            for el_circ in range(0, ncirc):
                ix_ab = (el_ab-1) + numpy.r_[0:2]
                ix_c = (el_circ-1) + numpy.r_[0:2]
                ix_r = (el_r-1) + numpy.r_[0:2]

                for xi3 in numpy.r_[0:2]:
                    for xi2 in numpy.r_[0:2]:
                        for xi1 in numpy.r_[0:2]:
                            el_nodes[4*xi3+2*xi2+xi1, :] = coords[:,ix_ab[xi2], ix_c[xi1], ix_r[xi3]]

                ei = ei + 1
                nn = (ei-1)*8 + numpy.r_[0:8]
                naive_elem[ei,:] = nn
                naive_nodes[nn,:] = el_nodes
                temp = numpy.r_[0:8]/4.0
                for i in temp:
                    naive_epi1endo0[nn] = (el_r-1)/nr + numpy.floor(i)/nr

    # Get rid of duplicate nodes
    nodes = naive_nodes[:, 0:3]
    temp = len(nodes)

    [nodes, reversenodemap, nodemap] = numpy.unique(nodes, True, True)
    elem = nodemap(naive_elem)
    print 'removed duplicates '+str(temp)+' -> '+str(len(nodes)) + '\n'

    UV = naive_nodes[reversenodemap, 4:5]

    

