__author__ = 'zwan145'

import numpy

def RewritePythonScript(filenamer, filenamew, elems):
    try:
        fidr = open(filenamer, 'r')
        fidw = open(filenamew, 'w')
    except IOError:
        print 'ERROR: Problem with reading or writing python scripts'
        return

    n = 0
    data = fidr.readline()
    n = n + 1
    while (n<276):
        if (n==62):
            temp = 'elems = ['+str(elems[0])+', '+str(elems[1])+','+str(elems[2])+']\n'
        else:
            temp = data
        fidw.write(str(temp))
        data = fidr.readline()
        n = n + 1

def BCFaces(numElem):
    backFaceNodes = [0,0,0,0]

    backFaceNodes[0] = 1
    backFaceNodes[1] = (numElem+1)+1
    backFaceNodes[2] = (numElem+1)*2 + 1
    backFaceNodes[3] = (numElem+1)*3 + 1

    bottomFaceNodes = numpy.linspace(1,(numElem+1)*2,(numElem+1)*2)

    return backFaceNodes, bottomFaceNodes

def ExtractInterpt(numElem):
    # Read displacement data file
    filename = 'problem1_tricubic_'+str(numElem)+'elem_displacement.exdata'
    fid = open(filename, 'r')
    for i in range(1,11):
        junk = fid.readline()
    ptsX = []
    ptsY = []
    ptsZ = []
    temp = fid.readline()
    while temp!='':
        junk = fid.readline()
        temp = fid.readline()
        ptsX.append(float(temp.split()[0]))
        temp = fid.readline()
        ptsY.append(float(temp.split()[0]))
        temp = fid.readline()
        ptsZ.append(float(temp.split()[0]))
        temp = fid.readline()

    ptsX = numpy.array(ptsX)
    ptsY = numpy.array(ptsY)
    ptsZ = numpy.array(ptsZ)
    disp_pts = [ptsX, ptsY, ptsZ]
    disp_pts = numpy.array(disp_pts)
    return disp_pts

def ExtractTff(numElem):
    Tff_pts = []
    filename = 'problem1_tricubic_'+str(numElem)+'elem_stress2PK.exdata'
    fid = open(filename, 'r')
    for i in range(1,17):
        junk = fid.readline()

    temp = fid.readline()
    while temp != '':
        junk = fid.readline()
        temp = fid.readline()
        Tff_pts.append(float(temp.split()[0]))
        for i in range(1,6):
            junk = fid.readline()
        temp = fid.readline()
    Tff_pts = numpy.array(Tff_pts)
    return Tff_pts


