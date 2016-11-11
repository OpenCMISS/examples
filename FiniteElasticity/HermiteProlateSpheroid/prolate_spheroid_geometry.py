""" Code for generating a cutoff prolate spheroid mesh to
    represent the left ventricle of the heart.
    Uses a rectangular Cartesian coordinate system and supports
    linear Lagrange, quadratic Lagrange and cubic Hermite interpolation.
"""

from collections import defaultdict
import itertools
import numpy as np
from numpy import pi, sin, cos, sinh, cosh

from opencmiss.iron import iron


class ProlateSpheroid(object):
    """ A prolate spheroid mesh where xi_1 corresponds to the theta direction,
        xi_2 corresponds to the mu direction, and xi_3 corresponds to the lambda direction.
        The nodes at the apex have DOFs constrained to collapse the apex faces.
    """
    def __init__(self, focus, internalLambda, externalLambda, cutoffAngle,
            numElements,
            endocardiumFibreAngle, epicardiumFibreAngle, sheetAngle,
            interpolations=None):
        self.focus = focus
        self.internalLambda = internalLambda
        self.externalLambda = externalLambda
        self.cutoffAngle = cutoffAngle

        if interpolations is not None:
            self.interpolations = interpolations
        else:
            self.interpolations = ['quadratic', 'linear']

        # Circumferential, longitudinal, transmural
        self.numElements = numElements

        self.endocardiumFibreAngle = endocardiumFibreAngle
        self.epicardiumFibreAngle = epicardiumFibreAngle
        self.sheetAngle = sheetAngle

        self._calculateNodes()
        self._calculateElements()

    def elements(self, component):
        """ Returns a list of elements, where each element
            is a tuple containing a list of node numbers
            and a boolean indicating whether the element is collapsed.
            Component is either "linear" or "quadratic"
        """
        return self.componentElements[component]

    def nodes(self):
        """ Node positions
        """
        return self.nodePositions

    def nodeValues(self):
        """ Node positions
        """
        return self._nodeValues

    def componentNodes(self, component):
        """ Node numbers in the given component
        """
        allNodes = np.concatenate(
                [element[0] for element in self.elements(component)])
        return set(int(node) for node in allNodes)

    def meshComponent(self, interpolation):
        components = dict(
                (interp, i)
                for i, interp in
                enumerate(self.interpolations, 1))
        return components[interpolation]

    def totalNumElements(self, component):
        return len(self.componentElements[component])

    def numNodes(self):
        return len(self.nodePositions)

    def nodeAtPosition(self, position):
        return self.positionNode[position]

    def nodeGroup(self, group):
        """ Returns nodes in the named surface/group
        """
        if 'quadratic' in self.interpolations:
            elFactor = 2
        else:
            elFactor = 1
        if group == 'base':
            return [n for p, n in self.positionNode.iteritems()
                    if p[1] == self.numElements[1] * elFactor]
        elif group in ('internal', 'endocardium'):
            return [n for p, n in self.positionNode.iteritems()
                    if p[2] == 0]
        elif group in ('external', 'epicardium'):
            return [n for p, n in self.positionNode.iteritems()
                    if p[2] == self.numElements[2] * elFactor]
        else:
            raise ValueError("Invalid node group: %s" % group)

    def generateMesh(self, region):
        bases = self.setupBases(region)

        # Start the creation of a mesh in the region, setting
        # the number of mesh elements and number of mesh components
        # There are two mesh components, one quadratic and one linear
        mesh = iron.Mesh()
        mesh.CreateStart(1, region, 3)
        mesh.NumberOfComponentsSet(len(self.interpolations))
        mesh.NumberOfElementsSet(self.totalNumElements(self.interpolations[0]))

        # Create nodes in the region, setting the total number required
        # for the prolate spheroid geometry
        meshNodes = iron.Nodes()
        meshNodes.CreateStart(region, self.numNodes())
        meshNodes.CreateFinish()

        # Create mesh component elements for each of the
        # linear and quadratic mesh components
        for meshComponent, interpolation in enumerate(self.interpolations, 1):
            meshElements = iron.MeshElements()
            # Set the default basis that has no collapsed nodes
            basis = bases[(interpolation, False)]
            meshElements.CreateStart(mesh, meshComponent, basis)
            for elementNum, element in enumerate(self.elements(interpolation), 1):
                elementNodes, collapsed = element
                if collapsed:
                    basis = bases[(interpolation, collapsed)]
                    meshElements.BasisSet(elementNum, basis)
                meshElements.NodesSet(elementNum, elementNodes)
            meshElements.CreateFinish()

        mesh.CreateFinish()
        return mesh

    def setGeometry(self, geometricField):
        decomposition = iron.Decomposition()
        geometricField.MeshDecompositionGet(decomposition)
        geometricInterpolation = self.interpolations[0]
        compNodeNumber = iron.ComputationalNodeNumberGet()
        geometricMeshComponent = 1
        # Set the geometric field parameters from the prolate spheroid geometry
        for nodeNum, values in enumerate(self.nodeValues(), 1):
            if decomposition.NodeDomainGet(nodeNum, geometricMeshComponent) == compNodeNumber:
                versionNumber = 1
                for component in range(3):
                    componentValues = values[component]
                    if geometricInterpolation != 'cubic_hermite':
                        componentValues = componentValues[:1]
                    for derivativeNumber, value in componentValues:
                        # Note: Python uses 0-based indexing, OpenCMISS uses 1-based
                        geometricField.ParameterSetUpdateNode(
                                iron.FieldVariableTypes.U,
                                iron.FieldParameterSetTypes.VALUES,
                                versionNumber, derivativeNumber, nodeNum, component + 1,
                                value)
        # After modifying the geometric field, ParameterSetUpdateStart/Finish
        # must be called so that the field data at ghost nodes is synchronised
        # across all processes. The UpdateFinish routine could be called later
        # on just before the geometric values are actually needed, to avoid
        # blocking.
        geometricField.ParameterSetUpdateStart(
                            iron.FieldVariableTypes.U,
                            iron.FieldParameterSetTypes.VALUES)
        geometricField.ParameterSetUpdateFinish(
                            iron.FieldVariableTypes.U,
                            iron.FieldParameterSetTypes.VALUES)

    def setFibres(self, fibreField):
        endocardiumFibreAngle = self.endocardiumFibreAngle
        epicardiumFibreAngle = self.epicardiumFibreAngle
        sheetAngle = self.sheetAngle
        geometricMeshComponent = 1
        compNodeNumber = iron.ComputationalNodeNumberGet()

        hasQuadratic = self.interpolations[0] == 'quadratic'
        if hasQuadratic:
            elFactor = 2
        else:
            elFactor = 1

        # fibre, imbrication, sheet
        epiAngles = np.array([epicardiumFibreAngle, 0.0, sheetAngle])
        endoAngles = np.array([endocardiumFibreAngle, 0.0, sheetAngle])
        lmbdaDerivs = (epiAngles - endoAngles) / self.numElements[2]

        decomposition = iron.Decomposition()
        fibreField.MeshDecompositionGet(decomposition)
        # Set fibre angles at nodes
        def setAngles(pos, angles):
            nodeNumber = self.nodeAtPosition(pos)
            if decomposition.NodeDomainGet(nodeNumber, geometricMeshComponent) == compNodeNumber:
                version = 1
                for component, angle in enumerate(angles, 1):
                    derivative = 1
                    fibreField.ParameterSetUpdateNode(
                        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                        version, derivative, nodeNumber,
                        component, angle)

        # Only loop over linear component nodes
        for k in range(0, self.numElements[2] * elFactor + 1, elFactor):
            transmuralPos = float(k) / float(self.numElements[2] * elFactor)
            fibreAngle = endocardiumFibreAngle + transmuralPos * (
                    epicardiumFibreAngle - endocardiumFibreAngle)
            for j in range(0, self.numElements[1] * elFactor + 1, elFactor):
                for i in range(0, self.numElements[0] * elFactor, elFactor):
                    if j == 0:
                        setAngles((i, j, k), [0.0, 0.0, 0.0])
                    else:
                        setAngles((i, j, k), [fibreAngle, 0.0, sheetAngle])

        fibreField.ParameterSetUpdateStart(
                            iron.FieldVariableTypes.U,
                            iron.FieldParameterSetTypes.VALUES)
        fibreField.ParameterSetUpdateFinish(
                            iron.FieldVariableTypes.U,
                            iron.FieldParameterSetTypes.VALUES)

    def indicesToPosition(self, indices):
        nodeNum = self.nodeAtPosition(indices)
        return np.array(self.nodes()[nodeNum - 1])

    def setupBases(self, region):
        """ Set up required bases in an OpenCMISS region
        """
        # Number of Gauss points used for integration, depends on
        # maximum interpolation degree and has to be the same
        # for all bases.
        interpolationGaussRequired = {
                'cubic_hermite': 4,
                'quadratic': 3,
                'linear': 2,
                }
        numberOfGaussXi = max(
                interpolationGaussRequired[interp]
                for interp in self.interpolations)

        # Define bases
        self.bases = {}
        def makeBasis(userNumber, interpolationName, collapsed):
            """ Helper function for creating a basis of the
                given interpolation type and collapsed state
            """
            cmissInterpolations = {
                    'cubic_hermite': iron.BasisInterpolationSpecifications.CUBIC_HERMITE,
                    'quadratic': iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE,
                    'linear': iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE,
                    }
            interpolation = cmissInterpolations[interpolationName]
            basis = iron.Basis()
            basis.CreateStart(basisUserNumber)
            basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
            basis.NumberOfXiSet(3)
            basis.InterpolationXiSet([interpolation] * 3)
            basis.QuadratureNumberOfGaussXiSet([numberOfGaussXi] * 3)
            # Have to enable face basis evaluation so we can
            # use pressure boundary conditions
            basis.QuadratureLocalFaceGaussEvaluateSet(True)
            if collapsed:
                basis.CollapsedXiSet([
                    iron.BasisXiCollapse.XI_COLLAPSED,
                    iron.BasisXiCollapse.COLLAPSED_AT_XI0,
                    iron.BasisXiCollapse.NOT_COLLAPSED,
                    ])
            basis.CreateFinish()
            # Save basis into a dictionary for easy
            # retrieval later
            self.bases[interpolationName, collapsed] = basis

        for i, interpolation in enumerate(self.interpolations):
            # Define standard basis
            basisUserNumber = i + 1
            makeBasis(basisUserNumber, interpolation, False)

        return self.bases

    def finaliseBases(self):
        for basis in self.bases.values():
            basis.Finalise()
        self.bases = {}

    def _calculateElements(self):
        includeQuadratic = 'quadratic' in self.interpolations
        def calcNodeAxes(i, j, k):
            """ Returns node indices to use in each xi direction
            """
            if includeQuadratic:
                # Account for loop around in theta:
                finalIPos = (0 if i == (self.numElements[0] - 1)
                        else i * 2 + 2)
                linearAxes = [
                            [k * 2, k * 2 + 2],
                            [j * 2, j * 2 + 2],
                            [i * 2, finalIPos]]
                nodeAxes = {
                        'cubic_hermite': linearAxes,
                        'linear': linearAxes,
                        'quadratic': [
                            [k * 2, k * 2 + 1, k * 2 + 2],
                            [j * 2, j * 2 + 1, j * 2 + 2],
                            [i * 2, i * 2 + 1, finalIPos]]}
            else:
                # Account for loop around in theta:
                finalIPos = (0 if i == (self.numElements[0] - 1)
                        else i + 1)
                linearAxes = [
                            [k, k + 1],
                            [j, j + 1],
                            [i, finalIPos]]
                nodeAxes = {
                        'linear': linearAxes,
                        'cubic_hermite': linearAxes}
            return nodeAxes

        self.componentElements = defaultdict(list)
        # Loop from endocardium to epicardium
        for k in xrange(self.numElements[2]):
            # Loop from apex to base
            for j in xrange(self.numElements[1]):
                # Loop circumferentially
                for i in xrange(self.numElements[0]):
                    nodeAxes = calcNodeAxes(i, j, k)
                    for component in self.interpolations:
                        positions = nodeAxes[component]
                        collapsed = False
                        nodes = [self.positionNode[ni, nj, nk]
                                for nk, nj, ni in
                                    itertools.product(*positions)]
                        self.componentElements[component].append(
                                (nodes, collapsed))

    def _calculateNodes(self):
        """ Calculates node positions as well as
            a map from position indices to node numbers
        """
        includeQuadratic = 'quadratic' in self.interpolations
        if includeQuadratic:
            elFactor = 2
        else:
            elFactor = 1
        lambdaPositions = np.linspace(
                self.internalLambda, self.externalLambda,
                self.numElements[2] * elFactor + 1)
        muPositions = np.linspace(
                0.0, self.cutoffAngle,
                self.numElements[1] * elFactor + 1)
        thetaPositions = np.linspace(
                0.0, 2.0 * pi,
                self.numElements[0] * elFactor + 1)
        # Account for wrap around, last nodes same as first:
        thetaPositions = thetaPositions[:-1]

        self.nodePositions = []
        self._nodeValues = []
        self.positionNode = {}
        nodeNumber = itertools.count(1)

        # Loop from endocardium to epicardium
        for k in xrange(self.numElements[2] * elFactor + 1):
            # Apex nodes at mu=0, have multiple apex nodes at same
            # position but different theta, get constrained to be
            # mapped together later:
            # Loop from apex to base
            for j in xrange(self.numElements[1] * elFactor + 1):
                # Loop circumferentially
                for i in xrange(self.numElements[0] * elFactor):
                    self.nodePositions.append(xyz(self.focus,
                        lambdaPositions[k], muPositions[j], thetaPositions[i]))
                    self.positionNode[(i, j, k)] = next(nodeNumber)
                    # list of x, y, z derivatives
                    # each set is no_deriv, s1, s2, s1_s2, s3, s1_s3, s2_s3, s1_s2_s3
                    derivatives = self.calculateDerivatives(
                        lambdaPositions[k], muPositions[j], thetaPositions[i])
                    if j == 0:
                        # At apex, so xi_1 derivs are zero
                        for component in range(3):
                            derivatives[component][1] = (2, 0.0) # S1
                            derivatives[component][3] = (4, 0.0) # S1_S2
                            derivatives[component][5] = (6, 0.0) # S1_S3
                            derivatives[component][7] = (8, 0.0) # S1_S2_S3
                    self._nodeValues.append(derivatives)

    def calculateDerivatives(self, lmbda, mu, theta):
        """ Calculate x, y, z derivatives wrt xi directions,
            for use with unit scaling.
        """
        return [[(derivNumber, self.xXiDerivative(lmbda, mu, theta, xCoord, derivNumber))
                    for derivNumber in range(1, 9)]
                for xCoord in range(3)]

    def xXiDerivative(self, lmbda, mu, theta, xCoordNum, derivativeNumber):
        """ Calculate dx_i/dxi at prolate spheroidal coordinat lmbda, mu, theta
            where i = xCoordNum and is referenced from 0.
            derivativeNumber gives xi coord to calc deriv wrt to, can be
            a double or triple derivative, is one of iron.GlobalDerivativeConstants
        """
        focus = self.focus
        pCoords = (lmbda, mu, theta)
        if derivativeNumber == iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV:
            return xyz(focus, lmbda, mu, theta)[xCoordNum]
        elif derivativeNumber == iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1:
            return (self.xProlateDerivative(xCoordNum, (0, ), pCoords) *
                    self.prolateXiDerivative(0, 0))
        elif derivativeNumber == iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2:
            return (self.xProlateDerivative(xCoordNum, (1, ), pCoords) *
                    self.prolateXiDerivative(1, 1))
        elif derivativeNumber == iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2:
            return (self.xProlateDerivative(xCoordNum, (0, 1), pCoords) *
                    self.prolateXiDerivative(0, 0) *
                    self.prolateXiDerivative(1, 1))
        elif derivativeNumber == iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3:
            return (self.xProlateDerivative(xCoordNum, (2, ), pCoords) *
                    self.prolateXiDerivative(2, 2))
        elif derivativeNumber == iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3:
            return (self.xProlateDerivative(xCoordNum, (0, 2), pCoords) *
                    self.prolateXiDerivative(0, 0) *
                    self.prolateXiDerivative(2, 2))
        elif derivativeNumber == iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3:
            return (self.xProlateDerivative(xCoordNum, (1, 2), pCoords) *
                    self.prolateXiDerivative(1, 1) *
                    self.prolateXiDerivative(2, 2))
        elif derivativeNumber == iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3:
            return (self.xProlateDerivative(xCoordNum, (0, 1, 2), pCoords) *
                    self.prolateXiDerivative(0, 0) *
                    self.prolateXiDerivative(1, 1) *
                    self.prolateXiDerivative(2, 2))
        else:
            raise ValueError("Invalid derivative number: {0:d}".format(derivativeNumber))

    def xProlateDerivative(self, xCoordNum, prolateCoordNums, prolateCoords):
        """ Returns d x_i / d ._j where i is xCoordNum (indexed from zero)
            and "." is one of theta, mu or lambda, corresponding to
            prolateCoordNum from zero
            prolateCoordNums is a tuple of derivatives, to allow second and third
            partial derivatives.
            Prolate coords is current position in [theta, mu, lambda]
        """
        lmbda, mu, theta = prolateCoords
        d = self.focus
        if xCoordNum == 0:
            # x = d sinh(lmbda) sin(mu) cos(theta)
            if prolateCoordNums == (0, ):
                # d x / d theta
                return - d * sinh(lmbda) * sin(mu) * sin(theta)
            elif prolateCoordNums == (1, ):
                # d x / d mu
                return d * sinh(lmbda) * cos(mu) * cos(theta)
            elif prolateCoordNums == (2, ):
                # d x / d lambda
                return d * cosh(lmbda) * sin(mu) * cos(theta)
            elif prolateCoordNums == (0, 1):
                # d x^2 / d theta d mu
                return - d * sinh(lmbda) * cos(mu) * sin(theta)
            elif prolateCoordNums == (0, 2):
                # d x^2 / d theta d lmbda
                return - d * cosh(lmbda) * sin(mu) * sin(theta)
            elif prolateCoordNums == (1, 2):
                # d x^2 / d mu d lmbda
                return d * cosh(lmbda) * cos(mu) * cos(theta)
            elif prolateCoordNums == (0, 1, 2):
                # d x^3 / d theta d mu d lmbda
                return - d * cosh(lmbda) * cos(mu) * sin(theta)
            else:
                raise ValueError("Invalid prolate coords: {0}".format(prolateCoordNums))
        elif xCoordNum == 1:
            # y = d sinh(lmbda) sin(mu) sin(theta)
            if prolateCoordNums == (0, ):
                # d y / d theta
                return d * sinh(lmbda) * sin(mu) * cos(theta)
            elif prolateCoordNums == (1, ):
                # d y / d mu
                return d * sinh(lmbda) * cos(mu) * sin(theta)
            elif prolateCoordNums == (2, ):
                # d y / d lambda
                return d * cosh(lmbda) * sin(mu) * sin(theta)
            elif prolateCoordNums == (0, 1):
                # d y^2 / d theta d mu
                return d * sinh(lmbda) * cos(mu) * cos(theta)
            elif prolateCoordNums == (0, 2):
                # d y^2 / d theta d lmbda
                return d * cosh(lmbda) * sin(mu) * cos(theta)
            elif prolateCoordNums == (1, 2):
                # d y^2 / d mu d lmbda
                return d * cosh(lmbda) * cos(mu) * sin(theta)
            elif prolateCoordNums == (0, 1, 2):
                # d y^3 / d theta d mu d lmbda
                return d * cosh(lmbda) * cos(mu) * cos(theta)
            else:
                raise ValueError("Invalid prolate coords: {0}".format(prolateCoordNums))
        elif xCoordNum == 2:
            # z = d cosh(lmbda) cos(mu)
            if prolateCoordNums == (0, ):
                # d z / d theta
                return 0.0
            elif prolateCoordNums == (1, ):
                # d z / d mu
                return - d * cosh(lmbda) * sin(mu)
            elif prolateCoordNums == (2, ):
                # d z / d lambda
                return d * sinh(lmbda) * cos(mu)
            elif prolateCoordNums == (0, 1):
                # d z^2 / d theta d mu
                return 0.0
            elif prolateCoordNums == (0, 2):
                # d z^2 / d theta d lmbda
                return 0.0
            elif prolateCoordNums == (1, 2):
                # d z^2 / d mu d lmbda
                return - d * sinh(lmbda) * sin(mu)
            elif prolateCoordNums == (0, 1, 2):
                # d z^3 / d theta d mu d lmbda
                return 0.0
            else:
                raise ValueError("Invalid prolate coords: {0}".format(prolateCoordNums))

    def prolateXiDerivative(self, prolateCoordNum, xiCoordNum):
        """ Returns d . / d xi_i where i is xiCoordNum (indexed from zero)
            and "." is one of theta, mu or lambda, corresponding to
            prolateCoordNum from zero
        """
        prolateCoord = ['theta', 'mu', 'lambda'][prolateCoordNum]
        xiCoord = xiCoordNum + 1
        if prolateCoord == 'theta' and xiCoord == 1:
            return 2.0 * pi / float(self.numElements[0])
        elif prolateCoord == 'mu' and xiCoord == 2:
            return self.cutoffAngle / float(self.numElements[1])
        elif prolateCoord == 'lambda' and xiCoord == 3:
            return (self.externalLambda - self.internalLambda) / float(self.numElements[2])
        else:
            # All other derivatives are zero, as the xi coordinates
            # are aligned with the prolate speroidal coordinates
            return 0.0

    def constrainedNodes(self):
        """ Return sets of nodes that are constrained to be equal.
            These are found at the heart apex.
        """
        includeQuadratic = 'quadratic' in self.interpolations
        if includeQuadratic:
            elFactor = 2
        else:
            elFactor = 1
        nodeSets = []
        # Loop radially
        for k in xrange(self.numElements[2] * elFactor + 1):
            nodeSets.append([])
            # Longitudinal position at apex:
            j = 0
            # Loop circumferentially
            for i in xrange(self.numElements[0] * elFactor):
                nodeNumber = self.nodeAtPosition((i, j, k))
                nodeSets[-1].append(nodeNumber)
        return nodeSets


def volume(f, lambdaRange, muRange, thetaRange):
    """ Compute the volume within a prolate-spheroid region
    """
    from scipy.integrate import tplquad

    def jacobian(theta, mu, lmbda):
        """ Compute sqrt(det(G_ij)), where G_ij is the covariant
            metric tensor for the transformation from prolate-spheroid
            coordinates to rectangular cartesian
        """
        return f ** 3 * sinh(lmbda) * sin(mu) * (sinh(lmbda) ** 2 + sin(mu) ** 2)

    # Note: Function to integrate has arguments reversed compared
    # to order of limits in tplquad arguments
    vol, err = tplquad(jacobian,
        lambdaRange[0], lambdaRange[1],
        lambda l: muRange[0], lambda _: muRange[1],
        lambda l, m: thetaRange[0], lambda l, m: thetaRange[1],
        epsabs=1.0e-10, epsrel=1.0e-10)

    return vol


def wallThickness(f, mu, lambdaInternal, lambdaExternal):
    """ Compute the wall thickness at the given mu position,
        from lambdaRange[0] to lambdaRange[1]
    """
    internal = xyz(f, lambdaInternal, mu, 0.0)
    external = xyz(f, lambdaExternal, mu, 0.0)
    return np.linalg.norm(external - internal)


def xyz(focus, lmbda, mu, theta):
    """ Convert from prolate-spheroid coordinates to rectangular cartesian
    """
    x = focus * sinh(lmbda) * sin(mu) * cos(theta)
    y = focus * sinh(lmbda) * sin(mu) * sin(theta)
    z = focus * cosh(lmbda) * cos(mu)
    return np.array([x, y, z])
