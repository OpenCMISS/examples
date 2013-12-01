#!/usr/bin/env python

import sys, os
sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

import numpy
import math
import scipy
from scipy.special import jn, jn_zeros



def PoiseuilleAxialVelocity(pressure,viscosity,length,r,R):
    ''' returns analytic solution to a parabolic velocity profile based
    for Poiseuille flow'''

    uAxial = -pressure/(4.*viscosity*length)*(R**2-r**2)
    return(uAxial);


def WomersleyAxialVelocity(t,pOffset,amplitude,R,r,period,viscosity,alpha,length):
    """ Computes analytic value for axial velocity assuming a Womersley profile """

    zeroTolerance = 1e-6
    omega = 2.*math.pi/period
    gamma = 1j**(3./2.)*alpha
    uWomComplex = 1j*(amplitude*(R**2))/(viscosity*(alpha**2))*(1-(jn(0,gamma*r/R))/(jn(0,gamma)))*math.exp(1)**(1j*omega*t)
    uAxial = uWomComplex.real
    if abs(pOffset) > zeroTolerance:
        PoiseuilleAxialVelocity(uOffset,pOffset,viscosity,length,r,R)
        uAxial += uOffset
    return(uAxial);

