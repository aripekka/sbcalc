from __future__ import division, print_function
import numpy as np

def johann_error(X,R,th,energy=None):
    '''
    Calculates the Johann error for given x-coordinates (dispersive direction)
    in angle or energy scale, depending whether energy is given.

    Input:
    X = coordinates along the dispersion direction (X=0 coincides with the
        center of the crystal)
    R = Bending radius
    th = Incidence (Bragg's) angle in degrees
    energy = energy of photons in units of preference

    The units of X and R don't matter as long as they are same.
    '''
    if energy == None:
        return X**2/(2*R**2*np.tan(np.radians(th)))
    else:
        return -X**2/(2*R**2*np.tan(np.radians(th))**2)*energy
