from __future__ import division, print_function

import numpy as np
from scipy.constants.codata import physical_constants

import xraylib

def energy(crystal_str,hkl,deg):
    '''
    Calculates the photon energy of the Bragg reflection using kinematic
    approximation for given crystal, reflection, and Bragg angle. Return
    energy in keV.
    '''
    crystal = xraylib.Crystal_GetCrystal(crystal_str)

    hc = physical_constants['Planck constant in eV s'][0]*physical_constants['speed of light in vacuum'][0]*1e6 #in keV*nm
    d = xraylib.Crystal_dSpacing(crystal,*hkl)*1e-1 #in nm

    th=np.radians(deg)

    return hc/(2*d*np.sin(th)) #in keV

def angle(crystal_str,hkl,energy):
    '''
    Calculates the Bragg angle for given photon energy using kinematic
    approximation for given crystal, reflection, and Bragg angle. energy is
    to be given in keV, returns angle in degrees.
    '''
    #constants
    crystal = xraylib.Crystal_GetCrystal(crystal_str)

    hc = physical_constants['Planck constant in eV s'][0]*physical_constants['speed of light in vacuum'][0]*1e6 #in keV*nm
    d = xraylib.Crystal_dSpacing(crystal,*hkl)*1e-1 #in nm

    if not hc/(2*d*energy) > 1:
        th0 = np.arcsin(hc/(2*d*energy))
    else:
        print('Given energy below the backscattering energy!')
        print('Setting theta to 90 deg.')
        return 90.0

    return np.degrees(th0)
