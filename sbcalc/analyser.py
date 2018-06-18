from __future__ import division, print_function
import numpy as np

from pyTTE import takagitaupin
from pyTTE.deformation import isotropic_plate, anisotropic_plate

from . import lateral_deformation
from . import crystal
from .johann_error import johann_error

class Analyser:
    '''
    Analyser object contains all the necessary information needed compute
    the analyser diffraction curve.
    '''

    def __init__(self,crystal_str,hkl,R_bend,thickness,**kwargs):
        '''
        crystal_str = String of the crystal material e.g. 'Si'
        hkl = Miller indeces [h,k,l]
        R_bend = bending radius in meters
        thickness = wafer_thickness in microns

        Other arguments depend on the type of analyser used. The type of the
        analyser is inferred from the set of key word arguments given as follows:

        CIRCULAR WAFER:
        diameter = diameter of the wafer in mm

        RECTANGULAR WAFER:
        a, b = side lengths of the wafer in mm

        STRIP-BENT:
        diameter = diameter of the analyser in mm
        strip_width = width of the strip in mm


        Anisotropic or isotropic theory is also chosen on the basis of given
        parameters:

        ANISOTROPIC:
        S = the compliance matrix in units 1/GPa

        ISOTROPIC:
        nu = Poisson ratio
        E = Young's modulus
        '''

        self.crystal_str = crystal_str
        self.hkl = hkl
        self.R_bend = R_bend
        self.thickness = thickness

        self.analyser_type = None
        self.scantype = None
        self.photon_energy = None
        self.incidence_angle = None
        self.polarization = None

        self.scan_vector = None
        self.TT_reflectivity = None
        self.lateral_shifts = None
        self.resolution_function = None

        #The effective Poisson ratios for the case of the anisotropic model
        self.effective_poisson_1D_TT = None

        #Init the analyser type
        if 'a' in kwargs and 'b' in kwargs:
            self.analyser_type = 'rectangular'
            self.a = kwargs['a']
            self.b = kwargs['b']
        elif 'diameter' in kwargs:
            if 'strip_width' in kwargs:
                self.analyser_type = 'stripbent'
                raise Exception('Strip-bent analyser not implemented yet!')
            else:
                self.analyser_type = 'circular'
                self.diameter = kwargs['diameter']
        else:
            raise Exception('Invalid analyser type input!\n'+\
                            'Check that the keyword arguments are consistent.')

        #Init the theory
        if 'S' in kwargs:
            if not self.analyser_type == 'circular':
                Exception('Anisotropic theory compatible only with circular wafer!')
            self.theory = 'anisotropic'
            self.compliance_matrix = kwargs['S']
        elif 'nu' in kwargs and 'E' in kwargs:
            self.theory = 'isotropic'
            self.poisson_ratio = kwargs['nu']
            self.young_modulus = kwargs['E']
        else:
            raise Exception('Invalid elastic constant input!\n'+\
                            'Check that the keyword arguments are consistent.')

        #compute the stress and strain fields
        if self.analyser_type == 'circular':
            if self.theory == 'anisotropic':
                self.X, self.Y, self.stress, self.strain = \
                lateral_deformation.anisotropic_circular(self.R_bend*1e3,
                                                         self.diameter,
                                                         self.compliance_matrix)
            else:
                self.X, self.Y, self.stress, self.strain = \
                lateral_deformation.isotropic_circular(self.R_bend*1e3,
                                                       self.diameter,
                                                       self.poisson_ratio,
                                                       self.young_modulus)
        elif self.analyser_type == 'rectangular':
            self.X, self.Y, self.stress, self.strain = \
            lateral_deformation.isotropic_rectangular(self.R_bend*1e3,
                                                   self.a, self.b,
                                                   self.poisson_ratio,
                                                   self.young_modulus)
        else:
            raise Exception('Strip-bent analyser not implemented yet!')

    def compute_TT(self,scan,**kwargs):
        '''
        Computes the 1D Takagi-Taupin curve for symmetric Bragg reflection.

        Input:

        scan = scan vector either in angle (arc sec) or energy (meV), depending
               whether kwarg 'energy' or 'theta' is defined, respectively.

        energy = energy of incident photons in keV
        OR
        theta = Incidence angle of photons in degrees
        polarization = 'sigma' or 'pi'. If not given 'sigma' is assumed
        '''

        if 'energy' in kwargs:
            self.scantype = 'angle'
            constant = kwargs['energy']
            self.photon_energy = kwargs['energy']
            self.incidence_angle = None
        elif 'theta' in kwargs:
            self.scantype = 'energy'
            constant = kwargs['theta']
            self.incidence_angle = kwargs['theta']
            self.photon_energy = None
        else:
            print("Warning! Either 'energy' or 'theta has to be given'.")
            return

        if 'polarization' in kwargs:
            if kwargs['polarization'] == 'sigma' or kwargs['polarization'] == 'pi':
                self.polarization = kwargs['polarization']
            else:
                print("Warning! Polarization has to be either 'sigma' or 'pi'.")
                return
        else:
                self.polarization = 'sigma'

        asymmetry = 0
        if self.theory == 'isotropic':
            ujac = isotropic_plate(self.R_bend,self.R_bend,self.poisson_ratio,self.thickness*1e-6)
        else:
            ujac = anisotropic_plate(self.R_bend,self.R_bend,self.compliance_matrix,self.thickness*1e-6)
            strain_grad = ujac(0,1-0.5*self.thickness*1e-6)[1][1]
            poisson_ratio = strain_grad/(2+strain_grad)
            self.effective_poisson_1D_TT = poisson_ratio
            print('Effective Poisson ratio in 1D TT-solving: '+str(poisson_ratio))

        R,T = takagitaupin('energy',scan,constant,self.polarization,self.crystal_str,self.hkl,asymmetry,self.thickness,ujac)

        self.scan_vector = scan
        self.TT_reflectivity = R

    def compute_lateral_shifts(self,include_johann_error=True):
        '''
        Computes the energy or angle shifts owing to the lateral strain.
        See doi:10.1107/S1600576716010402 for details.
        Determines the shift domain on the basis of 1D Takagi-Taupin run.
        The energy shits are calculated in meV, angle in arc sec
        '''

        if self.scantype == 'energy':
            E0 = crystal.energy(self.crystal_str,self.hkl,self.incidence_angle)*1e6
            self.lateral_shifts = -E0*self.strain['zz']
            if include_johann_error:
                jerr = johann_error(self.X,self.R_bend*1e3,self.incidence_angle,E0)
                self.lateral_shifts = self.lateral_shifts + jerr
        elif self.scantype == 'angle':
            th = crystal.angle(self.crystal_str,self.hkl,self.photon_energy)
            self.lateral_shifts = -self.strain['zz']*np.tan(np.radians(th))*206264.80625
            if include_johann_error:
                jerr = johann_error(self.X,self.R_bend*1e3,th)
                self.lateral_shifts = self.lateral_shifts + jerr*206264.80625
        else:
            print('Warning! No scantype found! Call compute_TT() first.')

    def compute_resolution_function(self,mask=None):
        '''
        Computes the resolution function of the crystal analyser by convolving
        the 1D TT curve with the calculated lateral shifts. If not calculated
        before, the function will call compute_lateral_shifts() with Johan error
        included.

        mask = Boolean matrix with same shape as the lateral_shift matrix. The
               area masked by 'False' entries are outcluded from the convolution.
        '''

        #apply possible mask and remove Nan values
        if not mask == None:
            shifts = self.lateral_shifts[mask]
        else:
            shifts = self.lateral_shifts
        shifts = shifts[np.logical_not(np.isnan(shifts))]

        interp = np.interp

        #create a new energy or angle axis for the convolved spectrum
        tt_axis = self.scan_vector
        I_tt = self.TT_reflectivity
        N = np.ceil((np.max(tt_axis)+np.max(shifts)-(np.min(tt_axis)+np.min(shifts)))/(np.max(tt_axis)-np.min(tt_axis))*tt_axis.size)
        new_axis = np.linspace(np.min(tt_axis)+np.min(shifts),np.max(tt_axis)+np.max(shifts),N)
        y_new = np.zeros(new_axis.shape)

        for shift in shifts.reshape((-1,)):
            y_new = y_new + interp(new_axis,tt_axis+shift,I_tt,left=0,right=0)

        #Normalization to the number of grid points
        y_new = y_new/np.logical_not(np.isnan(self.lateral_shifts))[:].sum()

        self.resolution_function = (new_axis,y_new)

    def convolve(self,*args):
        '''
        Convolves the computed resolution function with given function
        (two arguments x,y, where x has to be in units of meV or arc sec,
        matching the resolution function scale) OR with a gaussian (one
        argument, fwhm in the units of meV or arcsec)

        Returns the convolved resolution function
        '''

        if len(args) == 1:
            def gauss(x,x0,fwhm):
                return np.exp(-np.log(2)*((x-x0)/fwhm*2)**2)

            Ibw = gauss(self.resolution_function[0],np.mean(self.resolution_function[0]),args[0])

            x_conv = self.resolution_function[0].copy()
            y_conv = np.convolve(self.resolution_function[1],Ibw,'same')

            #normalize to same area as the original
            y_conv = y_conv*np.trapz(self.resolution_function[1],x_conv)/np.trapz(y_conv,x_conv)

            return x_conv, y_conv
        else:
            x = args[0]; y = args[1]
            x = x-np.mean(x)+np.mean(self.resolution_function[0])

            Ibw = np.interp(self.resolution_function[0],x,y,left=0,right=0)

            x_conv = self.resolution_function[0].copy()
            y_conv = np.convolve(self.resolution_function[1],Ibw,'same')

            #normalize to same area as the original
            y_conv = y_conv*np.trapz(self.resolution_function[1],x_conv)/np.trapz(y_conv,x_conv)

            return x_conv, y_conv
