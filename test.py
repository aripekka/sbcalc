import matplotlib.pyplot as plt
import numpy as np
from sbcalc import Analyser


def transform_to_polar(X,Y,Txx,Txy,Tyy):
    PHI = np.arctan2(Y,X)
    Urr = np.cos(PHI)**2*Txx + 2*np.cos(PHI)*np.sin(PHI)*Txy + np.sin(PHI)**2*Tyy
    Urphi = -np.cos(PHI)*np.sin(PHI)*Txx + (np.cos(PHI)**2 - np.sin(PHI)**2)*Txy + np.cos(PHI)*np.sin(PHI)*Tyy
    Uphiphi = np.sin(PHI)**2*Txx - 2*np.cos(PHI)*np.sin(PHI)*Txy + np.cos(PHI)**2*Tyy

    return Urr, Urphi, Uphiphi

def test_isotropic_circular():

    anal_iso = Analyser('Si',[6,6,0],1,150,diameter=100,E=165,nu=0.27)

    X,Y,stress,strain = anal_iso.X, anal_iso.Y, anal_iso.stress, anal_iso.strain

    plt.figure()
    plt.subplot(231)
    plt.pcolor(X,Y,stress['xx'])
    plt.title('$\sigma_{xx}$')
    plt.subplot(232)
    plt.pcolor(X,Y,stress['xy'])
    plt.title('$\sigma_{xy}$')
    plt.subplot(233)
    plt.pcolor(X,Y,stress['yy'])
    plt.title('$\sigma_{yy}$')
    plt.subplot(234)
    plt.pcolor(X,Y,stress['rr'])
    plt.title('$\sigma_{rr}$')
    plt.subplot(235)
    plt.pcolor(X,Y,stress['rphi'])
    plt.title('$\sigma_{r\phi}$')
    plt.subplot(236)
    plt.pcolor(X,Y,stress['phiphi'])
    plt.title('$\sigma_{\phi\phi}$')

    plt.figure()
    plt.subplot(231)
    plt.pcolor(X,Y,strain['xx'])
    plt.title('$u_{xx}$')
    plt.subplot(232)
    plt.pcolor(X,Y,strain['xy'])
    plt.title('$u_{xy}$')
    plt.subplot(233)
    plt.pcolor(X,Y,strain['yy'])
    plt.title('$u_{yy}$')
    plt.subplot(234)
    plt.pcolor(X,Y,strain['rr'])
    plt.title('$u_{rr}$')
    plt.subplot(235)
    plt.pcolor(X,Y,strain['rphi'])
    plt.title('$u_{r\phi}$')
    plt.subplot(236)
    plt.pcolor(X,Y,strain['phiphi'])
    plt.title('$u_{\phi\phi}$')

    #test tensor transfromations
    srr, srphi, sphiphi = transform_to_polar(X,Y,stress['xx'],stress['xy'],stress['yy'])
    urr, urphi, uphiphi = transform_to_polar(X,Y,strain['xx'],strain['xy'],strain['yy'])

    plt.figure()
    plt.subplot(231)
    plt.pcolor(X,Y,strain['rr']-urr)
    plt.title('$\Delta u_{rr}$')
    plt.subplot(232)
    plt.pcolor(X,Y,strain['rphi']-urphi)
    plt.title('$\Delta u_{r\phi}$')
    plt.subplot(233)
    plt.pcolor(X,Y,strain['phiphi']-uphiphi)
    plt.title('$\Delta u_{\phi\phi}$')
    plt.subplot(234)
    plt.pcolor(X,Y,stress['rr']-srr)
    plt.title('$\Delta \sigma_{rr}$')
    plt.subplot(235)
    plt.pcolor(X,Y,stress['rphi']-srphi)
    plt.title('$\Delta \sigma_{r\phi}$')
    plt.subplot(236)
    plt.pcolor(X,Y,stress['phiphi']-sphiphi)
    plt.title('$\Delta \sigma_{\phi\phi}$')


def test_anisotropic_circular():

    S = np.zeros((6,6))

    #The elastic matrix for isotropic crystal
    S[0,0] = 1
    S[1,1] = 1
    S[2,2] = 1

    S[0,1] = -0.27
    S[0,2] = -0.27
    S[1,2] = -0.27
    S[1,0] = -0.27
    S[2,0] = -0.27
    S[2,1] = -0.27

    S[3,3] = 2*(1+0.27)
    S[4,4] = 2*(1+0.27)
    S[5,5] = 2*(1+0.27)

    S = S/165

    anal_iso = Analyser('Si',[6,6,0],1,150,diameter=100,E=165,nu=0.27)
    anal_aniso = Analyser('Si',[6,6,0],1,150,diameter=100,S=S)
    X,Y,stress,strain = anal_iso.X, anal_iso.Y, anal_iso.stress, anal_iso.strain
    X,Y,stress2,strain2 = anal_aniso.X, anal_aniso.Y, anal_aniso.stress, anal_aniso.strain

    plt.figure()
    plt.subplot(231)
    plt.pcolor(X,Y,stress2['xx']-stress['xx'])
    plt.title('$\Delta \sigma_{xx}$')
    plt.subplot(232)
    plt.pcolor(X,Y,stress2['xy']-stress['xy'])
    plt.title('$\Delta \sigma_{xy}$')
    plt.subplot(233)
    plt.pcolor(X,Y,stress2['yy']-stress['yy'])
    plt.title('$\Delta \sigma_{yy}$')
    plt.subplot(234)
    plt.pcolor(X,Y,stress2['rr']-stress['rr'])
    plt.title('$\Delta \sigma_{rr}$')
    plt.subplot(235)
    plt.pcolor(X,Y,stress2['rphi']-stress['rphi'])
    plt.title('$\Delta \sigma_{r\phi}$')
    plt.subplot(236)
    plt.pcolor(X,Y,stress2['phiphi']-stress['phiphi'])
    plt.title('$\Delta \sigma_{\phi\phi}$')


    plt.figure()
    plt.subplot(131)
    plt.pcolor(X,Y,strain['zz'])
    plt.title('Isotropic $u_{zz}$')
    plt.colorbar()
    plt.subplot(132)
    plt.pcolor(X,Y,strain2['zz'])
    plt.title('Ansotropic $u_{zz}$')
    plt.colorbar()
    plt.subplot(133)
    plt.pcolor(X,Y,strain2['zz']-strain['zz'])
    plt.title('$\Delta u_{zz}$')

    plt.colorbar()
    plt.show()
