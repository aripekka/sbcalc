from sbcacl.lateral_deformation import *
import matplotlib.pyplot as plt

def test_isotropic_circular():

    X,Y,stress,strain = isotropic_circular(1,0.1,0.27,165)

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
    plt.pcolor(X,Y,strain['zz'])
    plt.title('$u_{zz}$')
    plt.colorbar()
    plt.show()

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

    X,Y,stress,strain = isotropic_circular(1,0.1,0.27,165)
    X,Y,stress2,strain2 = anisotropic_circular(1,0.1,S)

    plt.figure()
    #plt.subplot(231)
    #plt.pcolor(X,Y,stress['xx'])
    #plt.title('$\sigma_{xx}$')
    #plt.subplot(232)
    #plt.pcolor(X,Y,stress['xy'])
    #plt.title('$\sigma_{xy}$')
    #plt.subplot(233)
    #plt.pcolor(X,Y,stress['yy'])
    #plt.title('$\sigma_{yy}$')
    plt.subplot(234)
    plt.pcolor(X,Y,stress2['rr']-stress['rr'])
    plt.title('$\sigma_{rr}$')
    plt.subplot(235)
    plt.pcolor(X,Y,stress2['rphi']-stress['rphi'])
    plt.title('$\sigma_{r\phi}$')
    plt.subplot(236)
    plt.pcolor(X,Y,stress2['phiphi']-stress['phiphi'])
    plt.title('$\sigma_{\phi\phi}$')


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
