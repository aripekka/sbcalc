import matplotlib.pyplot as plt
import numpy as np
from analyser import Analyser

def test_circular_analyser():
    #isotropic model
    anal_iso = Analyser('Si',[6,6,0],1,150,diameter=100,E=165,nu=0.27)
    anal_iso.compute_TT(np.linspace(-600,600,100),theta=88)

    #anisotropic model
    S = np.zeros((6,6))

    #The elastic matrix for isotropic crystal
    S[0,0] = 1; S[1,1] = 1; S[2,2] = 1

    S[0,1] = -0.27; S[0,2] = -0.27; S[1,2] = -0.27
    S[1,0] = -0.27; S[2,0] = -0.27; S[2,1] = -0.27

    S[3,3] = 2*(1+0.27); S[4,4] = 2*(1+0.27); S[5,5] = 2*(1+0.27)

    S = S/165

    anal_aniso = Analyser('Si',[6,6,0],1,150,diameter=100,S=S)
    anal_aniso.compute_TT(np.linspace(-600,600,100),theta=88)


    plt.plot(anal_iso.scan_vector,anal_iso.TT_reflectivity,label='Isotropic')
    plt.plot(anal_aniso.scan_vector,anal_aniso.TT_reflectivity,label='Anisotropic')
    plt.legend()

    #lateral strain
    plt.figure()
    plt.subplot(311)

    anal_iso.compute_lateral_shifts(include_johann_error=False)
    plt.pcolor(anal_iso.X,anal_iso.Y,anal_iso.lateral_shifts)

    plt.subplot(312)
    old_lateral = anal_iso.lateral_shifts
    anal_iso.compute_lateral_shifts(include_johann_error=True)
    plt.pcolor(anal_iso.X,anal_iso.Y,anal_iso.lateral_shifts)

    plt.subplot(313)
    plt.pcolor(anal_iso.X,anal_iso.Y,anal_iso.lateral_shifts-old_lateral)

    #resolution function
    plt.figure()
    anal_iso.compute_resolution_function()
    plt.plot(*anal_iso.resolution_function)

    #convolution
    plt.figure()
    conv = anal_iso.convolve(300)
    plt.plot(*anal_iso.resolution_function)
    plt.plot(*conv)

    plt.show()
