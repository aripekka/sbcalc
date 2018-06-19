from __future__ import division, print_function
import numpy as np

N_GRID = 100

def isotropic_circular(R,L,nu=0.27,E=165):
    '''
    Computes the lateral stress and strain tensors for isotropic circular
    SBCA.

    Input:
    R = bending radius
    L = analyser diameter
    nu = Poisson's ratio
    E = Young's modulus

    The units of R and L do not matter, as long as they are the same. Same for
    E but the choice fixes the units of the computed stress tensor

    Output:
    X, Y = 2D matrices giving the coordinates of the surface points
    stress = The components of the stress tensor
    strain = The components of the strain tensor

    Note: stress and strain matrices contain nan for the points outside the
    crystal surface.
    '''
    #generate the grid shifts
    x=np.linspace(-L/2,L/2,N_GRID)
    X,Y=np.meshgrid(x,x)

    stress = {}
    strain = {}

    stress['xx'] = -E/(16*R**2)*(X**2 + 3*Y**2 -L**2/4)
    stress['yy'] = -E/(16*R**2)*(3*X**2 + Y**2 -L**2/4)
    stress['xy'] = E/(8*R**2)*X*Y
    stress['yx'] = stress['xy']

    stress['rr'] = E/(16*R**2)*(L**2/4-X**2-Y**2)
    stress['phiphi'] = E/(16*R**2)*(L**2/4-3*X**2-3*Y**2)
    stress['rphi'] = np.zeros(X.shape)
    stress['phir'] = stress['rphi']

    strain['xx'] = ((1-nu)*L**2/4-(1-3*nu)*X**2-(3-nu)*Y**2)/(16*R**2)
    strain['yy'] = ((1-nu)*L**2/4-(1-3*nu)*Y**2-(3-nu)*X**2)/(16*R**2)
    strain['xy'] = (1+nu)/(8*R**2)*X*Y
    strain['yx'] = strain['xy']

    strain['rr'] = 1/(16*R**2)*((1-nu)*L**2/4-(1-3*nu)*(X**2+Y**2))
    strain['phiphi'] = 1/(16*R**2)*((1-nu)*L**2/4-(3-nu)*(X**2+Y**2))
    strain['rphi'] = np.zeros(X.shape)
    strain['phir'] = strain['rphi']

    strain['zz'] = nu/(4*R**2)*(X**2+Y**2-L**2/8)

    for k in stress:
        stress[k][X**2+Y**2 > L**2/4] = np.nan
    for k in strain:
        strain[k][X**2+Y**2 > L**2/4] = np.nan

    return X,Y,stress,strain

def anisotropic_circular(R,L,S):
    '''
    Computes the lateral stress and strain tensors for anisotropic circular
    SBCA.

    Input:
    R = bending radius
    L = analyser diameter
    S = compliance matrix

    The units of R and L do not matter, as long as they are the same. Same for
    S but the choice fixes the units of the computed stress tensor

    Output:
    X, Y = 2D matrices giving the coordinates of the surface points
    stress = The components of the stress tensor
    strain = The components of the strain tensor

    Note: stress and strain matrices contain nan for the points outside the
    crystal surface.
    '''

    x=np.linspace(-L/2,L/2,N_GRID)
    X,Y=np.meshgrid(x,x)

    r_squared = X**2+Y**2
    phi = np.arctan2(Y,X)

    stress = {}
    strain = {}

    D = 1/(2*R**2*(3*(S[0,0]+S[1,1])+2*S[0,1]+S[5,5]))

    stress['xx'] = D*(L**2/4-X**2-3*Y**2)
    stress['yy'] = D*(L**2/4-3*X**2-Y**2)
    stress['xy'] = 2*D*X*Y
    stress['yx'] = stress['xy']

    stress['rr'] = D*(L**2/4-r_squared)
    stress['phiphi'] = D*(L**2/4-3*r_squared)
    stress['rphi'] = np.zeros(X.shape)
    stress['phir'] = stress['rphi']

    #shorthand notation
    uzzaux1 = (S[2,0]+S[2,1])*L**2/4
    uzzaux2 = 2*(S[2,0]+S[2,1])
    uzzaux3 = np.sqrt((S[2,1]-S[2,0])**2+S[2,5]**2)

    beta = np.arctan2(S[2,5],(S[2,1]-S[2,0]))

    strain['zz'] = D*(uzzaux1 - (uzzaux2+uzzaux3*np.sin(2*phi+beta))*r_squared)

    eff_poisson = -4*(S[2,0]+S[2,1])/(3*(S[0,0]+S[1,1])+2*S[0,1]+S[5,5])
    print('Effective Poisson ratio in lateral strain: ', str(eff_poisson))

    for k in stress:
        stress[k][X**2+Y**2 > L**2/4] = np.nan
    for k in strain:
        strain[k][X**2+Y**2 > L**2/4] = np.nan

    return X,Y,stress,strain

def isotropic_rectangular(R,a,b,nu=0.27,E=165):
    '''
    Computes the lateral stress and strain tensors for isotropic rectangular
    SBCA.

    Input:
    R = bending radius
    a, b = side lengths of the rectangle
    nu = Poisson's ratio
    E = Young's modulus

    The units of R and L do not matter, as long as they are the same. Same for
    E but the choice fixes the units of the computed stress tensor

    Output:
    X, Y = 2D matrices giving the coordinates of the surface points
           X-axis is parallel with the side a and y-axis with b
    stress = The components of the stress tensor
    strain = The components of the strain tensor

    Note: stress and strain matrices contain nan for the points outside the
    crystal surface.
    '''
    #resize grid points according to the linear dimensions so the separation
    #of the grid points is equal in x and y directions
    if a > b:
        x=np.linspace(-a/2,a/2,N_GRID)
        y=np.linspace(-b/2,b/2,int(np.round(N_GRID*b/a)))
    else:
        y=np.linspace(-b/2,b/2,N_GRID)
        x=np.linspace(-a/2,a/2,int(np.round(N_GRID*a/b)))

    X,Y=np.meshgrid(x,y)

    stress = {}
    strain = {}

    g = 8 + 10*((a/b)**2+(b/a)**2) + (1-nu)*((a/b)**2-(b/a)**2)**2

    stress['xx'] = E/(g*R**2) * (a**2/12-X**2 + ((1+nu)/2 + 5*(a/b)**2 + (1-nu)/2*(a/b)**4)*(b**2/12-Y**2))
    stress['yy'] = E/(g*R**2) * (b**2/12-Y**2 + ((1+nu)/2 + 5*(b/a)**2 + (1-nu)/2*(b/a)**2)*(a**2/12-X**2))
    stress['xy'] = 2*E/(g*R**2)*X*Y

    strain['zz'] = nu/(g*R**2) * (((3+nu)/2+5*(b/a)**2+(1-nu)/2*(b/a)**4)*(X**2 - a**2/12)+\
                                  ((3+nu)/2+5*(a/b)**2+(1-nu)/2*(a/b)**4)*(Y**2 - b**2/12))

    return X,Y,stress,strain
