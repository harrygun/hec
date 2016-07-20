import numpy as np
import cosmology.frw as frw
import scipy.integrate as intg


_default_Hz_normalization_='E(z)'
#_default_Hz_normalization_='km/s/Mpc'




def Hz(p, z, normalize_type=_default_Hz_normalization_):
    ''' ->> set default normalization  <<- '''

    if hasattr(p, 'Hz_normalization'):
        normalize_type=p.Hz_normalization

    return p.dist.H(z, normalize_type=normalize_type)

def mHz(p, z, normalize_type=_default_Hz_normalization_):
    ''' ->> set default normalization  <<- '''

    if hasattr(p, 'Hz_normalization'):
        normalize_type=p.Hz_normalization

    return p.dist.mH(z, normalize_type=normalize_type)






''' ->> variance sig2 <<- '''

def wd(k, R, window_type):
    if window_type=='Gauss':
        return np.exp(-(k*R)**2./2.)
    elif window_type=='Tophat':
        kr=k*R
        return 3.*(np.sin(kr)-kr*np.cos(kr) )/kr**3 

def sigma2_integrate(p, R=0., window_type='Gauss'):
    sf = lambda lgk: p.pk(10.**lgk)*wd(10.**lgk, R, window_type)**2*10.**(3.*lgk)\
                        *np.log(10.)/2./np.pi**2

    s2=intg.quad(sf, -4., 3., limit=np.int(1e8) )
    return s2


def sig2(p, z, R=0., window_type='Gauss'):
    s2=sigma2_integrate(p, R=R, window_type=window_type)
    return  p.pk.D1(z)**2*s2[0]

