import numpy as np
import cosmology.frw as frw



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





