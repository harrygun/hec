import scipy as sp
import numpy as np
import pynbody as pn
import sympy.mpmath as mpmath

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import cosmology.power  as power
import genscript.read as read
import genscript.myplot as mpl
import genscript.read as rd
import misc

import ode_solver as myode




def bi(ai):
    # ->> t: time list <<- #

    #b=np.zeros(3)
    idxall=np.arange(3)

    b=[]
    for i in range(3):
        idx=np.delete(idxall, i)
        print 'calculating RD at', idx[0], idx[1], i

        #b.append(2./3.*(mpmath.elliprd(ai[idx[0]], ai[idx[1]], ai[i])*np.prod(ai)-1.))
        b.append(2./3.*mpmath.elliprd(ai[idx[0]]**2., ai[idx[1]]**2., ai[i]**2.)*np.prod(ai))

    return b






def dynamic_elc(p, a, var):
    ''' ->> ellipsoidal collapse dynamics to be solved <<- 
        p:    prog controller, include cosmological routines 
	var:  ODE integrator variables 
    '''
    #l1, l2, l3, l1p, l2p, l3p=var
    l = var[:3]
    lp= var[3:]

    z=1./a-1.
    omz=p.pk.om0z(z)
    olz=p.pk.omLz(z)
    #mH=misc.mHz(p, z)
    qt=omz/2.-olz


    # ->> 
    dnr=a**3./np.prod(l)-1.
    bj=bi(l)
    l_ext=p.pk.D1(z)*np.([ ]) 

    #dyn_l=lambda idx: -mH**2*(olz+omz*l[idx]*(1./3.+dnr/3.+ dnr*bj[idx]/2.+l_ext[idx]))
    dyn_l = lambda idx: (1.+qt)*lp[idx]+(olz-omz/2.)*l[idx] \
                        -3./2.*omz*(bj[idx]*dnr[idx]/2.+l_ext[idx])*l[idx]  

    # ->> derivative <<- #
    l_p = lp
    l_pp= [dyn_l(i) for i in range(3)]

    return l_p+l_pp






def get_elliptraj_one(a, dynvar):
    # ->> ODE solver for ONE given initial condition <<- #

     




    return





def get_elliptraj(a, dynvar, var_type='nu_e_p'):
    ''' main routine to get ellipsoidal collapse trajectories '''

    if var_type=='nu_e_p':
        #->> change to eigenvalues <<- #
	F, e, p = dynvar

        l3=F/3.*(1.+3.*e+p)
        l2=F/3.*(1.-2.*p)
        l1=F/3.*(1.-3.*e+p)

        dynvar=[l1, l2, l3]


    print 'dynvar shape:', dynvar.shape
     
    #->> 
    for i in range(len(dynvar.shape[0])):
        for j in range(len(dynvar.shape[1])):
            for k in range(len(dynvar.shape[2])):
	        #->> 




    return  





