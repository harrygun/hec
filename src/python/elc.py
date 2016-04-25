import scipy as sp
import numpy as np
import pynbody as pn
import sympy.mpmath as mpmath
import pylab as pl

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
        #print 'calculating RD at', idx[0], idx[1], i

        #b.append(2./3.*(mpmath.elliprd(ai[idx[0]], ai[idx[1]], ai[i])*np.prod(ai)-1.))
        b.append(np.float(2./3.*mpmath.elliprd(ai[idx[0]]**2., ai[idx[1]]**2., ai[i]**2.)*np.prod(ai)))

    return b






def dynamic_elc(p, a, var, other_args):
    ''' ->> ellipsoidal collapse dynamics to be solved <<- 
        p:    prog controller, include cosmological routines 
	var:  ODE integrator variables 
    '''
    #l1, l2, l3, l1p, l2p, l3p=var
    #print 'inside dynamic_elc: ', var.shape, len(other_args), other_args

    if len(var)!=6: 
        raise Exception

    li=np.array(other_args)
    a_i, da_i = var[:3], var[3:]

    z=1./a-1.
    omz=p.pk.om0z(z)
    olz=p.pk.omLz(z)
    #mH=misc.mHz(p, z)
    qt=omz/2.-olz

    # ->> 
    dnr=a**3./np.prod(a_i)-1.
    bj=bi(a_i)
    lambda_ext=p.pk.D1(z)*(li-np.sum(li)/3.)

    dyn_ai = lambda idx: (1.+qt)*da_i[idx]+(olz-omz/2.)*a_i[idx]-3./2.*omz*\
                        (bj[idx]*dnr/2.+lambda_ext[idx])*a_i[idx]  

    # ->> derivative <<- #
    ai_p = da_i
    ai_pp= [dyn_ai(i) for i in range(3)]

    return ai_p+ai_pp






def get_elliptraj_one(p, a, lambda_i):
    # ->> ODE solver for ONE given initial condition <<- #
    # lambda_i: initial eigenvalues evaluated at z=0, so initial value = D(t_i) lambda_i
    # R:        initial Lagrangian size of the patch (NOT using)

    #->> set initial conditions <<- #
    #raise Exception('set IC for dot{a}.')
    a0=a[0] 
    z0=1./a0-1.
    D0=p.pk.D1(z0)
    f0=p.pk.f(z0)

    # ->> initial condition <<- #
    a_i = [a0*(1.-D0*lambda_i[i]) for i in range(3)]
    da_i = [a_i[i]-a0*D0*f0*lambda_i[i] for i in range(3)] 
    varl_0=a_i+da_i

    # ->> arguments <<- #
    other_args=list(lambda_i)
    all_args=tuple([p]+other_args)

    #->> return the result from ODE solver <<- #
    return myode.ode_solver_general(dynamic_elc, a0, varl_0, a, args=all_args)




def shape_to_eigval(F, e, p):
    # ->> convert {F, e, p} to eigenvalues {lambda_i} <<- #
    
    l3=F/3.*(1.+3.*e+p)
    l2=F/3.*(1.-2.*p)
    l1=F/3.*(1.-3.*e+p)
    
    return [l1, l2, l3]



def eigval_to_shape(l1, l2, l3):
    #->> convert {l1, l2, l3} into {F, e, p} <<- #
    F=l1+l2+l3 
    e=(l1-l3)/2./F
    p=(l1-2.*l2+l3)/2./F

    return [F, e, p]




'''
def get_elliptraj(a, dynvar, var_type='nu_e_p'):
    # main routine to get ellipsoidal collapse trajectories #

    if var_type=='nu_e_p':
        #->> change to eigenvalues <<- #
	F, e, p = dynvar

        l3=F/3.*(1.+3.*e+p)
        l2=F/3.*(1.-2.*p)
        l1=F/3.*(1.-3.*e+p)

        varl_0=[l1, l2, l3]
    else:
        varl_0=dynvar


    print 'varl_0 shape:', varl_0.shape
    traj=np.zeros(l1.shape)
     
    #->> 
    #for i in range(  ):
        #traj[i]=get_elliptraj_one(p, a, varl_0)

    return  
'''





''' ->> some testing routines <<- '''
def elltraj_test(p, a):

    #F, ee, pp=1., 0.5, 0.2
    F, ee, pp=1., 0., 0.
    l=shape_to_eigval(F, ee, pp)
    print 'testing lambda:', l

    traj=get_elliptraj_one(p, a, l)
    print 'traj:', traj
    print 'traj shape:', traj.shape
 

    eigval=eigval_to_shape(traj[:,0], traj[:,1], traj[:,2])
    print 'eigval', eigval

    # ->> plot <<- #
    pl.plot(a, eigval[0])    
    pl.show()
     

    return






