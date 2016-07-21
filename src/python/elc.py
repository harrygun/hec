import scipy as sp
import numpy as np
import pynbody as pn
import sympy.mpmath as mpmath
import scipy.integrate as integ
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


_mx_step_=1000

# ->> freezing fraction <<- #
f_freeze=0.18




def bi(ai):
    # ->> t: time list <<- #

    #b=np.zeros(3)
    idxall=np.arange(3)

    b=[]
    for i in range(3):
        idx=np.delete(idxall, i)
        #print 'calculating RD at', idx[0], idx[1], i
        #b.append(2./3.*(mpmath.elliprd(ai[idx[0]], ai[idx[1]], ai[i])*np.prod(ai)-1.))

        b.append(np.float64(2./3.*mpmath.elliprd(ai[idx[0]]**2., ai[idx[1]]**2., ai[i]**2.)*np.prod(ai)))

    return b




def dynamic_elc(var, lna, p, other_args):
    ''' ->> ellipsoidal collapse dynamics to be solved <<- 
        p:    prog controller, include cosmological routines 
	var:  ODE integrator variables 
    '''
    a=np.exp(lna)

    if len(var)!=6: 
        raise Exception

    li=np.array(other_args[:3])
    freeze=other_args[3:]

    a_i, da_i = list(var[:3]), list(var[3:])

    z=1./a-1.
    omz=p.pk.om0z(z)
    olz=p.pk.omLz(z)
    qt=omz/2.-olz

    # ->> 
    dnr=a**3./np.prod(a_i)-1.
    bj=bi(a_i)
    lambda_ext=p.pk.D1(z)*(li-np.sum(li)/3.)

    dyn_ai = lambda idx: (1.+qt)*da_i[idx]+(olz-omz/2.)*a_i[idx]-3./2.*omz*\
                        (bj[idx]*dnr/2.+lambda_ext[idx])*a_i[idx]  
    #dyn_ai = lambda idx: (1.+qt)*da_i[idx]+(olz-omz/2.)*a_i[idx]-3./2.*omz*(dnr/3.)*a_i[idx]  

    # ->> derivative <<- #
    ai_p = da_i
    ai_pp= [dyn_ai(i) for i in range(3)]

    for i in range(3):
        if (a_i[i]<=f_freeze*a):
            freeze[i]=True

        if (freeze[i]==True):
            ai_p[i]=0.
            ai_pp[i]=0.

    return ai_p+ai_pp






def get_elliptraj_one(p, a, lambda_i):
    # ->> ODE solver for ONE given initial condition <<- #
    # lambda_i: initial eigenvalues evaluated at z=0, so initial value = D(t_i) lambda_i
    # R:        initial Lagrangian size of the patch (NOT using)


    #->> set initial conditions <<- #
    a0=a[0] 
    z0=1./a0-1.
    D0=p.pk.D1(z0)
    f0=p.pk.f(z0)

    #print 'z0, D0=', z0, D0


    # ->> initial condition <<- #
    a_i = [a0*(1.-D0*lambda_i[i]) for i in range(3)]
    da_i = [a_i[i]-a0*D0*f0*lambda_i[i] for i in range(3)] 
    varl_0=a_i+da_i
    #print 'IC:', varl_0

    # ->> arguments <<- #
    freeze=[False]*3
    other_args=list(lambda_i)+freeze
    all_args=tuple([p]+other_args)

    #->> return the result from ODE solver <<- #
    lna=np.log(a)

    return integ.odeint(dynamic_elc, varl_0, lna, args=(p, other_args)) #, mxstep=_mx_step_)





def shape_to_eigval(F, e, p):
    # ->> convert {F, e, p} to eigenvalues {lambda_i} <<- #
    
    l3=F/3.*(1.+3.*e+p)
    l2=F/3.*(1.-2.*p)
    l1=F/3.*(1.-3.*e+p)
    
    return [l1, l2, l3]


def delta_to_eigval_sc(F):
    l=1.-(1./(1.+F))**(1./3.)
    return [l]*3


def eigval_to_shape(l1, l2, l3):
    #->> convert {l1, l2, l3} into {F, e, p} <<- #
    F=l1+l2+l3 
    e=(l1-l3)/2./F
    p=(l1-2.*l2+l3)/2./F

    return [F, e, p]



def lambda_comving(l1, l2, l3, a):
    # ->> return comoving eigvalues <<- #
    return np.array([l1/a, l2/a, l3/a])


def lambda_to_rho(l1, l2, l3, a):
    # ->> lambda to density contrast <<- #
    return  a**3./(l1*l2*l3)-1.

def clambda_to_rho(cl1, cl2, cl3):
    # ->> comoving lambda to density contrast <<- #
    return  1./(cl1*cl2*cl3)-1.

''' --------------------------------------------------------
              ->>   some testing routines   <<- 
    --------------------------------------------------------
'''
def elltraj_test(p, a):

    #->> 
    ln_a=np.linspace(-4.6, 0., 500)
    a=np.exp(ln_a)

    F, ee, pp=0.5, 0.2, 0.1
    #F, ee, pp=1.686, 0., 0.


    l=shape_to_eigval(F, ee, pp)
    print 'testing lambda:', l

    traj=get_elliptraj_one(p, a, l)
    #print 'traj:', traj
    print 'traj shape:', traj.shape
 
    shape=eigval_to_shape(traj[:,0], traj[:,1], traj[:,2])
    #print 'shape', shape


    # ->> plot <<- #
    nplt, ncol = 3, 3
    fig, ax=mpl.mysubplots(nplt,ncol_max=ncol,subp_size=5., gap_size=0.5, return_figure=True)


    for i in range(3):
        #ax[i].plot(a, shape[i])    
        ax[i].plot(a, traj[:,i]/a, 'k-')    
        #ax[i].plot(a, traj[:,i+3], 'b-')    

    pl.show()
     

    return






