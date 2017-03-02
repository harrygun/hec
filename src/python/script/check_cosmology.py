import os
import numpy as np
import pylab as pl
import copy
import scipy.interpolate as sitp

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import genscript.myarray as mar
import cosmology.cosparameter as cospar
import cosmology.power as power

import misc as misc
import elc as elc
import trajs_generator as gtraj






param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'cosmology_parameter_fname': 'parameters/cosparameter.cfg',
    'cosmology_parameter_sec': 'Cosmology_Parameters',
    'a_init': 1e-3,
    'smooth_R': 10.,
    'smooth_type': 'Gauss', 
    'smooth_R_list_type':  'linear'
    }

prog_control={
    #'do_testing': False, 
    #-------------------------------#
    'traj_generator_type':  'testing',
    #-------------------------------#
    'e_list':  None, 
    'p_list':  None, 
    #-------------------------------#
    }





if __name__=='__main__':

    # ->> initialization <<- #
    sec='General'
    init_dict=myDict(prog_control)+myDict(param_dict)
    p=pc.prog_init(section=sec, init_cosmology=True, **init_dict)

    root='../../workspace/result/'
     
    print 'Cospar:', p.cp.H0, p.cp.h, p.cp.omem, p.cp.omeb, p.cp.omec
    print 'D1:', p.pk.D1(1e3-1), p.pk.D1(0)

    zi=1.e3-1

    cp0=copy.copy(p.cp)
    cp0.omec=0.26
    cp0.omex=1.-cp0.omec
    pk0=power.PowerSpectrum(cp0)

    cpm=copy.copy(p.cp)
    cpm.omec=0.21
    cpm.omex=1.-cpm.omec
    pkm=power.PowerSpectrum(cpm)

    cpp=copy.copy(p.cp)
    cpp.omec=0.31
    cpp.omex=1.-cpp.omec
    pkp=power.PowerSpectrum(cpp)

    print pk0.D1(1e3-1), pkm.D1(1e3-1), pkp.D1(1e3-1)



    # read data #
    f0=np.load('../../workspace/result/cosdep/lcdm_026/traj_spherical_collapse.npz')
    fm=np.load('../../workspace/result/cosdep/lcdm_021/traj_spherical_collapse.npz')
    fp=np.load('../../workspace/result/cosdep/lcdm_031/traj_spherical_collapse.npz')

    a=f0['a']
    traj0=f0['traj']
    trajm=fm['traj']
    trajp=fp['traj']

    rho0, _rhom, _rhop = 1+traj0[:,:,1], 1+trajm[:,:,1], 1+trajp[:,:,1]
    r0, rm, rp = traj0[:,:,0], trajm[:,:,0], trajp[:,:,0]

    rhom=sitp.interp1d(_rhom[:,0], _rhom[:,-1],kind='linear',fill_value='extrapolate')(rho0[:,0])
    rhop=sitp.interp1d(_rhop[:,0], _rhop[:,-1],kind='linear',fill_value='extrapolate')(rho0[:,0])
    
    # 
    #pl.plot(rho0[:,-1], rhop[:,-1]/rho0[:,-1]-1, 'k-')
    #pl.plot(rho0[:,-1], rhom[:,-1]/rho0[:,-1]-1, 'r-')

    pl.plot(rho0[:,-1], rhop/rho0[:,-1]-1, 'k-', label=r'$\Omega_m=0.31$', lw=2)
    pl.plot(rho0[:,-1], rhom/rho0[:,-1]-1, 'r-', label=r'$\Omega_m=0.21$', lw=2)

    #pl.plot(r0[:,-1], rp[:,-1]/r0[:,-1], 'k-')
    #pl.plot(r0[:,-1], rm[:,-1]/r0[:,-1], 'r-')

    pl.axvline(x=1., color='k', ls='--')
    pl.legend(loc='upper center', fontsize=17)


    pl.xlim([0,5])
    pl.ylim([-0.15, 0.15])

    pl.xlabel(r'$\rho_m(\Omega_{m,fid})$', size=20)
    pl.ylabel(r'$\rho_m(\Omega_{m})/\rho_m(\Omega_{m, fid})-1$', size=20)

    pl.show()


    
    # ->> The End <<- #
    p.finalize()


