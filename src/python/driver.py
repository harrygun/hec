import os
import numpy as np
import pylab as pl

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import genscript.myarray as mar
import cosmology.cosparameter as cospar

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

    # ------------->> generating trajectories <<--------------- #
    # ->> traj_type_list=['testing', 'spherical_collapse',  <<- #
    # ->>                 '2D_ellipsoidal_collapse', ]      <<- #
    # ------------->> generating trajectories <<--------------- #

    ai, af, na = p.a_init, 1., 200
    a=np.linspace(ai, af, na)

    #traj_type='testing'
    #traj_type='spherical_collapse'
    traj_type='ellipsoidal_collapse_single_ep'

    # ->> some further info <<- #
    pardict={}

    if traj_type=='ellipsoidal_collapse_single_ep':
        sig=np.sqrt(misc.sig2(p, 0., R=p.smooth_R, window_type=p.smooth_type) )
        #eg, pg = mar.meshgrid(p.e_list, p.p_list)
	print 'doing ellipsoidal_collapse_single_ep, sig=', sig, 'e/p=', p.e_list, p.p_list

	pardict={'sig': sig, 'e': p.e_list, 'p': p.p_list, }


    # ->> now run the trajectories generator <<- #
    traj=gtraj.generate_trajs(p, traj_type, a=a, **pardict)


    ''' ->> now save data <<- '''
    if ((type(traj)==list)|(type(traj)==np.ndarray))&(mpi.rank0):
	fname_dat=root+'traj_'+traj_type+'.dat'
	fname_dat_a=root+'traj_'+traj_type+'_a.dat'

	fname_npz=root+'traj_'+traj_type
	fname_txt=root+'traj_'+traj_type+'.txt'

	print 'traj shape:', traj.shape

	#->> bindary data <<- #
	traj.tofile(fname_dat)
	a.tofile(fname_dat_a)

	#->> numpy arrays <<- #
	np.savez(fname_npz, traj=traj, a=a)


    
    # ->> The End <<- #
    p.finalize()


