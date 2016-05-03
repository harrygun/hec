import os
import numpy as np
import pylab as pl

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import genscript.myarray as mar

import elc as elc
import trajs_generator as gtraj







param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'a_init': 1e-2,
    'smooth_R': 0,
    'smooth_type': 'Gauss', 
    'smooth_R_list_type':  'linear'
    }

prog_control={
    #'do_testing': False, 
    #-------------------------------#
    #'do_spherical_collapse':    True,
    #-------------------------------#
    'traj_generator_type':  'testing'
    }





if __name__=='__main__':

    # ->> initialization <<- #
    init_dict=myDict(prog_control)+myDict(param_dict)
    p=pc.prog_init(**init_dict)

    root='../../workspace/result/'
     


    # ------------->> generating trajectories <<--------------- #
    # ->> traj_type_list=['testing', 'spherical_collapse',  <<- #
    # ->>                 '2D_ellipsoidal_collapse', ]      <<- #
    # ------------->> generating trajectories <<--------------- #

    ai, af, na = 0.01, 1., 200
    a=np.linspace(ai, af, na)

    #traj_type='testing'
    traj_type='spherical_collapse'

    traj=gtraj.generate_trajs(p, traj_type, a=a)


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


