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

    #traj_type='testing'
    traj_type='spherical_collapse'
    traj=gtraj.generate_trajs(p, traj_type)


    ''' ->> now save data <<- '''
    if (traj!=None)&(mpi.rank0):
	fname=root+'traj.dat'
	traj.tofile(fname)



    
    # ->> The End <<- #
    p.finalize()



