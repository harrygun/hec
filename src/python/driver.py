import os
import numpy as np
import pylab as pl

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi

import elc as elc





param_dict={
    'power_spectrum_fname': '/home/xwang/workspace/general-data/power/fiducial_matterpower.dat',
    'a_init': 1e-2,
    'smooth_R': 0,
    'smooth_type': 'Gauss', 
    'smooth_R_list_type':  'linear'
    }

prog_control={
    'do_testing': False, 
    #-------------------------------#
    #-------------------------------#
    }





if __name__=='__main__':

    # ->> initialization <<- #
    init_dict=myDict(prog_control)+myDict(param_dict)
    p=pc.prog_init(**init_dict)

    root='../../workspace/result/'






    '''->> calculate ellipsoidal collapse model <<-'''
    ai, af, na = 0.01, 1., 500
    a=np.linspace(ai, af, na)


    #->> dynvar:  list of rho and e, p <<- #
    rhor, er, pr = [], [], []
    dynvar=[rhor, er, pr]

    var_type='nu_e_p'    #eigenvalues
    elc.get_elliptraj(a, dynvar, var_type=var_type)

     





    
    # ->> The End <<- #
    p.finalize()



