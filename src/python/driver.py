import os
import numpy as np
import pylab as pl

import genscript.progcontrol as pc
from genscript.extendclass import *
import genscript.mpiutil as mpi
import genscript.myarray as mar

import elc as elc






def ep_list(erange, elen, p_step=0.5):
    '''erange:   the range of ellipticity 
       elen:     the length of ellipticity list 
    '''
    _el_=np.linspace(erange[0], erange[1], elen)

    elist, plist, pnlist = [], [], []
    for e in _el_:
        # ->> prolateness -e<=p<=e
        if e==0:
	    _e=[e]
            _p=_e
	    n_p=1
	else:
	    n_p=np.ceil(2.*e/p_step)+2
            _p=np.linspace(-e, e, n_p) 
            _e=[e]*n_p

	elist.append(_e)
	plist.append(_p)
        pnlist.append(n_p)

    return  [elist, plist, pnlist]






def get_hec_trajs(p, dynvar, a, var_type):

    #dynvar=[rho_lst, e_lst, p_lst, pn_lst]

    if p.mpi  :







    return








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


    ''' MPI initialization '''

    if p.mpi:
	calidx=mpi.mpirange(len(type_idx))
        print 'rank: ', mpi.rank, calidx, type_idx[calidx]
    else:
        calidx=range(len(type_idx))
        print calidx, type_idx[calidx]





    '''->> calculate ellipsoidal collapse model <<-'''
    ai, af, na = 0.01, 1., 500
    a=np.linspace(ai, af, na)


    #->> dynvar:  list of rho and e, p <<- #
    rho_lst=np.linspace(-1., 5., 101)
    #e_lst, p_lst, pn_lst=ep_list([0, 50.], 201, p_step=0.5)
    dynvar=[rho_list]+ep_list([0, 50.], 201, p_step=0.5)


    ''' ->> now calculate the trajectories <<- '''
    var_type='nu_e_p'    #eigenvalues
    #elc.get_elliptraj(a, dynvar, var_type=var_type)
    get_hec_trajs(p, dynvar, a, var_type)
     





    
    # ->> The End <<- #
    p.finalize()



