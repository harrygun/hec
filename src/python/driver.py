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
	    n_p=int(np.ceil(2.*e/p_step)+2)
            _p=list(np.linspace(-e, e, n_p) )
            _e=[e]*n_p

	elist=elist+_e
	plist=plist+_p

        pnlist.append(n_p)


    print 'plist.shape', np.array(plist).shape, np.array(elist).shape
    print 'pnlist: ', len(pnlist), np.sum(pnlist)


    return  [elist, plist, pnlist]






def get_hec_trajs(p, dynvar, a, var_type):

    if not var_type=='nu_e_p':
        raise Exception('only nu_e_p type of initial condition is supported right now.')
    
    rho_lst, e_lst, p_lst, pn_lst=dynvar
    rho, ell = mar.meshgrid(rho_lst, e_lst)
    rho, pro = mar.meshgrid(rho_lst, p_lst)

    print 'shape:', rho.shape, ell.shape, pro.shape, ', e_lst len:', len(e_lst)


    # ->> mpi idx <<- #
    idx_fulllist=range(len(rho_lst))
    if p.mpi:
        idx=np.array_split(idx_fulllist, mpi.size)[mpi.rank]
    else:
        idx=idx_fulllist

    print 'rank: ', mpi.rank, idx
    quit()


    # ->> now run <<- #
    for i in idx:
	for j in range(len(e_lst)):

	    # ->> convert to eigenvalues first <<- #
            lamb=elc.shape_to_eigval(rho[i,j], ell[i,j], pro[i,j])

            # ->> 
	    traj=[]
	    _traj=np.zeros(e_lst.shape)
	    _traj[j]=elc.get_elliptraj_one(p, a, lamb):

	    traj.append(_traj)





    # ->> gather <<- #
    if mpi.size>1:
        sg_ = np.array(mpi.world.gather(traj, root=0))
    
        if mpi.rank0:
            sg=np.concatenate((sg_[0], sg_[1]), axis=0)

            for i in range(2, mpi.size):
                sg=np.concatenate((sg, sg_[i]), axis=0)

        else:
            sg=None
    
        sg=mpi.world.bcast(sg, root=0)

        return sg


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



    '''->> calculate ellipsoidal collapse model <<-'''
    ai, af, na = 0.01, 1., 500
    a=np.linspace(ai, af, na)


    #->> dynvar:  list of rho and e, p <<- #
    rho_lst=np.linspace(-1., 5., 101)
    #e_lst, p_lst, pn_lst=ep_list([0, 50.], 201, p_step=0.5)
    dynvar=[rho_lst]+ep_list([0, 50.], 201, p_step=0.5)


    ''' ->> now calculate the trajectories <<- '''
    var_type='nu_e_p'    #eigenvalues
    #elc.get_elliptraj(a, dynvar, var_type=var_type)
    get_hec_trajs(p, dynvar, a, var_type)
     





    
    # ->> The End <<- #
    p.finalize()



