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





def gather_mpi(traj):
    # ->> gather <<- #

    if mpi.size>1:
        print 'gather all results from all threads'
        _trajall = np.array(mpi.world.gather(traj, root=0))
    
        if mpi.rank0:
            trajall=np.concatenate((_trajall[0], _trajall[1]), axis=0)
            for i in range(2, mpi.size):
                trajall=np.concatenate((trajall, _trajall[i]), axis=0)
        else:
            trajall=None
    
        trajall=mpi.world.bcast(trajall, root=0)
        return np.array(trajall)
    else:
        return np.array(traj)



''' --------------------------------------------------------
              ->>   trajectories generator <<- 
    --------------------------------------------------------
'''
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
        traj=[]

        _traj=np.zeros(e_lst.shape)
	for j in range(len(e_lst)):
	    # ->> convert to eigenvalues first <<- #
            lamb=elc.shape_to_eigval(rho[i,j], ell[i,j], pro[i,j])
	    _traj[j]=elc.get_elliptraj_one(p, a, lamb)

        traj.append(_traj)


    # ->> gather <<- #
    if mpi.size>1:
        _trajall = np.array(mpi.world.gather(traj, root=0))
    
        if mpi.rank0:
            trajall=np.concatenate((_trajall[0], _trajall[1]), axis=0)
            for i in range(2, mpi.size):
                trajall=np.concatenate((trajall, _trajall[i]), axis=0)
        else:
            trajall=None
    
        trajall=mpi.world.bcast(trajall, root=0)
        return np.array(trajall)
    else:
        return np.array(traj)





def get_sc_trajs(p, dynvar, a, var_type):
    ''' ->> spherical collapse trajectories <<- '''

    if not var_type=='rho_only':
        raise Exception('only `rho_only` IC is supported here.')
    
    rho=dynvar
    print '`get_sc_trajs`: rho shape:', rho.shape


    # ->> MPI parallization idx <<- #
    idx_fulllist=range(len(rho))
    if p.mpi:
        idx=np.array_split(idx_fulllist, mpi.size)[mpi.rank]
    else:
        idx=idx_fulllist

    #print 'rank: ', mpi.rank, idx


    # ->> now run trajectory integrator ... <<- #
    traj=[]
    for i in idx:

	print '{0}-th SC traj: ic_rho={1}'.format(i, rho[i])

        # ->> convert to eigenvalues first <<- #
        lamb=elc.shape_to_eigval(rho[i], 0., 0.)
        _traj=elc.get_elliptraj_one(p, a, lamb)[:,:3]

	ctraj=elc.lambda_comving(_traj[:,0], _traj[:,1], _traj[:,2], a)
        delta=elc.clambda_to_rho(ctraj[0], ctraj[1], ctraj[2])
	dd=np.rollaxis(np.concatenate((ctraj, np.array([delta])), axis=0), 1)

        traj.append(dd)

    return gather_mpi(traj)








''' --------------------------------------------------------
              ->> generator of trajectories <<- 
    --------------------------------------------------------
'''
traj_type_list=['testing', 
                'spherical_collapse',
		'2D_ellipsoidal_collapse',
                ]

def generate_trajs(p, traj_type, a='default', para_boundary='default'):

    if not traj_type in traj_type_list:
        raise Exception('traj_type NOT supported.')


    '''->> calculate ellipsoidal collapse model <<-'''
    # ->> initialization <<- #
    if (a=='default'):
        ai, af, na = 0.01, 1., 200
        a=np.linspace(ai, af, na)


    # ->> perform some tests <<- #
    if traj_type=='testing':
        elc.elltraj_test(p, a)
        traj=None

    # ->> run spherical collapse <<- #
    if traj_type=='spherical_collapse':

        #->> dynvar:  list of rho and e, p <<- #
        rho_lst=np.linspace(-1., 1.686, 501)
        traj=get_sc_trajs(p, rho_lst, a, 'rho_only')

        print 'final SC traj shape:', traj.shape


    # ->> run ellipsoidal collapse parameter space <<- #
    if traj_type=='2D_ellipsoidal_collapse':
        raise Exception('EC NOT supported yet.')


    return traj




