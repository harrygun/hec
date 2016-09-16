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
    # ->> dynvar:  grid of dynamical variables <<- #

    if not var_type=='nu_e_p':
        raise Exception('only nu_e_p type of initial condition is supported right now.')
    if len(np.array(dynvar).shape)>2:
        raise Exception()
    
    rho, ell, pro = dynvar
    ndim=len(rho)
    print 'dynvar shape:', rho.shape, ell.shape, pro.shape, 'ndim=', ndim


    # ->> mpi idx <<- #
    idx_fulllist=range(ndim)
    if p.mpi:
        idx=np.array_split(idx_fulllist, mpi.size)[mpi.rank]
    else:
        idx=idx_fulllist

    print 'rank: ', mpi.rank, idx

    #idx=[488]


    # ->> now run <<- #
    traj=[]
    for i in idx:

        # ->> convert to eigenvalues first <<- #
        lamb=elc.shape_to_eigval(rho[i], ell[i], pro[i])
	print '{0}-th HEC traj: ic_rho={1}-{2}-{3}, lamb={4}'.format(i, rho[i], ell[i], pro[i], lamb)

        _traj=elc.get_elliptraj_one(p, a, lamb)[:,:3]
	ctraj=elc.lambda_comving(_traj[:,0], _traj[:,1], _traj[:,2], a)

        # ->> only return comoving lambda and delta <<- #
        delta=elc.clambda_to_rho(ctraj[0], ctraj[1], ctraj[2])
        dd=np.rollaxis(np.concatenate((np.array(ctraj), np.array([delta])), axis=0), 1)

        traj.append(dd)

    quit()


    return gather_mpi(traj)





def get_sc_trajs(p, dynvar, a, var_type):
    ''' ->> spherical collapse trajectories <<- '''
    #print 'rank: ', mpi.rank, idx

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


    # ->> now run trajectory integrator ... <<- #
    traj=[]
    for i in idx:

        # ->> convert to eigenvalues first <<- #
        lamb=elc.shape_to_eigval(rho[i], 0., 0.)

	print '{0}-th SC traj: ic_rho={1}, lamb={2}'.format(i, rho[i], lamb)

        _traj=elc.get_elliptraj_one(p, a, lamb)[:,:3]

	ctraj=elc.lambda_comving(_traj[:,0], _traj[:,1], _traj[:,2], a)
        delta=elc.clambda_to_rho(ctraj[0], ctraj[1], ctraj[2])
	dd=np.rollaxis(np.concatenate((np.array([ctraj[0]]), np.array([delta])), axis=0), 1)


        traj.append(dd)

    return gather_mpi(traj)







''' --------------------------------------------------------
              ->> generator of trajectories <<- 
    --------------------------------------------------------
'''
traj_type_list=['testing', 
                'spherical_collapse',
		'2D_ellipsoidal_collapse',
		'ellipsoidal_collapse_single_ep',
                ]

def generate_trajs(p, traj_type, a='default', para_boundary='default', **pardict):

    if not traj_type in traj_type_list:
        raise Exception('traj_type NOT supported.')


    '''->> calculate ellipsoidal collapse model <<-'''
    # ->> initialization <<- #
    if (a=='default'):
        ai, af, na = 0.001, 1., 200
        a=np.linspace(ai, af, na)


    # ->> perform some tests <<- #
    if traj_type=='testing':
        elc.elltraj_test(p, a)
        traj=None

    # ->> run spherical collapse <<- #
    if traj_type=='spherical_collapse':

        #->> dynvar:  list of rho and e, p <<- #
        #rho_lst=np.linspace(-1., 1.68, 500)
        rho_lst=np.linspace(-10., 1.68, 500)

        traj=get_sc_trajs(p, rho_lst, a, 'rho_only')

        print 'final SC traj shape:', traj.shape


    if traj_type=='ellipsoidal_collapse_single_ep':
        # ->> get some extra information <<- #
        try:
            sig=pardict['sig'], 
            _el, _pl = pardict['ev'][0], pardict['pv'][0]
	    #ext_flag=pardict['ext_flag']
	except:
	    raise Exception()


        # ->>  generating rho-e-p grid data <<- #
	# ->> e is actually e*nu=e*delta/sig <<- #

        rhol=np.linspace(-10., 1.68, 500) # which is actually delta_rho #
        el=_el/(rhol/sig)
        pl=_pl*np.ones(rhol.shape)

        traj=get_hec_trajs(p, [rhol, el, pl], a, 'nu_e_p')

        print 'final HEC traj shape:', traj.shape


    return traj




