import numpy as np
import pylab as pl

import genscript.cyth.cubicspline as intp











if __name__=='__main__':

    folder='../../workspace/result/'
    fn_1=folder+'1/traj_spherical_collapse.npz'
    fn_2=folder+'old/traj_spherical_collapse.npz'

    fn_cp=folder+'newcp/traj_spherical_collapse.npz'
    fn_cp_full=folder+'traj_spherical_collapse.npz'


    f1, f2=np.load(fn_1), np.load(fn_2)
    a1, traj1, a2, traj2=f1['a'], f1['traj'], f2['a'], f2['traj']


    fcp=np.load(fn_cp)
    a_cp, traj_cp=fcp['a'], fcp['traj']

    fcp_full=np.load(fn_cp_full)
    a_cp_full, traj_cp_full=fcp_full['a'], fcp_full['traj']
    print a_cp_full.shape, traj_cp_full.shape



    n=50
    dlist=np.linspace(-1, 1.6, n)

    for i in range(n): 
        #pl.plot(a2, traj2[int(i*500./float(n)),:,1], 'k-')

        pl.plot(a1, traj1[int(i*500./float(n)),:,1], 'k-')
        pl.plot(a_cp, traj_cp[int(i*500./float(n)),:,1], 'b--')


    n=100
    for i in range(n): 
        pl.plot(a_cp_full, traj_cp_full[int(i*500./float(n)),:,1], 'r:')


    pl.ylim([-1, 10]) 
    pl.show()
    

