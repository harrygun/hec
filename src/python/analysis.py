import numpy as np
import pylab as pl

import genscript.cyth.cubicspline as intp











if __name__=='__main__':

    folder='../../workspace/result/'
    fn_1=folder+'traj_spherical_collapse.npz'
    fn_2=folder+'old/traj_spherical_collapse.npz'


    f1, f2=np.load(fn_1), np.load(fn_2)
    a1, traj1, a2, traj2=f1['a'], f1['traj'], f2['a'], f2['traj']

    print traj1.shape, a1[0], a2[0]

    #quit()


    n=50
    dlist=np.linspace(-1, 1.6, n)

    for i in range(n): 

        pl.plot(a1, traj1[int(i*500./float(n)),:,1], 'k-')
        pl.plot(a2, traj2[int(i*500./float(n)),:,1], 'r:')

        #pl.plot(a2, traj1[int(i*500./float(n)),:,1])


    pl.show()
    

