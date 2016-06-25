import numpy as np
import pylab as pl












if __name__=='__main__':

    folder='../../workspace/result/'
    fn_1=folder+'traj_spherical_collapse.npz'
    fn_2=folder+'old/traj_spherical_collapse.npz'


    f1, f2=np.load(fn_1), np.load(fn_2)

    a1, traj1, a2, traj2=f1['a'], f1['traj'], f2['a'], f2['traj']

    print traj1.shape


