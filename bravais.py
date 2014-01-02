from __future__ import division
import numpy as np
from solvers.lattice_tools import reciprocal_vectors
from viewers import path_tools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


fcc_arr = np.array([[0,1,1],
                    [1,0,1],
                    [1,1,0]])

bcc_arr = np.array([[-1,1,1],
                    [1,-1,1],
                    [1,1,-1]])

pc_arr = np.array([[1,0,0],
                   [0,1,0],
                   [0,0,1]])

pt1_arr = np.array([[1,0,0],
                    [0,1,0],
                    [0,0,0.5]])

pt2_arr = np.array([[1,0,0],
                    [0,1,0],
                    [0,0,3]])

if __name__ == "__main__":
    
    from viewers import standard as viewer
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    viewer.plot_first_brillouin(bcc_arr,lattice='body-centered-cubic',
                                ax=ax,plot=False)
    viewer.plot_path(bcc_arr,lattice='body-centered-cubic',ax=ax,plot=False)
    plt.show()
    #testing the 2D version
    # square_vecs = np.array([[1,0],[0,1]])
    # hex_vecs = np.array([[1,0],[0.5,np.sqrt(3)/2]])
    # rec_vecs = np.array([[1,0],[0,0.6]])
    # cen_rec_vecs = np.array([[1,0],[0.5,0.6]])
    # obl_vecs = np.array([[0.3,0.5],[-0.1,0.7]])
    # print bravais_lattice(square_vecs)
    # print bravais_lattice(hex_vecs)
    # print bravais_lattice(rec_vecs)
    # print bravais_lattice(cen_rec_vecs)
    # print bravais_lattice(obl_vecs)
    
