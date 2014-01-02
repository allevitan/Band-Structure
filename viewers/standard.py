import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import itertools as it

from path_tools import path_defaults

def fermi_energy(Es, num_states):
    """Returns the fermi energy"""
    return np.sort(Es)[int(num_states/2)-1]

def view_band_structure(calculator, resolution):
    num_bands = calculator.filled_bands / 2 + 10 # of bands to show
    bands = np.zeros((num_bands,2*resolution))
    
    for i in range(0,resolution):
        bands[:,resolution-1-i] = calculator.at_k((calculator.lattice.r_vecs[0] * i) / (2 * (resolution-1)), num_bands=num_bands)
    for i in range(0,resolution):
        bands[:,resolution+i] = calculator.at_k(((calculator.lattice.r_vecs[0] + calculator.lattice.r_vecs[1]) * i) / (2 * (resolution-1)),num_bands=num_bands)


    plt.plot(bands.transpose() - fermi_energy(np.ravel(bands),
                                              2*resolution*calculator.filled_bands))
    plt.xlabel('Reciprocal Lattice Points')
    plt.ylabel('Energy (Rydbergs)')
    plt.axis([0,2*resolution,-3,3])
    plt.show()


def plot_first_brillouin(lattice, ax=None, plot=True):
    """
    plots the first brillouin zone of the lattice vectors 
    defined in l_vecs.
    lattice can either be a string defining the bravais lattice,
    or None, in which case the program will try it's best to
    figure out the lattice from the vectors
    If ax is given, it will plot the brillouin zone on that axis
    If not, it will make it's own.
    """
    
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        
    if lattice.bravais == 'face-centered-cubic':
        F = np.array([[0.25,0.5,0.75],
                      [0.5,0.25,0.75],
                      [0.75,0.25,0.5],
                      [0.75,0.5,0.25],
                      [0.5,0.75,0.25],
                      [0.25,0.75,0.5],
                      [0.25,0.5,0.75]])
        for signs in it.product((1,-1),repeat=3):
            points = np.array(signs) * np.sum(lattice.r_vecs
                                              * F[:,:,np.newaxis],axis=1)
            ax.plot(*zip(*list(points)),color='k')

    elif lattice.bravais == 'body-centered-cubic':
        F = np.array([[-0.5,0.5,0.5],
                      [0.25,0.25,0.25],
                      [0.5,0.5,-0.5],
                      [0.25,0.25,0.25],
                      [0.5,-0.5,0.5]]) #kinda hacky, but it works
        for signs in it.product((1,-1),repeat=3):
            points = np.array(signs) * np.sum(lattice.r_vecs 
                                              * F[:,:,np.newaxis],axis=1)
            ax.plot(*zip(*list(points)),color='k')

    elif lattice.bravais == 'primitive-cubic' \
         or lattice.bravais == 'primitive-tetragonal' \
         or lattice.bravais == 'primitive-orthorhombic':
        F = np.array([[-0.5,0.5,0.5],
                      [0.5,0.5,0.5],
                      [0.5,-0.5,0.5],
                      [0.5,0.5,0.5],
                      [0.5,0.5,-0.5]])
        for signs in it.product((1,-1),repeat=3):
            points = np.array(signs) * np.sum(lattice.r_vecs
                                              * F[:,:,np.newaxis],axis=1)
            ax.plot(*zip(*list(points)),color='k')

    for vec in list(lattice.r_vecs):
        ax.plot(*zip([0,0,0],list(vec)), color='b')
    
    if plot:
        plt.show()


def plot_path(lattice, path=None, ax=None, plot=True):
    """
    Plots a path in reciprocal space. path can be:
    None: the program will use a default path
    List of strings: will plot a path through stored high-symmetry points.
        'Gamma' is for the origin, other points just use the letter.
    List of lists: will plot a path through the points defined by the
        lists, where [x,y,z] defines the point x*b_1+y*b_2+z*b_3
    """
    
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        
    default_path_info = path_defaults[lattice.bravais]
    
    if path == None:
        path = default_path_info['path']
    
    for (i,point) in enumerate(path):
        if type(point) == type(''):
            if point == 'Gamma':
                path[i] = [0,0,0]
            else:
                path[i] = default_path_info['points'][point]

    path = np.array(path)
    path = np.sum(lattice.r_vecs * path[:,:,np.newaxis],axis=1)
    
    ax.plot(*zip(*path),color='r',marker='.')
    ax.annotate('sup',xy=(1,1),xytext=(-20,20))
    if plot:
        plt.show()
    

