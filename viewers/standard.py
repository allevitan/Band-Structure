import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from lattice_tools import fermi_energy, bravais_lattice, reciprocal_vectors
import itertools as it

def view_band_structure(calculator, resolution):
    num_bands = sum(zip(*calculator.atoms)[0]) / 2 + 10 # of bands to show
    bands = np.zeros((num_bands,2*resolution))
    
    for i in range(0,resolution):
        bands[:,resolution-1-i] = calculator.at_k((calculator.r_vecs[0] * i) / (2 * (resolution-1)), num_bands=num_bands)
    for i in range(0,resolution):
        bands[:,resolution+i] = calculator.at_k(((calculator.r_vecs[0] + calculator.r_vecs[1]) * i) / (2 * (resolution-1)),num_bands=num_bands)


    plt.plot(bands.transpose() - fermi_energy(np.ravel(bands),
                                    2*resolution*sum(zip(*calculator.atoms)[0])))
    plt.xlabel('Reciprocal Lattice Points')
    plt.ylabel('Energy (Rydbergs)')
    plt.axis([0,2*resolution,-3,3])
    plt.show()


def plot_first_brillouin(l_vecs, ax=None, plot=True):
    """
    plots the first brillouin zone of the lattice vectors 
    defined in r_vecs. If ax is given, it will plot it
    on that axis - if not, it will make it's own.
    """

    lattice = bravais_lattice(l_vecs)
    r_vecs = reciprocal_vectors(l_vecs)

    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        
    
    if lattice == 'face-centered-cubic':
        F = np.array([[0.25,0.5,0.75],
                      [0.5,0.25,0.75],
                      [0.75,0.25,0.5],
                      [0.75,0.5,0.25],
                      [0.5,0.75,0.25],
                      [0.25,0.75,0.5],
                      [0.25,0.5,0.75]])
        for signs in it.product((1,-1),repeat=3):
            points = np.array(signs) * np.sum(r_vecs * F[:,:,np.newaxis],axis=1)
            ax.plot(*zip(*list(points)),color='k')

    elif lattice == 'body-centered-cubic':
        F = np.array([[-0.5,0.5,0.5],
                      [0.25,0.25,0.25],
                      [0.5,0.5,-0.5],
                      [0.25,0.25,0.25],
                      [0.5,-0.5,0.5]]) #kinda hacky, but it works
        for signs in it.product((1,-1),repeat=3):
            points = np.array(signs) * np.sum(r_vecs * F[:,:,np.newaxis],axis=1)
            ax.plot(*zip(*list(points)),color='k')

    elif lattice == 'primitive-cubic' or lattice == 'primitive-tetragonal' \
         or lattice == 'primitive-orthorhombic':
        F = np.array([[-0.5,0.5,0.5],
                      [0.5,0.5,0.5],
                      [0.5,-0.5,0.5],
                      [0.5,0.5,0.5],
                      [0.5,0.5,-0.5]])
        for signs in it.product((1,-1),repeat=3):
            points = np.array(signs) * np.sum(r_vecs * F[:,:,np.newaxis],axis=1)
            ax.plot(*zip(*list(points)),color='k')

    for vec in list(r_vecs):
        ax.plot(*zip([0,0,0],list(vec)), color='b')
    
    if plot:
        plt.show()


