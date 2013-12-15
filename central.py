from __future__ import division
import numpy as np
import scipy.linalg as spla #scipy doesn't autoload it's linalg module
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools as it
from time import time

def reciprocal_vectors(l_vecs):
    """
    Finds the reciprocal vectors for a set of 3 arbitrary lattice
    vectors.
    """
    vol = unit_vol(l_vecs)
    r_vecs = np.zeros((3,3))
    for i in range(0,3):
        r_vecs[i,:] = 2 * np.pi * np.cross(l_vecs[i-2,:],l_vecs[i-1,:]) / vol
    return r_vecs

def unit_vol(l_vecs):
    return np.dot(l_vecs[0,:],np.cross(l_vecs[1,:],l_vecs[2,:]))

def make_basis(r_vecs,r):
    """
    Makes an array, each row of which is a vector in reciprocal
    space. Each of these vectors is the wavevector of a complex
    exponential in the basis. r is the radius of the sphere that
    we will constrain our wavevectors to remain within.
    """
    #note: find a better way to calculate n than just making it big
    n=10
    scalings = list(it.product(range(-n,n+1),repeat=3))
    basis = []
    for scaling in scalings:
        vector = sum((r_vecs[i,:] * scaling[i] for i in range(0,3)))
        if np.linalg.norm(vector) < r:
            basis.append(vector)
    return np.array(basis)

def first_brillouin(r_vecs,resolution):
    """
    Makes an array containing a bunch of wavevectors in the first
    Brillouin zone defined by the reciprocal lattice vectors
    """
    neighbors = [sum((r_vecs[i,:] * scaling[i] for i in range(0,3)))
                 for scaling in it.product((-1,0,1),repeat=3)
                 if scaling != (0,0,0)]
    neighbors = np.array(neighbors)
    scale_factor = np.amax(np.power(np.sum(np.power(r_vecs,2),axis=1),0.5))
    possibilities = it.product(range(-resolution,resolution+1),repeat=3)
    possibilities = np.array(tuple(possibilities)) * (scale_factor / 2) / resolution
    dists = np.power(np.sum(np.power(possibilities[:,np.newaxis,:]
                                     - neighbors,2),axis=2),0.5)
    c_dists = np.power(np.sum(np.power(possibilities,2),axis=1),0.5)
    closest = (np.amin(dists,axis=1) - c_dists) > 0
    return possibilities[np.where(closest)[0],:]
    
#The fourier transform of the coulomb potential about an atom, z=1
fourier_coulomb = lambda(H): -4*np.pi / np.sqrt(np.sum(np.power(H,2),axis=2))

#This will remove the infinities along the diagonal, which makes the 
#matrix useable but doesn't change anything important
def clean_diagonal(H):
    H[np.diag_indices(H.shape[0])] = 0
    return H

def cross_terms(atoms,basis, l_vecs):
    """
    Makes a matrix that contains all the cross terms in the
    hamiltonian. These terms do not depend on the wavevector,
    and by adding the free electon energies along the diagonal
    this matrix will become the hamiltonian at any wavevector.
    """

    locs = np.array([loc for (z,loc) in atoms])
    zs = np.array([z for (z,loc) in atoms])

    #fourier_factor is due to the atomic number and positions of atoms in
    #the unit cell, and modifies the fourier transform of the potential
    fourier_factor = lambda(H): np.sum(np.exp(-1j * np.sum(H[:,:,np.newaxis,:] * locs,axis=3)) * zs,axis=2)
    H0 = basis[:,np.newaxis,:] - basis
    with np.errstate(divide='ignore'):
        H0 = (fourier_factor(H0) * clean_diagonal(fourier_coulomb(H0))) 
    return H0 / unit_vol(l_vecs)

def electron_terms(rho,basis,l_vecs):
    """
    Makes the matrix of cross terms due to the potential from some
    approximation of the fourier transform of the electron density
    function, rho.
    """
    He = basis[:,np.newaxis,:] - basis
    by_vec = {tuple(vec[1:]): vec[0] for vec
              in np.concatenate((rho[np.newaxis,:],basis.transpose())).transpose()}
    #calculates the factors due to the electron 
    def rho_factor(H):
        He = np.zeros((H.shape[0],H.shape[1]))
        for i in range(0,H.shape[0]):
            for j in range(0, H.shape[1]):
                He[i,j] = by_vec.get(tuple(H[i,j]),0)
        return He
    
    with np.errstate(divide='ignore'):
        He = rho_factor(He) * clean_diagonal(-1*fourier_coulomb(He))
    return He 

def free_energies(k, basis):
    """
    Makes a diagonal matrix that contains the free electron energies
    for an electron at k + each wavevector in the basis.
    """
    return np.diag(np.sum(np.power(basis + k,2),axis=1))

def find_rhos(phis, basis, l_vecs):
    """
    Given a gigantic array of the fourier series' of the
    eigenfunctions of the hamiltonian (psi), find the convolution
    of psi(-k)^* with psi(k) - the result is a giant array with the fourier
    series' of the electron density functions due to the original psis.
    """
    #Okay this whole function is going to be stupid ugly. I'm going to
    #make an array of the tensor products of the fourier transforms, then
    #I'm going to grab out the components that match each wavevector
    grabby = basis[np.newaxis,:,:] + basis[:,np.newaxis,:]
    #we're subtracting the negative above
    product = np.conj(phis[:,np.newaxis,:]) * phis[:,:,np.newaxis]
    places = [np.where(np.any(grabby - basis[i,:],axis=2)==False)
             for i in range(0,basis.shape[0])]
    convolution = [np.sum(product[:,places[i][0],places[i][1]],axis=1)
                   for i in range(0,basis.shape[0])]
    return np.array(convolution).transpose() / unit_vol(l_vecs) 


def stupid_band_structure(l_vecs, atoms, resolution,r):
    """
    Calculates a widly inaccurate band structure just to give some
    insight into what's going on - electron/electron interactions are
    completely ignored.
    """
    r_vecs = reciprocal_vectors(l_vecs)
    t = time()
    basis =  make_basis(r_vecs,r)
    print "Basis made: ", basis.shape[0], " vectors in ", time() - t, "seconds"
    t = time()
    H0 = cross_terms(atoms,basis,l_vecs)
    print "Cross terms made: ", time() - t, "seconds"
    t = time()
    num_bands = sum(zip(*atoms)[0]) / 2 + 3 # of bands to show
    bands = np.zeros((num_bands,2*resolution))
    
    for i in range(0,resolution):
        H = H0 + free_energies((r_vecs[0] * i) / (2 * (resolution-1)), basis)
        bands[:,resolution-1-i] = np.sort(np.real(np.linalg.eigvalsh(H)))[0:num_bands]
    for i in range(0,resolution):
        H = H0 + free_energies(((r_vecs[0] + r_vecs[1]) * i) / (2 * (resolution-1)), basis)
        bands[:,resolution+i] = np.sort(np.real(np.linalg.eigvalsh(H)))[0:num_bands]
    print "Band structure calculated: ", (time() - t) / (2*resolution), 'seconds per wavevector'
    plt.plot(bands.transpose())
    plt.xlabel('Reciprocal Lattice Points')
    plt.ylabel('Energy (Rydbergs)')
    plt.axis([0,2*resolution,0,5])
    plt.show()


def band_structure(l_vecs, atoms, resolution,r,bri_res=5):
    """
    Calculates the band structure! Eventually...
    """

    r_vecs = reciprocal_vectors(l_vecs)
    basis =  make_basis(r_vecs,r)
    b_zone = first_brillouin(r_vecs,bri_res)
    H0 = cross_terms(atoms,basis,l_vecs)
    He = np.zeros(H0.shape)
    bands_filled = float(sum(zip(*atoms)[0])) / 2
    e_per_band = np.tile(2,np.ceil(bands_filled))
    if bands_filled != int(bands_filled):
        e_per_band[-1] = 1    
    bands_filled = np.ceil(bands_filled)
    bands_shown = bands_filled + 3 # of bands to show
    print "Initial setup complete"
    for i in range(0,10):
        H0e = H0 + He
        H = H0e + free_energies(np.array((0,0,0)), basis)
        (Es,phis) = np.linalg.eigh(H)
        (Es, phis) = (Es[:bands_filled],phis[:bands_filled,:])
        rhos = find_rhos(phis,basis,l_vecs) * e_per_band[:,np.newaxis]
        rho = np.sum(rhos, axis=0)
        He = electron_terms(rho,basis,l_vecs)
        print "Iteration", i+1, "Complete"
        

    bands = np.zeros((bands_shown,2*resolution))
    for i in range(0,resolution):
        H = H0 + He + free_energies((r_vecs[0] * i) / (2 * (resolution-1)), basis)
        bands[:,resolution-1-i] = np.sort(np.real(np.linalg.eigvalsh(H)))[0:bands_shown]
    for i in range(0,resolution):
        H = H0 + He + free_energies(((r_vecs[0] + r_vecs[1]) * i) / (2 * (resolution-1)), basis)
        bands[:,resolution+i] = np.sort(np.real(np.linalg.eigvalsh(H)))[0:bands_shown]

    plt.plot(bands.transpose())
    plt.xlabel('Reciprocal Lattice Points')
    plt.ylabel('Energy (Rydbergs)')
    plt.axis([0,2*resolution,0,5])
    plt.show()
        

if __name__=='__main__':
    carbon_l_const = 6.7 #bohr radii
    carbon_l_vecs = np.array(((0.5,0.5,0),
                       (0,0.5,0.5),
                       (0.5,0,0.5))) * carbon_l_const
    carbon_atoms = ((6,np.array((0,0,0))),(6,np.array((0.25,0.25,0.25))*carbon_l_const))
    silicon_l_const = 10.26 #bohr radii
    silicon_l_vecs = np.array(((0.5,0.5,0),
                       (0,0.5,0.5),
                       (0.5,0,0.5))) * silicon_l_const
    silicon_atoms = ((14,np.array((0,0,0))),(14,np.array((0.25,0.25,0.25))*silicon_l_const))

    def brillouin_zone_demo(l_vecs):
        r_vecs = reciprocal_vectors(l_vecs)
        b_zone = first_brillouin(r_vecs,7)
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax.scatter(b_zone[:,0],b_zone[:,1],b_zone[:,2])
        plt.show()
    
    #brillouin_zone_demo(carbon_l_vecs)

    #stupid_band_structure(silicon_l_vecs, silicon_atoms, 50,3)
    #stupid_band_structure(carbon_l_vecs, carbon_atoms, 50,3)
    band_structure(carbon_l_vecs, carbon_atoms, 50,4)
