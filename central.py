from __future__ import division
import numpy as np
import scipy.linalg as spla #scipy doesn't autoload it's linalg module
import matplotlib.pyplot as plt
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
    fourier_factor = lambda(H): np.sum(np.exp(1j * np.sum(H[:,:,np.newaxis,:] * locs,axis=3)) * zs,axis=2)
    fourier_coulomb = lambda(H): -4*np.pi / np.sqrt(np.sum(np.power(H,2),axis=2))
    #This will remove the infinities along the diagonal, which makes the 
    #matrix useable but doesn't change anything important
    def clean_diagonal(H):
        H[np.diag_indices(H.shape[0])] = 0
        return H
    H = basis[:,np.newaxis,:] - basis
    with np.errstate(divide='ignore'):
        H = (fourier_factor(H) * clean_diagonal(fourier_coulomb(H))) 
    return H / unit_vol(l_vecs)

def free_energies(k, basis):
    """
    Makes a diagonal matrix that contains the free electron energies
    for an electron at k + each wavevector in the basis.
    """
    return np.diag(np.sum(np.power(basis + k,2),axis=1))

def convolve_phis(phis, basis):
    """
    Given a gigantic array of the fourier series' of the
    eigenfunctions of the hamiltonian (psi), find the convolution
    of psi^* with psi - the result is a giant array with the fourier
    series' of the electron density functions due to the original psis.
    """
    #Okay this whole function is going to be stupid ugly. I'm going to
    #make an array of the tensor products of the fourier transforms, then
    #I'm going to grab out the components that match each wavevector
    grabby = basis[np.newaxis,:,:] - basis[:,np.newaxis,:]
    product = np.conj(phis[:,np.newaxis,:]) * phis[:,:,np.newaxis]
    places = [np.where(np.any(grabby - basis[i,:],axis=2)==False)
             for i in range(0,basis.shape[0])]
    print places[0]
    
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


def band_structure(l_vecs, atoms, resolution,r):
    """
    Calculates the band structure! Eventually...
    """

    r_vecs = reciprocal_vectors(l_vecs)
    basis =  make_basis(r_vecs,r)
    H0 = cross_terms(atoms,basis,l_vecs)
    bands_filled = sum(zip(*atoms)[0]) / 2
    bands_shown = bands_filled + 3 # of bands to show
    print "Initial setup complete"
    
    H = H0 + free_energies(np.array((0,0,0)), basis)
    (Es,phis) = np.linalg.eigh(H)
    (Es, phis) = (Es[:bands_filled],phis[:bands_filled,:])
    convolve_phis(phis,basis)



if __name__=='__main__':
    l_const = 6.7 #bohr radii
    l_vecs = np.array(((0.5,0.5,0),
                       (0,0.5,0.5),
                       (0.5,0,0.5))) * l_const
    atoms = ((6,np.array((0,0,0))),(6,np.array((0.25,0.25,0.25))*l_const))
    #stupid_band_structure(l_vecs, atoms, 20,5)
    band_structure(l_vecs, atoms, 20,5)
