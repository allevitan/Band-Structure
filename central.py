from __future__ import division
import numpy as np
import itertools as it

def reciprocal_vectors(l_vecs):
    """
    Finds the reciprocal vectors for a set of 3 arbitrary lattice
    vectors.
    """
    vol = np.dot(l_vecs[0,:],np.cross(l_vecs[1,:],l_vecs[2,:]))
    r_vecs = np.zeros((3,3))
    for i in range(0,3):
        r_vecs[i,:] = 2 * np.pi * np.cross(l_vecs[i-2,:],l_vecs[i-1,:])
    return r_vecs

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

def fourier_coulomb(k):
    """
    Calculates the fourier coefficient of the coulomb potential
    around an atom with atomic number 1 at wavevector k
    """
    return 4*np.pi/(np.linalg.norm(k)**2)

def cross_terms(atoms,basis):
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
    fourier_coulomb = lambda(H): 4*np.pi / np.sqrt(np.sum(np.power(H,2),axis=2))
    def clean_diagonal(H):
        H[np.diag_indices(H.shape[0])] = 0
        return H
    H = basis[:,np.newaxis,:] - basis
    H = fourier_factor(H) * clean_diagonal(fourier_coulomb(H))
    return H

def free_energies(k, basis):
    """
    Makes a diagonal matrix that contains the free electron energies
    for an electron at k + each wavevector in the basis.
    """
    return np.diag(np.sum(np.power(basis + k,2),axis=1))


if __name__=='__main__':
    l_vecs = np.array(((0.5,0.5,0),
                       (0,0.5,0.5),
                       (0.5,0,0.5)))
    atoms = ((6,np.array((0,0,0))),(6,np.array((0.25,0.25,0.25))))
    r_vecs = reciprocal_vectors(l_vecs)
    basis =  make_basis(r_vecs,10)
    H0 = cross_terms(atoms,basis)
