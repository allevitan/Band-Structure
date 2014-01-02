from __future__ import division
import numpy as np
import itertools as it

def make_basis(r_vecs,n):
    """
    Makes an array, each row of which is a vector in reciprocal
    space. Each of these vectors is the wavevector of a complex
    exponential in the basis. n is a lower bound on the number of
    vectors in the basis - it enforces a symmetric basis set.
    """

    farthest= int(np.ceil(np.power(n*6/np.pi,1/3)/2))
    scalings = np.array(list(it.product(range(-farthest,farthest+1),repeat=3)))
    vectors = np.sum(scalings[:,:,np.newaxis] * r_vecs,axis=1)
    idx = np.argsort(np.sum(np.power(vectors,2),axis=1))
    while np.linalg.norm(vectors[idx[n]]) - np.linalg.norm(vectors[idx[n-1]]) < np.linalg.norm(vectors[idx[n]]) / 1000:
        n+=1
    return vectors[idx[:n]]


#The fourier transform of the coulomb potential about an atom, z=1
fourier_coulomb = lambda(H): -4*np.pi / np.sqrt(np.sum(np.power(H,2),axis=2))


def clean_diagonal(H):
    """
    Removes the infinities along the diagonal, which makes the 
    matrix useable but doesn't change anything important
    """
    H[np.diag_indices(H.shape[0])] = 0
    return H


def cross_terms(lattice,basis):
    """
    Makes a matrix that contains all the cross terms in the
    hamiltonian. These terms do not depend on the wavevector,
    and by adding the free electon energies along the diagonal
    this matrix will become the hamiltonian at any wavevector.
    """

    H0 = basis[:,np.newaxis,:] - basis
    with np.errstate(divide='ignore'):
        H0 = sum((z*factors for (z,factors) 
                  in lattice.structure_factors(H0).iteritems())) \
             * clean_diagonal(fourier_coulomb(H0))
    return H0 / lattice.unit_vol


def free_energies(k, basis):
    """
    Makes a diagonal matrix that contains the free electron energies
    for an electron at k + each wavevector in the basis.
    """
    return np.diag(np.sum(np.power(basis + k,2),axis=1))


class Calculator():

    def __init__(self, lattice, basis_size=200):
        self.lattice = lattice
        self.basis = make_basis(self.lattice.r_vecs,basis_size)
        self.H0 = cross_terms(self.lattice,self.basis)
        self.filled_bands = sum([z*points.shape[0] for (z,points)
                                 in self.lattice.basis.iteritems()])
    
    def at_k(self,k,num_bands=None):
        num_bands = num_bands or (self.filled_bands * 3 / 4)
        H = self.H0 + free_energies(k, self.basis)
        return np.sort(np.real(np.linalg.eigvalsh(H)))[0:num_bands]

