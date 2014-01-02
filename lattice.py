import numpy as np

class Lattice:
    
    def __init__(self, l_vecs, basis):
        """
        l_vecs is a numpy array of the lattice vectors
        basis is a dictionary, with atomic numbers mapping to
        numpy arrays storing the location of all atoms with that
        atomic number inside the unit cell
        """
        self.l_vecs = l_vecs
        self.r_vecs = reciprocal_vectors(l_vecs)
        self.unit_vol = unit_vol(l_vecs)
        self.basis = basis
        self.bravais_lattice = bravais_lattice(l_vecs)
        self.space_group = space_group(l_vecs,basis)

    def in_lattice(self,vec):
        """
        Returns whether or not the vector is in the lattice
        """
        coeffs = np.abs(np.dot(np.transpose(self.r_vecs)/(2*np.pi),vec))
        return np.allclose(coeffs,np.floor(coeffs+0.5))
        
    def in_reciprocal_lattice(self,vec):
        """
        Returns whether or not the vector is in the reciprocal lattice
        """
        coeffs = np.abs(np.dot(np.transpose(self.l_vecs)*2*np.pi,vec))
        return np.allclose(coeffs,np.floor(coeffs+0.5))
    
    def structure_factors(self,ks):
        """
        Returns the structure factors for each set of equivalent
        points in basis at the wavevectors in ks.
        For speed, it doesn't check whether the ks are actually in the
        reciprocal lattice, and will give innacurate (i.e. nonzero)
        results when k is not in the reciprocal lattice
        """
        ks = np.expand_dims(ks,axis=len(ks.shape)-1)
        return {z: np.sum(np.exp(-1j*np.sum(ks * points,axis=len(ks.shape)-1)),
                          axis=len(ks.shape)-2)
                for (z,points) in self.basis.iteritems()}

def reciprocal_vectors(l_vecs):
    """
    Finds the reciprocal vectors for a set of 3 arbitrary lattice
    vectors.
    """
    #let's explain what this is doing: the reciprocal lattice vectors
    #are the row vectors of the inverse matrix of the lattice vectors
    #(modified by a factor of 2pi). We transpose the row-major array
    #of lattice vectors to make it the appropriate matrix, and leave
    #the output because it's by default stored in row-major order.
    return np.array(np.linalg.inv(l_vecs.transpose()) * 2 * np.pi)

def unit_vol(l_vecs):
    """
    returns the volume of a unit cell (or area, for a 2d lattice)
    """
    return np.abs(np.linalg.det(l_vecs)) #det(A) = det(A^T)    

def bravais_lattice(l_vecs):
    """
    returns the bravais lattice defined by the lattice vectors
    """
    pass

def space_group(l_vecs,basis):
    """
    returns the space group defined by the lattice and basis
    """
    pass


if __name__ == '__main__':
    l_vecs = np.array([[0.5,0.5,0],
                       [0.5,0,0.5],
                       [0,0.5,0.5]])
    basis = {6:np.array([[0,0,0],
                         [0.25,0.25,0.25]])}

    lattice = Lattice(l_vecs,basis)

