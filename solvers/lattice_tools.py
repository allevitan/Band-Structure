#General purpose tools that are useful for manipulating lattices

import numpy as np

def bravais_lattice(l_vecs):
    """
    finds the bravais lattice of the lattice vectors
    """

    #just fake it for now

    #lattice = "primitive-cubic"
    #lattice = "body-centered-cubic"
    #lattice = "face-centered-cubic" 
    lattice = "primitive-tetragonal"

    return lattice

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
    r_vecs = np.array(np.linalg.inv(l_vecs.transpose()) * 2 * np.pi)
    return r_vecs

def unit_vol(l_vecs):
    """
    returns the volume of a unit cell (or area, for a 2d lattice)
    """
    return np.linalg.det(l_vecs) #det(A) = det(A^T)    

def fermi_energy(Es, num_states):
    """Returns the fermi energy"""
    return np.sort(Es)[int(num_states/2)-1]

