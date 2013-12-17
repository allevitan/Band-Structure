#General purpose tools that are useful for manipulating lattices

import numpy as np

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

def fermi_energy(Es, num_states):
    """Returns the fermi energy"""
    return np.sort(Es)[int(num_states/2)-1]

