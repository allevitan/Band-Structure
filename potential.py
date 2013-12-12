from __future__ import division
import numpy as np
import itertools as it
from time import time

def point_potential(point, atoms):
    """
    Calculates the potential at a point due soley to the ionic cores 
    in the "atoms" list.

    Inputs
    ------------------------------------------------
    point: a vector describing the point of interest
    atoms: a tuple of a vector and an array. The vector contains
    atomic numbers, and the array contains lattice points.

    Outputs the potential at thet point in rydbergs. They are
    not defined with respect to any particular zero, but multiple
    calls on the same set of atoms will return answers defined
    from the same zero.

    All distance units are in multiples of the bohr radius.
    """
    displacements = (atoms[1] - point)
    distances = np.sum(np.abs(displacements)**2,axis=1)**(1./2)
    return -np.sum(atoms[0] / distances)


def make_atoms(l_consts, l_vecs, atoms, n=7, radius=50):
    """
    Takes a relatively normal representation of crystal structure and
    generates a small lattice, represented the way that point_potential
    wants it.
    """
    l_consts = np.array(l_consts)
    l_vecs = [np.array(l_vec) * l_consts  for l_vec in l_vecs]
    (Zs, cores) = zip(*atoms)
    cores = [np.array(core) * l_consts for core in cores]
    coordinates = it.product(range(-n,n+1),range(-n,n+1),range(-n,n+1))
    origins = (l_vecs[0]*coord[0] + l_vecs[1]*coord[1] + l_vecs[2]*coord[2]
               for coord in coordinates)
    #Make the atoms approximately spherical, might help?
    origins = [origin for origin in origins if np.linalg.norm(origin) < radius]
    cores = np.array([core + origin for core in cores for origin in origins])
    Zs = Zs * int(cores.shape[0] / len(Zs))
    return (np.array(Zs),cores)


def all_points(l_consts, l_vecs, atoms, dim, n=4):
    """
    Generates a discrete potential function
    
    """
    t = time()
    atoms = make_atoms(l_consts,l_vecs, atoms, n=n)
    l_consts = np.array(l_consts)
    l_vecs = [np.array(l_vec) * l_consts  for l_vec in l_vecs]
    
    #locations is 5D array, the first 3 dimensions are mirroring positions in
    #the final output array, the 4th dimension is a singleton for broadcasting
    #purposes, and the 5th contains a vector of position at that location
    locations = [[[ [a*l_vecs[0]/dim + b*l_vecs[1]/dim + c*l_vecs[2]/dim]
                    for a in range(0,dim)]
                  for b in range(0,dim)]
                 for c in range(0,dim)]
    locations = np.array(locations)
    displacements = np.zeros((dim,dim,dim) + atoms[1].shape)
    displacements[:,:,:] = atoms[1]
    displacements -= locations
    distances = np.sum(np.abs(displacements)**2,axis=4)**(1./2)
    potential = -np.sum(atoms[0] / distances, axis=3) 
    return potential


l_consts = (10.3,10.3,10.3)
l_vecs = ((0.5,0.5,0),
          (0,0.5,0.5),
          (0.5,0,0.5))
atoms = ((6,(0,0,0)),(6,(0.25,0.25,0.25)))

print all_points(l_consts,l_vecs,atoms,10,n=4)
