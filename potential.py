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

    Outputs the potential at thet point in rydbergs.

    All distance units are in multiples of the bohr radius.
    """
    
    displacements = (atoms[1].transpose() - point).transpose()
    distances = np.sum(np.abs(displacements)**2,axis=0)**(1./2)
    return -np.sum(atoms[0] / distances)


def make_atoms(l_consts, l_vecs, atoms, n=4):
    """
    Takes a relatively normal representation of crystal structure and
    generates a small lattice, represented the way that point_potential
    wants it.
    """
    l_consts = np.array(l_consts)
    l_vecs = [np.array(l_vec) * l_consts  for l_vec in l_vecs]
    l_vecs += [-l_vec for l_vec in l_vecs]
    (Zs, cores) = zip(*atoms)
    cores = [np.array(core) * l_consts for core in cores]
    
    #Here I've gotta figure out how to build a list of all unique lattice
    #sites within n - currently has anisotropic behavior.
    origins = [np.array((0,0,0)), np.array((5.15,5.15,0))]
    for i in range(1,n+1):
        for combo in it.product(origins, l_vecs):
            scombo = sum(combo)
            if not any((np.all(scombo == origin) for origin in origins)):
                origins.append(scombo)
    cores = np.array([core + origin for core in cores for origin in origins]).transpose()
    Zs = Zs * (cores.shape[1] / len(Zs))
    return (Zs,cores)



l_consts = (10.3,10.3,10.3)
l_vecs = ((0.5,0.5,0),
          (0,0.5,0.5),
          (0.5,0,0.5))
atoms = ((6,(0,0,0)),(6,(0.25,0.25,0.25)))

all_the_atoms = make_atoms(l_consts,l_vecs,atoms,n=2)

point = np.array((0.1,0.2,0.3))
atoms = (np.array((6,6,6)), 
         np.array(((0,0,0),
                   (5,5,0),
                   (0,5,5))).transpose())
t = time()
print point_potential(point,all_the_atoms)
print time() - t
