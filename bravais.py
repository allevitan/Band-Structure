from __future__ import division
import numpy as np


def bravais_lattice(l_vecs):
    """
    finds the bravais lattice of the lattice vectors
    """
    if l_vecs.shape == (2,2):
        lattice = "oblique"
        dists = [np.linalg.norm(l_vecs[0]),np.linalg.norm(l_vecs[1])]
        cos = np.dot(l_vecs[0],l_vecs[1]) / (dists[0]*dists[1])
        if np.allclose(dists[0],dists[1]):
            if cos==0:
                lattice = "square"
            elif np.allclose(np.abs(cos),0.5):
                lattice = "hexagonal"
        else:
            if cos==0:
                lattice = "rectangular"
            elif np.allclose(dists[0]*cos*2,dists[1]) \
                 or np.allclose(dists[1]*cos*2,dists[0]):
                lattice = "centered-rectangular"                

    if l_vecs.shape == (3,3):
        lattice = "triclinic"
        dists = [np.linalg.norm(l_vecs[0]),np.linalg.norm(l_vecs[1]),
                 np.linalg.norm(l_vecs[1])]
        dorder = np.argsort(dists)
        cosines = [np.dot(l_vecs[0],l_vecs[1]) / (dists[0]*dists[1]),
                   np.dot(l_vecs[1],l_vecs[1]) / (dists[0]*dists[2]),
                   np.dot(l_vecs[2],l_vecs[1]) / (dists[0]*dists[0])]
        corder = np.argsort(cosines)

        if np.allclose(np.array([0,0,0]),cosines):
            if np.allclose(dists[dorder[0]],dists[dorder[2]]):
                lattice = "primitive-cubic"
            elif np.allclose(dists[0],dists[1]) \
                 or np.allclose(dists[1],dists[2]):
                lattice = "primitive-tetragonal"
            else:
                lattice = "primitive-orthorhombic"

        elif np.allclose(np.array([0,0]),cosines[dorder[:2]]):
            if np.allclose(dists[dorder[0]],dists[dorder[2]]):
            else:
                lattice = "primitive-monoclinic"
        
    return lattice



if __name__ == "__main__":
    
    #testing the 2D version
    square_vecs = np.array([[1,0],[0,1]])
    hex_vecs = np.array([[1,0],[0.5,np.sqrt(3)/2]])
    rec_vecs = np.array([[1,0],[0,0.6]])
    cen_rec_vecs = np.array([[1,0],[0.5,0.6]])
    obl_vecs = np.array([[0.3,0.5],[-0.1,0.7]])
    print bravais_lattice(square_vecs)
    print bravais_lattice(hex_vecs)
    print bravais_lattice(rec_vecs)
    print bravais_lattice(cen_rec_vecs)
    print bravais_lattice(obl_vecs)
