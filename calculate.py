import json
import numpy as np

from viewers import standard as viewer
#import bravais
from lattice import Lattice

import solvers.central
solver_dict = {"central":solvers.central}

#For now, just edit this line. Elegance will come later.
filename = 'Example_Inputs/Diamond.json'

with open(filename) as param_file:
    params = json.load(param_file)

#Format the input in a standard way
l_consts = (params['lattice constants']['a'],
            params['lattice constants']['b'],
            params['lattice constants']['c'])
basis = {int(z): np.array(points)*l_consts
         for (z,points) in params['atoms'].iteritems()}

l_vecs = np.array([np.array(vec)*l_consts for vec 
                    in params['primitive transforms']])

#Generate a calculator object for the system
lattice = Lattice(l_vecs,basis)
solver = solver_dict[params['solver']]
calculator = solver.Calculator(lattice, **params.get('options',{}))


viewer.view_band_structure(calculator,50)
