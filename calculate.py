import json
import numpy as np

import viewers

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
atoms = [(atom['Z'], np.array(atom['loc'])*l_consts)
         for atom in params['atoms']]
l_vecs = np.array([np.array(vec)*l_consts for vec 
                    in params['primitive transforms']])

#Generate a calculator object for the system
solver = solver_dict[params['solver']]
calculator = solver.Calculator(atoms, l_vecs, **params.get('options',{}))


viewers.view_band_structure(calculator,50)
