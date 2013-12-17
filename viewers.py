import matplotlib.pyplot as plt
import numpy as np
from solvers.lattice_tools import fermi_energy

def view_band_structure(calculator, resolution):
    num_bands = sum(zip(*calculator.atoms)[0]) / 2 + 10 # of bands to show
    bands = np.zeros((num_bands,2*resolution))
    
    for i in range(0,resolution):
        bands[:,resolution-1-i] = calculator.at_k((calculator.r_vecs[0] * i) / (2 * (resolution-1)), num_bands=num_bands)
    for i in range(0,resolution):
        bands[:,resolution+i] = calculator.at_k(((calculator.r_vecs[0] + calculator.r_vecs[1]) * i) / (2 * (resolution-1)),num_bands=num_bands)


    plt.plot(bands.transpose() - fermi_energy(np.ravel(bands),
                                    2*resolution*sum(zip(*calculator.atoms)[0])))
    plt.xlabel('Reciprocal Lattice Points')
    plt.ylabel('Energy (Rydbergs)')
    plt.axis([0,2*resolution,-3,3])
    plt.show()
