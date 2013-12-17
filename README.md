Band-Structure is the skeleton of a band structure calculating tool.

It currently implements a central equation solver using a bare coulomb potential - not the quickest or most accurate technique in the book by any means. The program is designed to allow new solvers to be quickly developed and integrated into the program.

The program is used by preparing a json file that describes the lattice and specifies what solver to use. Then, edit calculate.py to point to your input file and run "$python calculate.py."

Band Structure!