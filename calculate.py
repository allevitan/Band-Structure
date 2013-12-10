import json

#For now, just edit this line. Elegance will come later.
filename = 'silicon.json'


with open(filename) as param_file:
    params = json.load(param_file)

