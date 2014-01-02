path_defaults = {
    'face-centered-cubic': {
        'points': {
            'L': [0.5,0.5,0.5],
            'K': [0.75,0.375,0.375],
            'W': [0.75,0.5,0.25],
            'X': [0.5,0.5,0],
            'U': [0.625,0.625,0.25]
        },
        'path': ['Gamma','X','W','K','Gamma','L','U','W','L','K']
    },
    'body-centered-cubic': {
        'points': {
            'P': [0.25,0.25,0.25],
            'N': [0.5,0,0],
            'H': [0.5,0.5,-0.5]
        },
        'path':['Gamma','H','N','Gamma','P','H']
    },
    'primitive-cubic': {
        'points':{
            'R':[0.5,0.5,0.5],
            'X':[0,0.5,0],
            'M':[0.5,0.5,0]
        },
        'path': ['Gamma','X','M','Gamma','R','X']
    },
    'primitive-tetragonal': {
        'points': {
        },
        'path': []
    },
    'primitive-orthorhombic': {
        'points': {
        },
        'path': []
    }
}
